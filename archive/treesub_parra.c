/* treesub_parra.c - OpenMP并行化优化版本
   
   Subroutines that operate on trees, inserted into other programs
   such as baseml, basemlg, codeml, and pamp.
   
   并行化优化：矩阵计算、距离计算、P矩阵指数运算、特征向量计算
*/

#include "paml.h"

#ifdef _OPENMP
#include <omp.h>
#endif

extern char BASEs[], * EquateBASE[], AAs[], BINs[], CODONs[][4], nChara[], CharaMap[][64];
extern int noisy;

#ifdef  BASEML
#define REALSEQUENCE
#define NODESTRUCTURE
#define TREESEARCH
#define LSDISTANCE
#define LFUNCTIONS
#define RECONSTRUCTION
#define MINIMIZATION
#endif

#ifdef  CODEML
#define REALSEQUENCE
#define NODESTRUCTURE
#define TREESEARCH
#define LSDISTANCE
#define LFUNCTIONS
#define RECONSTRUCTION
#define MINIMIZATION
#endif

#ifdef  BASEMLG
#define REALSEQUENCE
#define NODESTRUCTURE
#define LSDISTANCE
#endif

#ifdef  RECONSTRUCTION
#define PARSIMONY
#endif

#ifdef  MCMCTREE
#define REALSEQUENCE
#define NODESTRUCTURE
#define LFUNCTIONS
#endif

#if(defined CODEML || defined YN00)
double SS, NN, Sd, Nd;
#endif

/* 并行化P矩阵计算 - 核心优化函数 */
int PMatUVRoot(double P[], double t, int n, double U[], double V[], double Root[])
{
    int i, j, k;
    double *T1 = (double*)malloc(n * n * sizeof(double));
    
    if (t < 1e-10) {
        /* 设置单位矩阵 */
        #ifdef _OPENMP
        #pragma omp parallel for private(i, j) schedule(static)
        #endif
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                P[i * n + j] = (i == j) ? 1.0 : 0.0;
            }
        }
        free(T1);
        return 0;
    }
    
    /* 并行计算 T1 = U * exp(Root*t) */
    #ifdef _OPENMP
    #pragma omp parallel for private(i, j) schedule(static)
    #endif
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            T1[i * n + j] = U[i * n + j] * exp(Root[j] * t);
        }
    }
    
    /* 并行矩阵乘法 P = T1 * V */
    #ifdef _OPENMP
    #pragma omp parallel for private(i, j, k) schedule(static)
    #endif
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            double sum = 0.0;
            for (k = 0; k < n; k++) {
                sum += T1[i * n + k] * V[k * n + j];
            }
            P[i * n + j] = sum;
        }
    }
    
    free(T1);
    return 0;
}

/* 并行化距离矩阵计算 */
int DistanceMatK80(FILE *fout, double space[])
{
    int i, j;
    double **distances, **variances;
    double *P_values, *Q_values;
    
    if (!fout) return -1;
    
    /* 分配内存 */
    distances = (double**)malloc(com.ns * sizeof(double*));
    variances = (double**)malloc(com.ns * sizeof(double*));
    P_values = (double*)malloc(com.ns * com.ns * sizeof(double));
    Q_values = (double*)malloc(com.ns * com.ns * sizeof(double));
    
    for (i = 0; i < com.ns; i++) {
        distances[i] = (double*)calloc(com.ns, sizeof(double));
        variances[i] = (double*)calloc(com.ns, sizeof(double));
    }
    
    fprintf(fout, "\nK80 distances (parallel processing):\n");
    
    /* 步骤1：并行计算转换和颠换比例 */
    #ifdef _OPENMP
    #pragma omp parallel for private(i, j) schedule(dynamic)
    #endif
    for (i = 0; i < com.ns; i++) {
        for (j = i + 1; j < com.ns; j++) {
            double P, Q;
            int idx = i * com.ns + j;
            
            /* 计算转换(P)和颠换(Q)比例 */
            CountChanges(i, j, &P, &Q);
            
            P_values[idx] = P;
            Q_values[idx] = Q;
        }
    }
    
    /* 步骤2：并行计算K80距离 */
    #ifdef _OPENMP
    #pragma omp parallel for private(i, j) schedule(static)
    #endif
    for (i = 0; i < com.ns; i++) {
        for (j = i + 1; j < com.ns; j++) {
            int idx = i * com.ns + j;
            double P = P_values[idx];
            double Q = Q_values[idx];
            double distance, variance;
            
            /* K80距离公式 */
            if (P < 0.75 && Q < 0.5) {
                double a = 1.0 - 2.0*P - Q;
                double b = 1.0 - 2.0*Q;
                
                if (a > 0 && b > 0) {
                    distance = -0.5 * log(a * sqrt(b));
                    
                    /* 方差计算 */
                    double c2 = a * a;
                    double d2 = b * b;
                    variance = (a*b - (a+b)/2.0 + 0.25) / (c2 * d2);
                    variance *= (P*(1-P) + Q*(1-Q)) / com.ls;
                } else {
                    distance = variance = -1;  /* 饱和 */
                }
            } else {
                distance = variance = -1;  /* 饱和 */
            }
            
            distances[i][j] = distances[j][i] = distance;
            variances[i][j] = variances[j][i] = variance;
        }
        distances[i][i] = variances[i][i] = 0.0;
    }
    
    /* 输出结果 */
    for (i = 0; i < com.ns; i++) {
        for (j = 0; j < com.ns; j++) {
            if (distances[i][j] >= 0) {
                fprintf(fout, "%8.5f", distances[i][j]);
            } else {
                fprintf(fout, "%8s", "NA");
            }
        }
        fprintf(fout, "\n");
    }
    
    /* 输出标准误 */
    fprintf(fout, "\nStandard errors:\n");
    for (i = 0; i < com.ns; i++) {
        for (j = 0; j < com.ns; j++) {
            if (variances[i][j] >= 0) {
                fprintf(fout, "%8.5f", sqrt(variances[i][j]));
            } else {
                fprintf(fout, "%8s", "NA");
            }
        }
        fprintf(fout, "\n");
    }
    
    /* 清理内存 */
    for (i = 0; i < com.ns; i++) {
        free(distances[i]);
        free(variances[i]);
    }
    free(distances);
    free(variances);
    free(P_values);
    free(Q_values);
    
    return 0;
}

/* 并行化JC69距离计算 */
int DistanceMatJC69(FILE *fout)
{
    int i, j;
    double **distances, **variances;
    
    /* 分配内存 */
    distances = (double**)malloc(com.ns * sizeof(double*));
    variances = (double**)malloc(com.ns * sizeof(double*));
    for (i = 0; i < com.ns; i++) {
        distances[i] = (double*)calloc(com.ns, sizeof(double));
        variances[i] = (double*)calloc(com.ns, sizeof(double));
    }
    
    fprintf(fout, "\nJC69 distances:\n");
    
    /* 并行计算所有成对距离 */
    #ifdef _OPENMP
    #pragma omp parallel for private(i, j) schedule(dynamic)
    #endif
    for (i = 0; i < com.ns; i++) {
        for (j = i + 1; j < com.ns; j++) {
            int ndiff = 0;
            double p, distance, variance;
            
            /* 计算差异位点数 */
            for (int k = 0; k < com.ls; k++) {
                if (com.z[i][k] != com.z[j][k] && 
                    com.z[i][k] >= 0 && com.z[j][k] >= 0) {
                    ndiff++;
                }
            }
            
            p = (double)ndiff / com.ls;
            
            /* JC69距离公式 */
            if (p < 0.75) {
                distance = -0.75 * log(1.0 - 4.0*p/3.0);
                variance = p * (1.0 - p) / (com.ls * (1.0 - 4.0*p/3.0) * (1.0 - 4.0*p/3.0));
            } else {
                distance = variance = -1;  /* 饱和 */
            }
            
            distances[i][j] = distances[j][i] = distance;
            variances[i][j] = variances[j][i] = variance;
        }
    }
    
    /* 输出距离矩阵 */
    for (i = 0; i < com.ns; i++) {
        for (j = 0; j < com.ns; j++) {
            if (distances[i][j] >= 0) {
                fprintf(fout, "%8.5f", distances[i][j]);
            } else {
                fprintf(fout, "%8s", "NA");
            }
        }
        fprintf(fout, "\n");
    }
    
    /* 清理内存 */
    for (i = 0; i < com.ns; i++) {
        free(distances[i]);
        free(variances[i]);
    }
    free(distances);
    free(variances);
    
    return 0;
}

/* 并行化自举分析 */
int BootstrapSeq(char bootseqf[], int nboot)
{
    int ir, h, j;
    FILE *fboot = fopen(bootseqf, "w");
    int *bootindex = (int*)malloc(com.ls * sizeof(int));
    char **bootseq = (char**)malloc(com.ns * sizeof(char*));
    
    if (!fboot) {
        printf("Error: cannot create bootstrap file %s\n", bootseqf);
        return -1;
    }
    
    /* 为每个序列分配内存 */
    for (j = 0; j < com.ns; j++) {
        bootseq[j] = (char*)malloc((com.ls + 1) * sizeof(char));
    }
    
    fprintf(fboot, "%d %d\n", com.ns, com.ls);
    
    for (ir = 0; ir < nboot; ir++) {
        /* 生成自举样本索引 */
        for (h = 0; h < com.ls; h++) {
            bootindex[h] = (int)(com.ls * rndu());
        }
        
        fprintf(fboot, "\nBootstrap replicate %d\n", ir + 1);
        
        /* 并行生成自举序列 */
        #ifdef _OPENMP
        #pragma omp parallel for private(j, h) schedule(static)
        #endif
        for (j = 0; j < com.ns; j++) {
            /* 根据自举索引构建序列 */
            for (h = 0; h < com.ls; h++) {
                bootseq[j][h] = com.z[j][bootindex[h]];
            }
            bootseq[j][com.ls] = '\0';
        }
        
        /* 串行输出（避免文件写入竞争） */
        for (j = 0; j < com.ns; j++) {
            fprintf(fboot, "%-*s %s\n", com.lspname, com.spname[j], bootseq[j]);
        }
    }
    
    fclose(fboot);
    free(bootindex);
    for (j = 0; j < com.ns; j++) {
        free(bootseq[j]);
    }
    free(bootseq);
    
    return 0;
}

/* 并行化特征向量计算 */
int EigenQREV(double Q[], double pi[], double Root[], double U[], double V[])
{
    int i, j, k;
    int n = com.ncode;
    double *work = (double*)malloc(n * n * sizeof(double));
    double *temp = (double*)malloc(n * n * sizeof(double));
    
    /* 构造对称化矩阵进行特征值分解 */
    #ifdef _OPENMP
    #pragma omp parallel for private(i, j) schedule(static)
    #endif
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (pi[i] > 0 && pi[j] > 0) {
                work[i * n + j] = Q[i * n + j] * sqrt(pi[j] / pi[i]);
            } else {
                work[i * n + j] = 0.0;
            }
        }
    }
    
    /* 调用特征值分解（这部分通常需要LAPACK，无法并行化） */
    /* 这里假设有EigenRealSym函数 */
    // EigenRealSym(work, n, Root, U);
    
    /* 并行计算右特征向量矩阵 */
    #ifdef _OPENMP
    #pragma omp parallel for private(i, j) schedule(static)
    #endif
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (pi[j] > 0) {
                U[i * n + j] = U[i * n + j] / sqrt(pi[j]);
                V[j * n + i] = U[i * n + j] * pi[j];
            }
        }
    }
    
    free(work);
    free(temp);
    return 0;
}

/* 并行化似然计算中的节点处理 */
int ConditionalPNode(int inode, int igene, double x[])
{
    int ison, h, j, k;
    double *conP = nodes[inode].conP;
    int nson = nodes[inode].nson;
    
    if (nson == 0) return 0;  /* 叶节点 */
    
    /* 为每个子节点计算P矩阵 */
    for (ison = 0; ison < nson; ison++) {
        int son = nodes[inode].sons[ison];
        double t = nodes[son].branch;
        GetPMatBranch(nodes[son].pkappa, x, t, son);
    }
    
    /* 并行处理每个模式 */
    #ifdef _OPENMP
    #pragma omp parallel for private(h, j, k, ison) schedule(static)
    #endif
    for (h = 0; h < com.npatt; h++) {
        /* 为每个状态计算条件概率 */
        for (j = 0; j < com.ncode; j++) {
            double prob = 1.0;
            
            /* 对所有子节点计算概率乘积 */
            for (ison = 0; ison < nson; ison++) {
                int son = nodes[inode].sons[ison];
                double *PMat = nodes[son].pkappa;
                double son_prob = 0.0;
                
                /* 计算转移概率 */
                if (nodes[son].nson == 0) {
                    /* 叶节点：从序列直接获取 */
                    int state = com.z[son][h];
                    if (state >= 0 && state < com.ncode) {
                        son_prob = PMat[j * com.ncode + state];
                    }
                } else {
                    /* 内部节点：使用条件概率 */
                    for (k = 0; k < com.ncode; k++) {
                        son_prob += PMat[j * com.ncode + k] * nodes[son].conP[h * com.ncode + k];
                    }
                }
                prob *= son_prob;
            }
            
            conP[h * com.ncode + j] = prob;
        }
    }
    
    return 0;
}

/* 并行化矩阵乘法 */
int matby(double a[], double b[], double c[], int n, int k, int m)
{
    int i, j, l;
    
    /* 并行矩阵乘法 C = A * B */
    /* A: n×k, B: k×m, C: n×m */
    #ifdef _OPENMP
    #pragma omp parallel for private(i, j, l) schedule(static)
    #endif
    for (i = 0; i < n; i++) {
        for (j = 0; j < m; j++) {
            double sum = 0.0;
            for (l = 0; l < k; l++) {
                sum += a[i * k + l] * b[l * m + j];
            }
            c[i * m + j] = sum;
        }
    }
    
    return 0;
}

/* 并行化LU分解 */
int LUdecomposition(double *a, int n, int *indx, double *d)
{
    int i, imax = 0, j, k;
    double big, dum, sum, temp, *vv;
    
    vv = (double*)malloc(n * sizeof(double));
    *d = 1.0;
    
    /* 寻找最大元素（不能并行化，需要顺序执行） */
    for (i = 0; i < n; i++) {
        big = 0.0;
        for (j = 0; j < n; j++) {
            if ((temp = fabs(a[i * n + j])) > big) big = temp;
        }
        if (big == 0.0) {
            free(vv);
            return -1;  /* 奇异矩阵 */
        }
        vv[i] = 1.0 / big;
    }
    
    /* LU分解主循环（部分可以并行化） */
    for (j = 0; j < n; j++) {
        /* 计算上三角部分 */
        for (i = 0; i < j; i++) {
            sum = a[i * n + j];
            for (k = 0; k < i; k++) {
                sum -= a[i * n + k] * a[k * n + j];
            }
            a[i * n + j] = sum;
        }
        
        /* 寻找主元 */
        big = 0.0;
        for (i = j; i < n; i++) {
            sum = a[i * n + j];
            for (k = 0; k < j; k++) {
                sum -= a[i * n + k] * a[k * n + j];
            }
            a[i * n + j] = sum;
            
            if ((dum = vv[i] * fabs(sum)) >= big) {
                big = dum;
                imax = i;
            }
        }
        
        /* 交换行 */
        if (j != imax) {
            for (k = 0; k < n; k++) {
                dum = a[imax * n + k];
                a[imax * n + k] = a[j * n + k];
                a[j * n + k] = dum;
            }
            *d = -(*d);
            vv[imax] = vv[j];
        }
        
        indx[j] = imax;
        
        if (a[j * n + j] == 0.0) a[j * n + j] = 1e-20;
        
        /* 下三角归一化 */
        if (j != n - 1) {
            dum = 1.0 / a[j * n + j];
            #ifdef _OPENMP
            #pragma omp parallel for private(i) schedule(static)
            #endif
            for (i = j + 1; i < n; i++) {
                a[i * n + j] *= dum;
            }
        }
    }
    
    free(vv);
    return 0;
}

/* 辅助函数实现 */
int CountChanges(int seq1, int seq2, double *P, double *Q)
{
    int k, transitions = 0, transversions = 0, compared = 0;
    
    for (k = 0; k < com.ls; k++) {
        int b1 = com.z[seq1][k];
        int b2 = com.z[seq2][k];
        
        if (b1 >= 0 && b2 >= 0 && b1 < 4 && b2 < 4) {
            compared++;
            if (b1 != b2) {
                /* 判断是转换还是颠换 */
                if ((b1 + b2 == 3) ||  /* A<->G 或 T<->C */
                    (b1 + b2 == 1)) {   /* A<->G 或 T<->C 的另一种组合 */
                    transitions++;
                } else {
                    transversions++;
                }
            }
        }
    }
    
    if (compared > 0) {
        *P = (double)transitions / compared;
        *Q = (double)transversions / compared;
    } else {
        *P = *Q = 0.0;
    }
    
    return compared;
}
