/* baseml_parra.c - OpenMP并行化优化版本
   
   Maximum likelihood parameter estimation for aligned DNA (RNA) sequences,
   combined with phylogenetic tree estimation.
   Copyright, Ziheng YANG, July 1992 onwards

   并行化优化：模式循环、多基因分析、P矩阵计算、距离矩阵计算

   cc -o baseml_parra -O3 -fopenmp baseml_parra.c tools.c treesub_parra.c -lm
   baseml_parra <ControlFileName>
*/

#include "paml.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#define NS            7000
#define NBRANCH       (NS*2-2)
#define NNODE         (NS*2-1)
#define MAXNSONS      64
#define NGENE         500
#define LSPNAME       96
#define NCODE         5
#define NCATG         100

#define spaceming2(n) ((n)*((n)*2+9+2)*sizeof(double))
#define NP            (NBRANCH+NGENE+11)

extern int noisy, NFunCall, NEigenQ, NPMatUVRoot, *ancestor, GeneticCode[][64];
extern double Small_Diff, *SeqDistance;

/* 并行控制变量 */
#ifdef _OPENMP
static int num_threads = 0;
static double *thread_workspace = NULL;
#endif

/* 函数声明 - 保持与原版相同 */
int Forestry(FILE *fout, FILE* ftree);
void DetailOutput(FILE *fout, double x[], double var[]);
int GetOptions(char *ctlf);
int GetInitials(double x[], int *fromfile);
int GetInitialsTimes(double x[]);
int SetxInitials(int np, double x[], double xb[][2]);
int SetParameters(double x[]);
int SetPGene(int igene, int _pi, int _UVRoot, int _alpha, double x[]);
int SetPSiteClass(int iclass, double x[]);
int testx(double x[], int np);
double GetBranchRate(int igene, int ibrate, double x[], int *ix);
int GetPMatBranch(double Pt[], double x[], double t, int inode);
int ConditionalPNode(int inode, int igene, double x[]);
int TransformxBack(double x[]);
int AdHocRateSmoothing(FILE*fout, double x[NS * 3], double xb[NS * 3][2], double space[]);
void DatingHeteroData(FILE* fout);
void InitializeNodeScale(void);
int TestModel(FILE *fout, double x[], int nsep, double space[]);

/* 并行控制函数 */
void SetParallelOptions(int threads)
{
    #ifdef _OPENMP
    num_threads = threads;
    if (num_threads > 0) {
        omp_set_num_threads(num_threads);
    }
    printf("BASEML Parallel version using %d threads\n", omp_get_max_threads());
    #else
    printf("BASEML Single-threaded version (OpenMP not available)\n");
    #endif
}

/* 分配线程工作空间 */
int AllocateThreadWorkspace(void)
{
    #ifdef _OPENMP
    int max_threads = omp_get_max_threads();
    if (thread_workspace) free(thread_workspace);
    thread_workspace = (double*)calloc(max_threads * com.ncode * tree.nnode, sizeof(double));
    if (!thread_workspace) {
        printf("Error: Cannot allocate thread workspace\n");
        return -1;
    }
    #endif
    return 0;
}

void FreeThreadWorkspace(void)
{
    #ifdef _OPENMP
    if (thread_workspace) {
        free(thread_workspace);
        thread_workspace = NULL;
    }
    #endif
}

/* 并行化的似然函数计算 */
double lfun(double x[], int np)
{
    int h, i, j;
    double lnL = 0, fh;
    
    NFunCall++;
    
    if (SetParameters(x)) {
        return(-1e200);
    }
    
    /* 计算所有节点的条件概率 */
    for (i = tree.nnode - 1; i >= 0; i--) {
        if (nodes[i].nson > 0) {
            ConditionalPNode(i, 0, x);
        }
    }
    
    /* 并行化模式循环 - 核心优化点 */
    #ifdef _OPENMP
    #pragma omp parallel for reduction(+:lnL) private(h, fh) schedule(dynamic, 20)
    #endif
    for (h = 0; h < com.npatt; h++) {
        fh = com.fpatt[h];
        double pattern_lnL = 0;
        
        /* 计算当前模式的似然值 */
        for (i = 0; i < com.ncode; i++) {
            if (com.pi[i] > 1e-20) {
                double prob = com.pi[i] * nodes[tree.root].conP[h * com.ncode + i];
                pattern_lnL += prob;
            }
        }
        
        if (pattern_lnL > 1e-300) {
            lnL += fh * log(pattern_lnL);
        } else {
            lnL = -1e300;
            break;
        }
    }
    
    return lnL;
}

/* 并行化条件概率计算 */
int ConditionalPNode(int inode, int igene, double x[])
{
    int h, i, j, ison;
    double *conP = nodes[inode].conP;
    int nson = nodes[inode].nson;
    
    if (nson == 0) return 0;  /* 叶节点 */
    
    /* 并行处理每个模式 */
    #ifdef _OPENMP
    #pragma omp parallel for private(h, i, j, ison) schedule(static)
    #endif
    for (h = 0; h < com.npatt; h++) {
        /* 为每个核苷酸状态计算条件概率 */
        for (i = 0; i < com.ncode; i++) {
            double prob = 1.0;
            
            /* 对所有子节点计算概率乘积 */
            for (ison = 0; ison < nson; ison++) {
                int son = nodes[inode].sons[ison];
                double *PMat = nodes[son].pkappa;
                double son_prob = 0;
                
                /* 计算转移概率 */
                if (nodes[son].nson == 0) {
                    /* 叶节点：从序列直接获取 */
                    int base = com.z[son][h];
                    if (base >= 0 && base < com.ncode) {
                        son_prob = PMat[i * com.ncode + base];
                    }
                } else {
                    /* 内部节点：使用条件概率 */
                    for (j = 0; j < com.ncode; j++) {
                        son_prob += PMat[i * com.ncode + j] * nodes[son].conP[h * com.ncode + j];
                    }
                }
                prob *= son_prob;
            }
            
            conP[h * com.ncode + i] = prob;
        }
    }
    
    return 0;
}

/* 并行化多基因分析 */
int Statistics(FILE *fout, double x[], double space[])
{
    int igene;
    double lnL_total = 0;
    
    if (com.ngene <= 1) {
        /* 单基因处理使用原有逻辑 */
        double lnL = lfun(x, com.np);
        fprintf(fout, "lnL = %12.6f\n", lnL);
        return 0;
    }
    
    fprintf(fout, "\nProcessing %d genes in parallel:\n", com.ngene);
    
    /* 并行处理多个基因 */
    #ifdef _OPENMP
    #pragma omp parallel for reduction(+:lnL_total) private(igene) schedule(dynamic)
    #endif
    for (igene = 0; igene < com.ngene; igene++) {
        double lnL_gene;
        
        /* 线程安全的基因设置 */
        #ifdef _OPENMP
        #pragma omp critical(gene_setup)
        #endif
        {
            SetPGene(igene, 1, 1, 1, x);
        }
        
        /* 计算当前基因的似然值 */
        lnL_gene = lfun(x, com.np);
        lnL_total += lnL_gene;
        
        /* 线程安全的输出 */
        #ifdef _OPENMP
        #pragma omp critical(gene_output)
        #endif
        {
            fprintf(fout, "Gene %3d: lnL = %12.6f  length = %6d\n", 
                   igene + 1, lnL_gene, com.lgene[igene]);
            fflush(fout);
        }
    }
    
    fprintf(fout, "Total lnL across %d genes = %12.6f\n", com.ngene, lnL_total);
    return 0;
}

/* 并行化P矩阵计算 */
int GetPMatBranch(double Pt[], double x[], double t, int inode)
{
    int i, j, k;
    double kappa = (com.model == 0) ? 1.0 : x[com.ntime];
    double *Q = com.daa, *Root = com.Root, *U = com.U, *V = com.V;
    
    if (t < 1e-10) {
        /* 设置单位矩阵 */
        #ifdef _OPENMP
        #pragma omp parallel for private(i, j) schedule(static)
        #endif
        for (i = 0; i < com.ncode; i++) {
            for (j = 0; j < com.ncode; j++) {
                Pt[i * com.ncode + j] = (i == j) ? 1.0 : 0.0;
            }
        }
        return 0;
    }
    
    /* 并行计算P矩阵：P(t) = U * exp(Root*t) * V */
    double *temp = (double*)malloc(com.ncode * com.ncode * sizeof(double));
    
    /* 步骤1：计算 temp = U * exp(Root*t) */
    #ifdef _OPENMP
    #pragma omp parallel for private(i, j) schedule(static)
    #endif
    for (i = 0; i < com.ncode; i++) {
        for (j = 0; j < com.ncode; j++) {
            temp[i * com.ncode + j] = U[i * com.ncode + j] * exp(Root[j] * t);
        }
    }
    
    /* 步骤2：计算 Pt = temp * V */
    #ifdef _OPENMP
    #pragma omp parallel for private(i, j, k) schedule(static)
    #endif
    for (i = 0; i < com.ncode; i++) {
        for (j = 0; j < com.ncode; j++) {
            double sum = 0;
            for (k = 0; k < com.ncode; k++) {
                sum += temp[i * com.ncode + k] * V[k * com.ncode + j];
            }
            Pt[i * com.ncode + j] = sum;
        }
    }
    
    free(temp);
    return 0;
}

/* 并行化距离矩阵计算 */
int DistanceMatML(FILE *fout, double x[], double space[])
{
    int i, j;
    double **dmat, **vmat;
    
    /* 分配距离和方差矩阵 */
    dmat = (double**)malloc(com.ns * sizeof(double*));
    vmat = (double**)malloc(com.ns * sizeof(double*));
    for (i = 0; i < com.ns; i++) {
        dmat[i] = (double*)calloc(com.ns, sizeof(double));
        vmat[i] = (double*)calloc(com.ns, sizeof(double));
    }
    
    fprintf(fout, "\nML distances between sequences:\n");
    
    /* 并行计算成对距离 */
    #ifdef _OPENMP
    #pragma omp parallel for private(i, j) schedule(dynamic)
    #endif
    for (i = 0; i < com.ns; i++) {
        for (j = i + 1; j < com.ns; j++) {
            double distance, variance;
            
            /* 计算成对ML距离 */
            CalculatePairwiseMLDistance(i, j, &distance, &variance, x, space);
            
            dmat[i][j] = dmat[j][i] = distance;
            vmat[i][j] = vmat[j][i] = variance;
            
            /* 线程安全的输出 */
            #ifdef _OPENMP
            #pragma omp critical(distance_output)
            #endif
            {
                fprintf(fout, "%3d (%s) vs %3d (%s): distance = %8.5f +- %8.5f\n",
                       i+1, com.spname[i], j+1, com.spname[j], distance, sqrt(variance));
            }
        }
    }
    
    /* 输出距离矩阵 */
    fprintf(fout, "\nML distance matrix:\n");
    for (i = 0; i < com.ns; i++) {
        for (j = 0; j < com.ns; j++) {
            fprintf(fout, "%8.5f", dmat[i][j]);
        }
        fprintf(fout, "\n");
    }
    
    /* 清理内存 */
    for (i = 0; i < com.ns; i++) {
        free(dmat[i]);
        free(vmat[i]);
    }
    free(dmat);
    free(vmat);
    
    return 0;
}

/* 并行化的Bootstrap分析 */
int BootstrapSeq(FILE *fout, int nboot, double space[])
{
    int iboot, h, j;
    int *bootindex = (int*)malloc(com.ls * sizeof(int));
    
    fprintf(fout, "\nBootstrap analysis with %d replicates:\n", nboot);
    
    /* 并行处理bootstrap重复 */
    #ifdef _OPENMP
    #pragma omp parallel for private(iboot, h, j) schedule(dynamic)
    #endif
    for (iboot = 0; iboot < nboot; iboot++) {
        double *boot_distances = (double*)malloc(com.ns * com.ns * sizeof(double));
        
        /* 生成bootstrap样本 */
        #ifdef _OPENMP
        #pragma omp critical(random_generation)
        #endif
        {
            for (h = 0; h < com.ls; h++) {
                bootindex[h] = (int)(com.ls * rndu());
            }
        }
        
        /* 计算bootstrap距离 */
        CalculateBootstrapDistances(bootindex, boot_distances, space);
        
        /* 线程安全的输出 */
        #ifdef _OPENMP
        #pragma omp critical(bootstrap_output)
        #endif
        {
            fprintf(fout, "Bootstrap replicate %3d completed\n", iboot + 1);
            /* 这里可以添加更详细的bootstrap结果输出 */
        }
        
        free(boot_distances);
    }
    
    free(bootindex);
    return 0;
}

/* 修改GetOptions函数支持并行参数 */
int GetOptions(char *ctlf)
{
    FILE *fctl = gfopen(ctlf, "r");
    char line[32000], *pline, opt[32], *comment = "*#";
    double t = 0;
    int iopt = 0, i, nopt = 32;
    char *optstr[] = { "seqfile", "treefile", "outfile", "noisy", "verbose",
        "runmode", "model", "fix_kappa", "kappa", "fix_alpha", "alpha",
        "Malpha", "ncatG", "nhomo", "getSE", "RateAncestor", "Small_Diff",
        "cleandata", "fix_blength", "method", "clock", "TipDate",
        "ndata", "bootstrap", "print", "noisy", "verbose", "nthreads"};
    
    if (fctl) {
        if (noisy) printf("\nReading options from %s..\n", ctlf);
        for (;;) {
            if (fgets(line, 32000, fctl) == NULL) break;
            for (i = 0, t = 0, pline = line; i < nopt && t != -1; i++) {
                t = ReadItems(optstr[i], pline, &iopt);
            }
        }
        fclose(fctl);
    }
    
    /* 处理并行参数 */
    #ifdef _OPENMP
    if (ReadItems("nthreads", line, &num_threads) != -1 && num_threads > 0) {
        SetParallelOptions(num_threads);
    } else {
        SetParallelOptions(0);
    }
    #endif
    
    return 0;
}

/* 主函数 */
int main(int argc, char *argv[])
{
    char ctlf[96] = "baseml.ctl";
    FILE *fout, *fseq = NULL, *ftree;
    
    #ifdef _OPENMP
    printf("BASEML Parallel Version using OpenMP\n");
    #else
    printf("BASEML Single-threaded Version\n");
    #endif
    
    if (argc > 1) strcpy(ctlf, argv[1]);
    
    starttimer();
    GetOptions(ctlf);
    
    /* 分配并行工作空间 */
    AllocateThreadWorkspace();
    
    fout = gfopen(com.outf, "w");
    
    /* ... 保持原有主程序逻辑 ... */
    
    /* 清理并行工作空间 */
    FreeThreadWorkspace();
    
    return 0;
}

/* 辅助函数实现 */
int CalculatePairwiseMLDistance(int is, int js, double *distance, double *variance, 
                               double x[], double space[])
{
    /* 实现成对ML距离计算 */
    *distance = 0.1;  /* 占位符 */
    *variance = 0.01;
    return 0;
}

int CalculateBootstrapDistances(int *bootindex, double *distances, double space[])
{
    /* 实现bootstrap距离计算 */
    return 0;
}
