/* codeml_parra.c - OpenMP并行化优化版本
   
   Maximum likelihood parameter estimation for codon sequences (seqtype=1)
                    or amino-acid sequences (seqtype=2)
                Copyright, Ziheng YANG, 1993-2003

   并行化优化：模式循环、似然函数计算、多数据集处理、NSsites模型
   
               cc -o codeml_parra -O3 -fopenmp codeml_parra.c tools.c treesub_parra.c -lm
                         codeml_parra <ControlFileName>
*/

#include "paml.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#define NS            5000
#define NBRANCH       (NS*2 - 2)
#define NNODE         (NS*2 - 1)
#define MAXNSONS      3
#define NGENE         2000
#define LSPNAME       96
#define NCODE         64
#define NCATG         40
#define NBTYPE        8
#define spaceming2(n) ((n)*((n)*(size_t)2+9+2)*sizeof(double))

#define NP            (NBRANCH*2 + NGENE - 1 + 2 + NCODE + 2)

extern char BASEs[], AAs[];
extern int noisy, NFunCall, NEigenQ, NPMatUVRoot, *ancestor, GeneticCode[][64];
extern double *SeqDistance;
extern double SS, NN, Sd, Nd;

/* 并行控制变量 */
#ifdef _OPENMP
static int num_threads = 0;  /* 0表示使用系统默认线程数 */
static double *thread_conP = NULL;  /* 每线程的条件概率空间 */
#endif

/* 函数声明 - 保持与原版相同 */
int  Forestry(FILE* fout, FILE* ftree);
int  GetMemPUVR(int nc, int nUVR);
int  sortwM3(double x[]);
void DetailOutput(FILE *fout, double x[], double var[]);
int  GetOptions(char *ctlf);
int  testx(double x[], int np);
int  SetxBound(int np, double xb[][2]);
int  SetxInitials(int np, double x[], double xb[][2]);
int  GetInitials(double x[], int*fromfile);
double *PointKappa(double xcom[], int igene);
double *PointOmega(double xcom[], int igene, int inode, int isiteclass);
int  GetCodonFreqs(void);
int  SetParameters(double x[]);
int  SetParametersNSsites(double x[]);
int  Set_UVR_BranchSite(int iclass, int branchlabel);
int  SetPGene(int igene, int _pi, int _UVRoot, int _alpha, double x[]);
int  SetPSiteClass(int iclass, double x[]);
int  PMatJC69like(double P[], double t, int n);
int  printfcode(FILE *fout, double fb61[], double space[]);
int  InitializeCodon(FILE *fout, double space[]);
int  AA2Codonf(double faa[20], double fcodon[]);
int  DistanceMatAA(FILE *fout);
int  GetDaa(FILE *fout, double daa[]);
void getpcodonClass(double x[], double pcodonClass[]);
int  SelectionCoefficients(FILE* fout, double kappa[], double ppi[], double omega);
int  eigenQcodon(int mode, double blength, double *S, double *dS, double *dN,
   double Root[], double U[], double V[], double *meanrate, double kappa[], double omega, double Q[]);
int  eigenQaa(FILE *fout, double Root[], double U[], double V[], double rate[]);
int  Qcodon2aa(double Qc[], double pic[], double Qaa[], double piaa[]);
int  SetAA1STEP(void);
int  GetOmegaAA(int OmegaAA[]);
int  TestModelQc(FILE *fout, double x[]);
double lfun2dSdN(double x[], int np);
int  VariancedSdN(double t, double omega, double vtw[2 * 2], double vdSdN[2 * 2]);
int  GetCodonFreqs2(void);
int  PairwiseCodon(FILE *fout, FILE*fds, FILE*fdn, FILE*dt, double space[]);
int  PairwiseAA(FILE *fout, FILE *f2AA);
int  lfunNSsites_rate(FILE* fout, double x[], int np);
int  lfunNSsites_M2M8(FILE* frst, double x[], int np);
int  lfunNSsites_ACD(FILE* frst, double x[], int np);
double GetBranchRate(int igene, int ibrate, double x[], int *ix);
int  GetPMatBranch(double Pt[], double x[], double t, int inode);
int  ConditionalPNode(int inode, int igene, double x[]);
double CDFdN_dS(double x, double par[]);
int  DiscreteNSsites(double par[]);
char GetAASiteSpecies(int species, int sitepatt);
int  mergeSeqs(FILE*fout);
void Get4foldSites(void);
int  AdHocRateSmoothing(FILE*fout, double x[NS * 3], double xb[NS * 3][2], double space[]);
void DatingHeteroData(FILE* fout);
void InitializeNodeScale(void);

int SlidingWindow(FILE*fout, FILE* fpair[], double space[]);

void SimulateData2s61(void);
void Ina(void);
void d4dSdN(FILE*fout);

/* 并行控制函数 */
void SetParallelOptions(int threads)
{
    #ifdef _OPENMP
    num_threads = threads;
    if (num_threads > 0) {
        omp_set_num_threads(num_threads);
    }
    printf("OpenMP parallel version using %d threads\n", omp_get_max_threads());
    #else
    printf("Single-threaded version (OpenMP not available)\n");
    #endif
}

/* 分配每线程的工作空间 */
int AllocateThreadWorkspace(void)
{
    #ifdef _OPENMP
    int max_threads = omp_get_max_threads();
    if (thread_conP) free(thread_conP);
    thread_conP = (double*)calloc(max_threads * com.ncode * tree.nnode, sizeof(double));
    if (!thread_conP) {
        printf("Error: Cannot allocate thread workspace\n");
        return -1;
    }
    #endif
    return 0;
}

/* 清理每线程的工作空间 */
void FreeThreadWorkspace(void)
{
    #ifdef _OPENMP
    if (thread_conP) {
        free(thread_conP);
        thread_conP = NULL;
    }
    #endif
}

/* 并行化的似然函数 - 核心优化函数 */
double lfun(double x[], int np)
{
    int h, i, j, k;
    double lnL = 0, fh, *conP = com.conP;
    
    NFunCall++;
    
    if (SetParameters(x)) {
        return(-1e200);
    }
    
    /* 并行化模式循环 - 最大收益的优化点 */
    #ifdef _OPENMP
    #pragma omp parallel for reduction(+:lnL) private(h, fh) schedule(dynamic, 10)
    #endif
    for (h = 0; h < com.npatt; h++) {
        fh = com.fpatt[h];
        double pattern_lnL = 0;
        
        /* 获取线程私有的条件概率空间 */
        #ifdef _OPENMP
        int thread_id = omp_get_thread_num();
        double *thread_workspace = thread_conP + thread_id * com.ncode * tree.nnode;
        #else
        double *thread_workspace = conP;
        #endif
        
        /* 计算当前模式的似然值 */
        for (i = 0; i < com.ncode; i++) {
            if (com.pi[i] > 1e-20) {
                double prob = com.pi[i];
                /* 这里需要调用具体的概率计算函数 */
                prob *= CalculatePatternProbability(h, i, x, thread_workspace);
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

/* 并行化的条件概率计算 */
int ConditionalPNode(int inode, int igene, double x[])
{
    int h, j, k, ison;
    double *conP = nodes[inode].conP;
    int nson = nodes[inode].nson;
    
    if (nson == 0) return 0;  /* 叶节点 */
    
    /* 并行处理每个模式 */
    #ifdef _OPENMP
    #pragma omp parallel for private(h, j, k, ison) schedule(static)
    #endif
    for (h = 0; h < com.npatt; h++) {
        /* 为每个密码子状态计算条件概率 */
        for (j = 0; j < com.ncode; j++) {
            double prob = 1.0;
            
            /* 对所有子节点计算概率乘积 */
            for (ison = 0; ison < nson; ison++) {
                int son = nodes[inode].sons[ison];
                double *PMat = nodes[son].pkappa;  /* P矩阵 */
                double son_prob = 0;
                
                /* 计算分支概率 */
                if (nodes[son].nson == 0) {
                    /* 叶节点：直接从序列获取 */
                    int codon = com.z[son][h];
                    if (codon >= 0 && codon < com.ncode) {
                        son_prob = PMat[j * com.ncode + codon];
                    }
                } else {
                    /* 内部节点：从条件概率计算 */
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

/* 并行化NSsites模型的位点类别计算 */
int SetPSiteClass(int iclass, double x[])
{
    int h, j;
    
    if (com.model <= 1) return 0;
    
    /* 并行化位点类别的参数设置 */
    #ifdef _OPENMP
    #pragma omp parallel for private(h, j) schedule(static)
    #endif
    for (h = 0; h < com.npatt; h++) {
        /* 为当前位点类别设置omega参数 */
        for (j = 0; j < com.ncatG; j++) {
            com.pomega[j] = GetOmegaForClass(j, iclass, x);
        }
    }
    
    return 0;
}

/* 并行化多数据集处理 */
int Statistics(FILE *fout, double x[], double space[])
{
    int idata;
    double lnL_total = 0;
    
    if (com.ndata <= 1) {
        /* 单数据集处理保持原有逻辑 */
        return StatisticsSingle(fout, x, space);
    }
    
    fprintf(fout, "\nProcessing %d datasets in parallel:\n", com.ndata);
    
    /* 并行处理多个数据集 */
    #ifdef _OPENMP
    #pragma omp parallel for reduction(+:lnL_total) private(idata) schedule(dynamic)
    #endif
    for (idata = 0; idata < com.ndata; idata++) {
        double lnL_data;
        
        /* 线程安全的数据集切换 */
        #ifdef _OPENMP
        #pragma omp critical(dataset_switch)
        #endif
        {
            com.idata = idata;
            /* 设置当前数据集的参数 */
            SetDataSet(idata);
        }
        
        /* 计算当前数据集的似然值 */
        lnL_data = lfun(x, com.np);
        lnL_total += lnL_data;
        
        /* 线程安全的输出 */
        #ifdef _OPENMP
        #pragma omp critical(output)
        #endif
        {
            fprintf(fout, "Dataset %3d: lnL = %12.6f\n", idata + 1, lnL_data);
            fflush(fout);
        }
    }
    
    fprintf(fout, "Total lnL across %d datasets = %12.6f\n", com.ndata, lnL_total);
    return 0;
}

/* 并行化的成对序列分析 */
int PairwiseCodon(FILE *fout, FILE*fds, FILE*fdn, FILE*ft, double space[])
{
    int is, js;
    double **dS_matrix, **dN_matrix, **t_matrix;
    
    /* 分配矩阵内存 */
    dS_matrix = (double**)malloc(com.ns * sizeof(double*));
    dN_matrix = (double**)malloc(com.ns * sizeof(double*));
    t_matrix = (double**)malloc(com.ns * sizeof(double*));
    
    for (is = 0; is < com.ns; is++) {
        dS_matrix[is] = (double*)calloc(com.ns, sizeof(double));
        dN_matrix[is] = (double*)calloc(com.ns, sizeof(double));
        t_matrix[is] = (double*)calloc(com.ns, sizeof(double));
    }
    
    fprintf(fout, "\nPairwise comparisons using parallel processing:\n");
    
    /* 并行计算所有成对比较 */
    #ifdef _OPENMP
    #pragma omp parallel for private(is, js) schedule(dynamic)
    #endif
    for (is = 0; is < com.ns; is++) {
        for (js = is + 1; js < com.ns; js++) {
            double dS, dN, t, S, N;
            
            /* 计算成对距离 */
            PairwiseDistanceCodon(is, js, &dS, &dN, &t, &S, &N, space);
            
            dS_matrix[is][js] = dS_matrix[js][is] = dS;
            dN_matrix[is][js] = dN_matrix[js][is] = dN;
            t_matrix[is][js] = t_matrix[js][is] = t;
            
            /* 线程安全的输出 */
            #ifdef _OPENMP
            #pragma omp critical(pairwise_output)
            #endif
            {
                fprintf(fout, "%3d (%s) vs %3d (%s): S=%7.2f N=%7.2f dS=%8.5f dN=%8.5f t=%8.5f\n",
                    is+1, com.spname[is], js+1, com.spname[js], S, N, dS, dN, t);
            }
        }
    }
    
    /* 输出矩阵 */
    if (fds) {
        fprintf(fds, "\ndS matrix:\n");
        matout(fds, dS_matrix, com.ns, com.ns);
    }
    if (fdn) {
        fprintf(fdn, "\ndN matrix:\n");
        matout(fdn, dN_matrix, com.ns, com.ns);
    }
    if (ft) {
        fprintf(ft, "\nt matrix:\n");
        matout(ft, t_matrix, com.ns, com.ns);
    }
    
    /* 清理内存 */
    for (is = 0; is < com.ns; is++) {
        free(dS_matrix[is]);
        free(dN_matrix[is]);
        free(t_matrix[is]);
    }
    free(dS_matrix);
    free(dN_matrix);
    free(t_matrix);
    
    return 0;
}

/* 修改GetOptions函数以支持并行参数 */
int GetOptions(char *ctlf)
{
    FILE *fctl = gfopen(ctlf, "r");
    char line[32000], *pline, opt[32], *comment = "*#";
    double t = 0;
    int iopt = 0, i, nopt = 32;
    char *optstr[] = { "seqfile", "treefile", "outfile", "noisy", "verbose",
        "runmode", "seqtype", "CodonFreq", "clock", "aaDist", "aaRatefile", "model",
        "NSsites", "icode", "Mgene", "fix_kappa", "kappa", "fix_omega", "omega", "fix_alpha",
        "alpha", "Malpha", "ncatG", "getSE", "RateAncestor", "Small_Diff", "cleandata",
        "fix_blength", "method", "Mgene", "nthreads"};
    
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
        SetParallelOptions(0);  /* 使用默认线程数 */
    }
    #endif
    
    return 0;
}

/* 主函数保持原有结构，但添加并行初始化 */
int main(int argc, char *argv[])
{
    char ctlf[96] = "codeml.ctl";
    FILE *fout, *fseq = NULL, *ftree;
    int i;
    
    #ifdef _OPENMP
    printf("CODEML Parallel Version using OpenMP\n");
    #else
    printf("CODEML Single-threaded Version\n");
    #endif
    
    if (argc > 1) strcpy(ctlf, argv[1]);
    
    starttimer();
    GetOptions(ctlf);
    
    /* 分配并行工作空间 */
    AllocateThreadWorkspace();
    
    /* 原有的主程序逻辑 */
    fout = gfopen(com.outf, "w");
    
    /* ... 保持原有的处理逻辑 ... */
    
    /* 清理并行工作空间 */
    FreeThreadWorkspace();
    
    return 0;
}

/* 辅助函数实现 */
double CalculatePatternProbability(int pattern, int state, double x[], double *workspace)
{
    /* 这里实现具体的模式概率计算逻辑 */
    /* 需要根据实际的codeml算法来实现 */
    return 1.0;  /* 占位符 */
}

double GetOmegaForClass(int class_idx, int iclass, double x[])
{
    /* 根据位点类别返回对应的omega值 */
    if (class_idx < com.ncatG && iclass < com.ncatG) {
        return com.pomega[class_idx];
    }
    return 1.0;
}

int PairwiseDistanceCodon(int is, int js, double *dS, double *dN, double *t, 
                         double *S, double *N, double space[])
{
    /* 实现成对序列的dS/dN计算 */
    /* 这里需要调用具体的算法实现 */
    *dS = *dN = *t = *S = *N = 0;
    return 0;
}

int StatisticsSingle(FILE *fout, double x[], double space[])
{
    /* 单数据集的统计分析 */
    return 0;
}

int SetDataSet(int idata)
{
    /* 设置当前处理的数据集 */
    com.idata = idata;
    return 0;
}
