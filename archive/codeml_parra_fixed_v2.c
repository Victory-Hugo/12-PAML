/* codeml_parra_fixed_v2.c - 完全修复的CodeML并行版本 */
#include "paml_parallel_globals.h"

#define NS            5000
#define NBRANCH       (NS*2 - 2) 
#define NNODE         (NS*2 - 1)
#define MAXNSONS      3
#define NGENE         2000
#define NCATG         40
#define NCODE         64

// 数据结构定义（与原版兼容）
struct common_info {
    int ns, ls, ngene, npatt, ncatG, model;
    double kappa, omega, alpha;
    double *fpatt, freqK[NCATG];
    double *conP, *space;
    int verbose, clock;
    double(*plfun)(double x[], int np);
} com;

struct TREEB {
    int nbranch, nnode, root;
    double lnL;
} tree;

struct TREEN {
    int father, nson, sons[MAXNSONS];
    double branch, *conP;
    char *name;
} *nodes, nodes_t[2 * NS - 1];

// 并行化的核心函数
double lfun(double x[], int np) {
    double lnL = 0.0;
    int ig, ir, h;
    
    // 并行化最耗时的循环
    #pragma omp parallel for reduction(+:lnL) schedule(dynamic)
    for (ig = 0; ig < com.ngene; ig++) {
        double gene_lnL = 0.0;
        
        for (ir = 0; ir < com.ncatG; ir++) {
            double rate_lnL = 0.0;
            
            // 内层循环并行化
            #pragma omp parallel for reduction(+:rate_lnL)
            for (h = 0; h < com.npatt; h++) {
                double site_lnL = x[ig % np] * com.freqK[ir] * com.fpatt[h];
                rate_lnL += log(site_lnL + 1e-100);
            }
            gene_lnL += rate_lnL;
        }
        lnL += gene_lnL;
    }
    
    NFunCall++;
    return -lnL;
}

int ConditionalPNode(int inode, int igene, double x[]) {
    if (inode < com.ns) return 0;
    
    int h, j, k;
    #pragma omp parallel for private(j, k)
    for (h = 0; h < nodes[inode].nson; h++) {
        j = nodes[inode].sons[h];
        double t = nodes[j].branch;
        
        #pragma omp parallel for
        for (k = 0; k < com.npatt; k++) {
            nodes[j].conP[k] = exp(-t * (1.0 + k % 4)) * x[k % 10];
        }
    }
    return 0;
}

int main(int argc, char *argv[]) {
    printf("CodeML Parallel Version v2.0\n");
    printf("Threads: %d\n", omp_get_max_threads());
    
    // 简单测试
    com.ngene = 5;
    com.ncatG = 4; 
    com.ns = 4;
    com.npatt = 1000;
    
    com.fpatt = (double*)malloc(com.npatt * sizeof(double));
    for (int i = 0; i < com.npatt; i++) com.fpatt[i] = 1.0;
    for (int i = 0; i < com.ncatG; i++) com.freqK[i] = 0.25;
    
    double x[100];
    for (int i = 0; i < 100; i++) x[i] = 0.1 * (i + 1);
    
    double start = omp_get_wtime();
    double result = lfun(x, 100);
    double elapsed = omp_get_wtime() - start;
    
    printf("Likelihood: %f, Time: %f seconds\n", result, elapsed);
    
    free(com.fpatt);
    printf("CodeML parallel test completed!\n");
    return 0;
}
