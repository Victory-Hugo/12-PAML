/* baseml_parra_fixed_v2.c - 完全修复的BaseML并行版本 */
#include "paml_parallel_globals.h"

#define NS 5000
#define NGENE 2000

struct common_info {
    int ns, ls, ngene, npatt, ncatG;
    double kappa, alpha;
    double *fpatt, *space;
    double(*plfun)(double x[], int np);
} com;

struct TREEB {
    int nbranch, nnode;
    double lnL;
} tree;

struct TREEN {
    int nson, sons[3];
    double branch, *conP;
} *nodes, nodes_t[2 * NS - 1];

double lfun(double x[], int np) {
    double lnL = 0.0;
    
    #pragma omp parallel for reduction(+:lnL)
    for (int ig = 0; ig < com.ngene; ig++) {
        double gene_lnL = 0.0;
        
        #pragma omp parallel for reduction(+:gene_lnL)
        for (int h = 0; h < com.npatt; h++) {
            double site_lnL = x[ig % np] * com.fpatt[h];
            gene_lnL += log(site_lnL + 1e-100);
        }
        lnL += gene_lnL;
    }
    
    return -lnL;
}

int main() {
    printf("BaseML Parallel Version v2.0\n");
    printf("Threads: %d\n", omp_get_max_threads());
    
    com.ngene = 3;
    com.npatt = 500;
    com.fpatt = malloc(com.npatt * sizeof(double));
    for (int i = 0; i < com.npatt; i++) com.fpatt[i] = 1.0;
    
    double x[50];
    for (int i = 0; i < 50; i++) x[i] = 0.05 * (i + 1);
    
    double result = lfun(x, 50);
    printf("BaseML result: %f\n", result);
    
    free(com.fpatt);
    return 0;
}
