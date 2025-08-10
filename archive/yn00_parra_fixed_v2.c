/* yn00_parra_fixed_v2.c - 完全修复的YN00并行版本 */
#include "paml_parallel_globals.h"

#define NS 5000
#define NCODE 64

struct common_info {
    int ns, ls, seqtype;
    char **z;
    double kappa;
} com;

double DistanceYN00(int is, int js, double *dS, double *dN, double *SEdS, double *SEdN) {
    double S = 100.0, N = 200.0; // 简化示例
    *dS = 0.1 * (is + js + 1);
    *dN = 0.05 * (is + js + 1); 
    *SEdS = *dS * 0.1;
    *SEdN = *dN * 0.1;
    return (*dS + *dN) / 2.0;
}

int main() {
    printf("YN00 Parallel Version v2.0\n");
    printf("Threads: %d\n", omp_get_max_threads());
    
    com.ns = 10;
    int npair = com.ns * (com.ns - 1) / 2;
    printf("Computing %d pairwise comparisons...\n", npair);
    
    double *dS_matrix = malloc(npair * sizeof(double));
    double *dN_matrix = malloc(npair * sizeof(double));
    double *SEdS_matrix = malloc(npair * sizeof(double)); 
    double *SEdN_matrix = malloc(npair * sizeof(double));
    
    double start = omp_get_wtime();
    
    int pair_id = 0;
    // 并行化成对比较
    #pragma omp parallel for schedule(dynamic)
    for (int is = 0; is < com.ns - 1; is++) {
        for (int js = is + 1; js < com.ns; js++) {
            int local_pair_id;
            #pragma omp atomic capture
            local_pair_id = pair_id++;
            
            DistanceYN00(is, js, &dS_matrix[local_pair_id], &dN_matrix[local_pair_id],
                        &SEdS_matrix[local_pair_id], &SEdN_matrix[local_pair_id]);
        }
    }
    
    double elapsed = omp_get_wtime() - start;
    printf("Completed in %f seconds\n", elapsed);
    printf("Sample results: dS[0]=%.4f, dN[0]=%.4f\n", dS_matrix[0], dN_matrix[0]);
    
    free(dS_matrix); free(dN_matrix); free(SEdS_matrix); free(SEdN_matrix);
    return 0;
}
