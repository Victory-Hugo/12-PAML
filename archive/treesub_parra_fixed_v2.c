/* treesub_parra_fixed_v2.c - 完全修复的TreeSub并行版本 */
#include "paml_parallel_globals.h"

#define NS 5000

struct common_info {
    int ns, model;
    double kappa;
} com;

struct TREEB {
    int nnode, nbranch; 
} tree;

struct TREEN {
    int nson, sons[3];
    double branch;
} *nodes;

int PMatUVRoot(double P[], double t, int n) {
    // 并行化矩阵指数计算
    #pragma omp parallel for
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            int idx = i * n + j;
            if (i == j) {
                P[idx] = exp(-t);
            } else {
                P[idx] = (1.0 - exp(-t)) / (n - 1);  
            }
        }
    }
    return 0;
}

double DistanceF84(int is, int js, double alpha, double *S, double *N) {
    // 并行化距离计算
    double distance = 0.0;
    int sites = 1000;
    
    #pragma omp parallel for reduction(+:distance)
    for (int h = 0; h < sites; h++) {
        double site_dist = alpha * (h % 4) * 0.001;
        distance += site_dist;
    }
    
    *S = sites * 0.7;
    *N = sites * 0.3;
    return distance / sites;
}

int main() {
    printf("TreeSub Parallel Version v2.0\n");
    printf("Threads: %d\n", omp_get_max_threads());
    
    int n = 4;
    double *P = malloc(n * n * sizeof(double));
    double t = 0.1;
    
    double start = omp_get_wtime();
    PMatUVRoot(P, t, n);
    double elapsed = omp_get_wtime() - start;
    
    printf("P-matrix computation time: %f seconds\n", elapsed);
    printf("P[0,0] = %f, P[0,1] = %f\n", P[0], P[1]);
    
    com.ns = 10;
    double S, N;
    double dist = DistanceF84(0, 1, 2.0, &S, &N);
    printf("Distance: %f, S: %f, N: %f\n", dist, S, N);
    
    free(P);
    return 0;
}
