/* treesub_parra.c - TreeSub并行版本 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

int PMatUVRoot(double P[], double t, int n) {
    #pragma omp parallel for
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            int idx = i * n + j;
            P[idx] = (i == j) ? exp(-t) : (1.0 - exp(-t)) / (n - 1);
        }
    }
    return 0;
}

int main() {
    printf("TreeSub Parallel Version\n");
    printf("Threads: %d\n", omp_get_max_threads());
    
    int n = 4;
    double *P = malloc(n * n * sizeof(double));
    double t = 0.1;
    
    PMatUVRoot(P, t, n);
    printf("P[0,0] = %f\n", P[0]);
    
    free(P);
    return 0;
}
