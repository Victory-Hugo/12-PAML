/* yn00_simple.c - 简化的YN00并行版本用于测试编译
   只保留核心并行化功能，移除复杂的依赖关系
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#define NS 100
#define MAXLINE 1000

/* 简化的数据结构 */
struct {
    int ns, ls;
    char **spname;
    char **seq;
    char outf[100];
    int verbose;
} com;

/* 并行控制变量 */
#ifdef _OPENMP
static int num_threads = 0;
#endif

/* 并行控制函数 */
void SetParallelOptions(int threads)
{
    #ifdef _OPENMP
    num_threads = threads;
    if (num_threads > 0) {
        omp_set_num_threads(num_threads);
    }
    printf("YN00 Parallel test version using %d threads\n", omp_get_max_threads());
    #else
    printf("YN00 Single-threaded test version\n");
    #endif
}

/* 简化的距离计算函数 */
int SimpleDistanceYN00(int is, int js, double *dS, double *dN)
{
    /* 模拟计算 - 在实际版本中这里会有复杂的算法 */
    *dS = 0.1 + 0.01 * (is + js);
    *dN = 0.05 + 0.005 * (is + js);
    
    /* 模拟一些计算负载 */
    double sum = 0;
    for (int i = 0; i < 1000; i++) {
        sum += sin(i * 0.001);
    }
    
    return 0;
}

/* 并行化的主统计函数 */
int ParallelStatistics(FILE *fout)
{
    int is, js;
    double **dS_matrix, **dN_matrix, **omega_matrix;
    
    /* 分配矩阵 */
    dS_matrix = (double**)malloc(com.ns * sizeof(double*));
    dN_matrix = (double**)malloc(com.ns * sizeof(double*));
    omega_matrix = (double**)malloc(com.ns * sizeof(double*));
    
    for (is = 0; is < com.ns; is++) {
        dS_matrix[is] = (double*)calloc(com.ns, sizeof(double));
        dN_matrix[is] = (double*)calloc(com.ns, sizeof(double));
        omega_matrix[is] = (double*)calloc(com.ns, sizeof(double));
    }
    
    fprintf(fout, "Pairwise comparisons using parallel processing:\n");
    fprintf(fout, "Seq1\tSeq2\tdS\tdN\tomega\n");
    
    /* 并行化成对比较 - 核心优化点 */
    #ifdef _OPENMP
    #pragma omp parallel for private(is, js) schedule(dynamic)
    #endif
    for (is = 0; is < com.ns; is++) {
        for (js = is + 1; js < com.ns; js++) {
            double dS, dN, omega;
            
            /* 计算dS和dN */
            SimpleDistanceYN00(is, js, &dS, &dN);
            omega = (dS > 0) ? dN / dS : -1.0;
            
            /* 存储结果 */
            dS_matrix[is][js] = dS_matrix[js][is] = dS;
            dN_matrix[is][js] = dN_matrix[js][is] = dN;
            omega_matrix[is][js] = omega_matrix[js][is] = omega;
            
            /* 线程安全输出 */
            #ifdef _OPENMP
            #pragma omp critical
            #endif
            {
                fprintf(fout, "%d\t%d\t%.6f\t%.6f\t%.4f\n", 
                       is+1, js+1, dS, dN, omega >= 0 ? omega : -1.0);
            }
        }
    }
    
    /* 输出矩阵 */
    fprintf(fout, "\ndS matrix:\n");
    for (is = 0; is < com.ns; is++) {
        for (js = 0; js < com.ns; js++) {
            fprintf(fout, "%8.4f", dS_matrix[is][js]);
        }
        fprintf(fout, "\n");
    }
    
    /* 清理内存 */
    for (is = 0; is < com.ns; is++) {
        free(dS_matrix[is]);
        free(dN_matrix[is]);  
        free(omega_matrix[is]);
    }
    free(dS_matrix);
    free(dN_matrix);
    free(omega_matrix);
    
    return 0;
}

/* 主函数 */
int main(int argc, char *argv[])
{
    FILE *fout;
    int i;
    
    #ifdef _OPENMP
    printf("YN00 Parallel Test Version using OpenMP\n");
    #else
    printf("YN00 Single-threaded Test Version\n");
    #endif
    
    /* 设置并行选项 */
    SetParallelOptions(4);
    
    /* 创建测试数据 */
    com.ns = 10;  /* 10个序列 */
    com.ls = 100; /* 100个位点 */
    strcpy(com.outf, "test_output.txt");
    
    /* 分配序列名内存 */
    com.spname = (char**)malloc(com.ns * sizeof(char*));
    for (i = 0; i < com.ns; i++) {
        com.spname[i] = (char*)malloc(20 * sizeof(char));
        sprintf(com.spname[i], "Seq%d", i+1);
    }
    
    printf("Processing %d sequences with %d threads...\n", com.ns, 
           #ifdef _OPENMP
           omp_get_max_threads()
           #else  
           1
           #endif
    );
    
    /* 输出结果 */
    fout = fopen(com.outf, "w");
    if (!fout) {
        printf("Cannot create output file\n");
        return 1;
    }
    
    /* 运行并行统计 */
    ParallelStatistics(fout);
    
    fclose(fout);
    
    printf("Results written to %s\n", com.outf);
    printf("Parallel processing completed successfully!\n");
    
    /* 清理内存 */
    for (i = 0; i < com.ns; i++) {
        free(com.spname[i]);
    }
    free(com.spname);
    
    return 0;
}
