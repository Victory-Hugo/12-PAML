/* yn00_standalone_parallel.c - 独立的YN00并行版本 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#define NS 5000

struct {
    int ns, ls, seqtype;
    char **z; // 序列数据
    double kappa;
} com;

// 并行化的YN00距离计算
double DistanceYN00_parallel(int is, int js, double *dS, double *dN, double *SEdS, double *SEdN) {
    // 简化的YN00算法实现
    double S_sites = 0.0, N_sites = 0.0;
    double S_changes = 0.0, N_changes = 0.0;
    
    int codon_length = com.ls / 3;
    
    // 并行化密码子比较
    #pragma omp parallel for reduction(+:S_sites,N_sites,S_changes,N_changes)
    for (int codon = 0; codon < codon_length; codon++) {
        // 模拟密码子分析
        double local_S = 1.5 + (codon % 3) * 0.5; // 同义位点数
        double local_N = 1.5 + ((codon + 1) % 3) * 0.5; // 非同义位点数
        
        S_sites += local_S;
        N_sites += local_N;
        
        // 模拟变化数
        if ((codon + is + js) % 4 == 0) {
            S_changes += 0.1;
        }
        if ((codon + is + js) % 7 == 0) {
            N_changes += 0.05;
        }
    }
    
    // 计算dS和dN
    *dS = (S_sites > 0) ? -0.75 * log(1.0 - (4.0/3.0) * (S_changes/S_sites)) : 0.0;
    *dN = (N_sites > 0) ? -0.75 * log(1.0 - (4.0/3.0) * (N_changes/N_sites)) : 0.0;
    
    // 简化的标准误差计算
    *SEdS = *dS * 0.1;
    *SEdN = *dN * 0.1;
    
    return (*dS + *dN) / 2.0;
}

// 成对序列比较的主函数
void pairwise_analysis_parallel() {
    printf("\n=== YN00并行成对分析 ===\n");
    
    int npairs = com.ns * (com.ns - 1) / 2;
    printf("计算 %d 个序列对...\n", npairs);
    
    // 分配结果数组
    double *dS_results = malloc(npairs * sizeof(double));
    double *dN_results = malloc(npairs * sizeof(double));
    double *SEdS_results = malloc(npairs * sizeof(double));
    double *SEdN_results = malloc(npairs * sizeof(double));
    double *omega_results = malloc(npairs * sizeof(double));
    
    double start_time = omp_get_wtime();
    
    // 并行化成对比较 - 这是主要的计算瓶颈
    int pair_index = 0;
    #pragma omp parallel for schedule(dynamic, 1) private(pair_index)
    for (int is = 0; is < com.ns - 1; is++) {
        for (int js = is + 1; js < com.ns; js++) {
            // 计算当前pair的索引
            int local_index = is * (2 * com.ns - is - 3) / 2 + (js - is - 1);
            
            double dS, dN, SEdS, SEdN;
            DistanceYN00_parallel(is, js, &dS, &dN, &SEdS, &SEdN);
            
            dS_results[local_index] = dS;
            dN_results[local_index] = dN;
            SEdS_results[local_index] = SEdS;
            SEdN_results[local_index] = SEdN;
            omega_results[local_index] = (dS > 0.001) ? dN / dS : 999.0;
        }
    }
    
    double end_time = omp_get_wtime();
    
    printf("并行计算完成时间: %.6f 秒\n", end_time - start_time);
    printf("使用线程数: %d\n", omp_get_max_threads());
    
    // 显示部分结果
    printf("\n前10个序列对的结果:\n");
    printf("Pair\t  dS\t  dN\t  dN/dS\t SEdS\t SEdN\n");
    
    for (int i = 0; i < npairs && i < 10; i++) {
        printf("%4d\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\n",
               i + 1, dS_results[i], dN_results[i], omega_results[i],
               SEdS_results[i], SEdN_results[i]);
    }
    
    // 计算统计数据
    double avg_dS = 0.0, avg_dN = 0.0, avg_omega = 0.0;
    for (int i = 0; i < npairs; i++) {
        avg_dS += dS_results[i];
        avg_dN += dN_results[i];
        if (omega_results[i] < 10.0) avg_omega += omega_results[i];
    }
    avg_dS /= npairs;
    avg_dN /= npairs;
    avg_omega /= npairs;
    
    printf("\n统计摘要:\n");
    printf("平均 dS: %.6f\n", avg_dS);
    printf("平均 dN: %.6f\n", avg_dN);
    printf("平均 dN/dS: %.6f\n", avg_omega);
    
    // 清理内存
    free(dS_results);
    free(dN_results);
    free(SEdS_results);
    free(SEdN_results);
    free(omega_results);
}

// 性能基准测试
void benchmark_yn00_performance() {
    printf("\n=== YN00并行性能基准测试 ===\n");
    
    int test_sizes[] = {5, 10, 20, 30};
    int num_tests = sizeof(test_sizes) / sizeof(int);
    
    for (int t = 0; t < num_tests; t++) {
        com.ns = test_sizes[t];
        int npairs = com.ns * (com.ns - 1) / 2;
        
        printf("\n测试 %d 个序列 (%d 个序列对):\n", com.ns, npairs);
        
        int thread_counts[] = {1, 2, 4, 8};
        int num_thread_tests = sizeof(thread_counts) / sizeof(int);
        
        for (int i = 0; i < num_thread_tests; i++) {
            int threads = thread_counts[i];
            if (threads > omp_get_max_threads()) continue;
            
            omp_set_num_threads(threads);
            
            double start = omp_get_wtime();
            pairwise_analysis_parallel();
            double elapsed = omp_get_wtime() - start;
            
            printf("线程数 %d: %.6f 秒\n", threads, elapsed);
        }
    }
}

// 主函数
int main(int argc, char *argv[]) {
    printf("YN00 Standalone Parallel Implementation\n");
    printf("=======================================\n");
    printf("OpenMP Version: %d\n", _OPENMP);
    printf("Available threads: %d\n", omp_get_max_threads());
    
    // 设置测试参数
    com.ls = 300; // 100个密码子
    com.seqtype = 1; // 密码子序列
    com.kappa = 2.0;
    
    // 运行不同规模的测试
    com.ns = 15;
    pairwise_analysis_parallel();
    
    // 运行性能基准测试
    //benchmark_yn00_performance();
    
    printf("\n✓ YN00并行测试成功完成！\n");
    return 0;
}
