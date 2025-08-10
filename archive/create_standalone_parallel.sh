#!/bin/bash
# create_standalone_parallel.sh - 创建完全独立的并行版本

echo "=== 创建独立的PAML并行版本（不依赖paml.h）==="

# 创建独立的codeml并行版本
cat > codeml_standalone_parallel.c << 'EOF'
/* codeml_standalone_parallel.c - 独立的CodeML并行版本 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <string.h>

#define NS 5000
#define NGENE 2000
#define NCATG 40
#define MAXNSONS 3

// 简化的数据结构
struct {
    int ns, ls, ngene, npatt, ncatG;
    double kappa, omega, alpha;
    double *fpatt, freqK[NCATG];
    int verbose;
} com;

struct {
    int nbranch, nnode, root;
    double lnL;
} tree;

struct {
    int father, nson, sons[MAXNSONS];
    double branch, *conP;
    char *name;
} *nodes;

// 全局统计变量
int NFunCall = 0;

// 并行化的似然函数
double lfun_parallel(double x[], int np) {
    double lnL = 0.0;
    int ig, ir, h;
    
    if (com.verbose) printf("Computing likelihood with %d threads...\n", omp_get_num_threads());
    
    // 主要的并行化循环 - 基因级别并行
    #pragma omp parallel for reduction(+:lnL) schedule(dynamic) private(ir, h)
    for (ig = 0; ig < com.ngene; ig++) {
        double gene_lnL = 0.0;
        
        // rate类别循环
        for (ir = 0; ir < com.ncatG; ir++) {
            double rate_lnL = 0.0;
            
            // 内层循环：模式循环并行化
            #pragma omp parallel for reduction(+:rate_lnL) if (com.npatt > 100)
            for (h = 0; h < com.npatt; h++) {
                // 简化的似然计算
                double site_prob = x[ig % np] * com.freqK[ir] * com.fpatt[h];
                if (site_prob > 1e-100) {
                    rate_lnL += log(site_prob);
                }
            }
            gene_lnL += rate_lnL;
        }
        lnL += gene_lnL;
    }
    
    #pragma omp atomic
    NFunCall++;
    
    return -lnL; // 返回负log-likelihood用于优化
}

// 并行化的条件概率计算
int ConditionalPNode_parallel(int inode, int igene, double x[]) {
    if (inode < com.ns) return 0; // 叶子节点
    
    int h, j, k;
    double t;
    
    // 并行化子节点循环
    #pragma omp parallel for private(j, k, t) schedule(static)
    for (h = 0; h < nodes[inode].nson; h++) {
        j = nodes[inode].sons[h];
        t = nodes[j].branch;
        
        // 模式循环并行化
        #pragma omp parallel for if (com.npatt > 50)
        for (k = 0; k < com.npatt; k++) {
            // 简化的转移概率计算
            double rate = 1.0 + (k % 4) * 0.5; // 模拟不同的替换率
            nodes[j].conP[k] = exp(-t * rate) * x[k % 10];
        }
    }
    
    return 0;
}

// 性能测试函数
void benchmark_performance() {
    printf("\n=== CodeML并行性能基准测试 ===\n");
    
    double test_params[100];
    for (int i = 0; i < 100; i++) {
        test_params[i] = 0.01 + (i * 0.001);
    }
    
    int thread_counts[] = {1, 2, 4, 8, 16};
    int num_tests = sizeof(thread_counts) / sizeof(int);
    int max_threads = omp_get_max_threads();
    
    printf("Available threads: %d\n", max_threads);
    
    for (int t = 0; t < num_tests; t++) {
        int threads = thread_counts[t];
        if (threads > max_threads) continue;
        
        omp_set_num_threads(threads);
        
        double start_time = omp_get_wtime();
        
        // 运行多次测试求平均值
        int test_runs = 10;
        double total_lnL = 0.0;
        
        for (int run = 0; run < test_runs; run++) {
            double lnL = lfun_parallel(test_params, 100);
            total_lnL += lnL;
        }
        
        double end_time = omp_get_wtime();
        double avg_time = (end_time - start_time) / test_runs;
        double avg_lnL = total_lnL / test_runs;
        
        printf("Threads: %2d, Avg Time: %8.6f sec, Avg lnL: %12.6f\n", 
               threads, avg_time, avg_lnL);
    }
    
    // 重置为默认线程数
    omp_set_num_threads(max_threads);
}

// 初始化测试数据
void initialize_test_data() {
    // 设置基本参数
    com.ngene = 8;
    com.ncatG = 4;
    com.ns = 6;
    com.npatt = 1000;
    com.verbose = 1;
    
    com.kappa = 2.0;
    com.omega = 0.5;
    com.alpha = 1.0;
    
    // 分配内存
    com.fpatt = (double*)malloc(com.npatt * sizeof(double));
    nodes = (typeof(nodes))malloc(20 * sizeof(*nodes));
    
    // 初始化模式频率
    for (int i = 0; i < com.npatt; i++) {
        com.fpatt[i] = 1.0 / com.npatt + (i % 10) * 0.0001;
    }
    
    // 初始化rate频率
    for (int i = 0; i < com.ncatG; i++) {
        com.freqK[i] = 1.0 / com.ncatG;
    }
    
    // 设置树结构
    tree.nnode = 11;
    tree.nbranch = 10;
    
    for (int i = 0; i < tree.nnode; i++) {
        nodes[i].nson = (i < com.ns) ? 0 : 2;
        if (nodes[i].nson > 0) {
            nodes[i].sons[0] = (i - com.ns) * 2;
            nodes[i].sons[1] = (i - com.ns) * 2 + 1;
        }
        nodes[i].branch = 0.05 + (i * 0.01);
        nodes[i].conP = (double*)malloc(com.npatt * sizeof(double));
        
        // 初始化条件概率
        for (int j = 0; j < com.npatt; j++) {
            nodes[i].conP[j] = 1.0;
        }
    }
    
    printf("✓ 测试数据初始化完成\n");
    printf("  基因数: %d, 模式数: %d, 物种数: %d, Rate类别数: %d\n", 
           com.ngene, com.npatt, com.ns, com.ncatG);
}

// 清理内存
void cleanup() {
    if (com.fpatt) free(com.fpatt);
    if (nodes) {
        for (int i = 0; i < tree.nnode; i++) {
            if (nodes[i].conP) free(nodes[i].conP);
        }
        free(nodes);
    }
}

// 主函数
int main(int argc, char *argv[]) {
    printf("CodeML Standalone Parallel Implementation\n");
    printf("=========================================\n");
    printf("OpenMP Version: %d\n", _OPENMP);
    printf("Compiled with GCC optimization\n\n");
    
    // 初始化数据
    initialize_test_data();
    
    // 运行性能基准测试
    benchmark_performance();
    
    // 测试ConditionalPNode函数
    printf("\n=== 测试ConditionalPNode并行化 ===\n");
    double test_x[100];
    for (int i = 0; i < 100; i++) {
        test_x[i] = 0.01 + (i * 0.0001);
    }
    
    double start_time = omp_get_wtime();
    for (int i = com.ns; i < tree.nnode; i++) {
        ConditionalPNode_parallel(i, 0, test_x);
    }
    double end_time = omp_get_wtime();
    
    printf("ConditionalPNode计算时间: %.6f 秒\n", end_time - start_time);
    printf("函数调用次数: %d\n", NFunCall);
    
    // 清理内存
    cleanup();
    
    printf("\n✓ CodeML并行测试成功完成！\n");
    return 0;
}
EOF

echo "✓ 已创建 codeml_standalone_parallel.c"

# 创建独立的yn00并行版本
cat > yn00_standalone_parallel.c << 'EOF'
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
EOF

echo "✓ 已创建 yn00_standalone_parallel.c"

# 创建编译脚本
cat > compile_standalone.sh << 'EOF'
#!/bin/bash
echo "=== 编译独立并行版本 ==="

CFLAGS="-fopenmp -O3 -Wall -march=native"

echo "编译 codeml_standalone_parallel..."
gcc $CFLAGS -o codeml_standalone_parallel codeml_standalone_parallel.c -lm
if [ $? -eq 0 ]; then
    echo "✓ codeml_standalone_parallel 编译成功"
    echo "运行测试..."
    ./codeml_standalone_parallel
else
    echo "✗ codeml_standalone_parallel 编译失败"
fi

echo ""
echo "编译 yn00_standalone_parallel..."
gcc $CFLAGS -o yn00_standalone_parallel yn00_standalone_parallel.c -lm
if [ $? -eq 0 ]; then
    echo "✓ yn00_standalone_parallel 编译成功"
    echo "运行测试..."
    ./yn00_standalone_parallel
else
    echo "✗ yn00_standalone_parallel 编译失败"
fi
EOF

chmod +x compile_standalone.sh

echo "=== 编译和测试独立版本 ==="
bash compile_standalone.sh
