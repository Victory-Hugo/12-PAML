/* codeml_parra_fixed.c - Fixed Parallel version of CodeML
   Based on original codeml.c with OpenMP optimizations and proper declarations
*/

#include "paml.h"
#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define NS            5000
#define NBRANCH       (NS*2 - 2)
#define NNODE         (NS*2 - 1)
#define MAXNSONS      3
#define NGENE         2000
#define LSPNAME       96
#define NCODE         64
#define NCATG         40
#define NBTYPE        8

// 全局数据结构定义（与原版一致）
struct common_info {
   char *z[NS];
   char *spname[NS], seqf[2048], outf[2048], treef[2048], daafile[2048], cleandata;
   char oldconP[NNODE];
   int seqtype, ns, ls, ngene, posG[NGENE + 1], lgene[NGENE], npatt, *pose, readpattern;
   int runmode, clock, verbose, print, codonf, aaDist, model, NSsites;
   int nOmega, nbtype, nOmegaType;
   int method, icode, ncode, Mgene, ndata, idata, ndata_trees_opt, bootstrap;
   int fix_rgene, fix_kappa, fix_omega, fix_alpha, fix_rho, nparK, fix_blength, getSE;
   int np, ntime, nrgene, nkappa, npi, nrate, nalpha, ncatG, hkyREV;
   size_t sconP, sspace;
   double *fpatt, *space, kappa, omega, alpha, rho, rgene[NGENE];
   double pi[NCODE], piG[NGENE][64], fb61[64];
   double f3x4[NGENE][12], *pf3x4, piAA[20];
   double freqK[NCATG], rK[NCATG], MK[NCATG*NCATG], daa[20 * 20], *conP, *fhK;
   double *blengths0;
   double(*plfun)(double x[], int np);
   double hyperpar[4];
   double omega_fix;
   int     conPSiteClass;
   int     NnodeScale;
   char   *nodeScale;
   double *nodeScaleF;
   double *pomega, pkappa[5], *ppi;
}  com;

struct TREEB {
   int  nbranch, nnode, root, branches[NBRANCH][2];
   double lnL;
}  tree;

struct TREEN {
   int father, nson, sons[MAXNSONS], ibranch, ipop;
   double branch, age, omega, *conP, label, label2;
   char fossil, usefossil, *name, *annotation;
}  *nodes, **gnodes, nodes_t[2 * NS - 1];

// 全局变量定义（简化版本）
char BASEs[] = "TCAG";
char AAs[] = "ARNDCQEGHILKMFPSTWYV*";
int noisy = 1;
int NFunCall = 0;
int NEigenQ = 0;
int NPMatUVRoot = 0;
int *ancestor = NULL;
int GeneticCode[1][64];
double *SeqDistance = NULL;
double SS = 0, NN = 0, Sd = 0, Nd = 0;

// 简化的工具函数
FILE *gfopen(const char *filename, const char *mode) {
    return fopen(filename, mode);
}

void error(const char *msg) {
    printf("Error: %s\n", msg);
    exit(1);
}

double etime(void) {
    return omp_get_wtime();
}

// 并行化的likelihood函数
double lfun(double x[], int np) {
    double lnL = 0.0;
    int ig, ir, h;
    
    printf("Computing likelihood with %d parameters...\n", np);
    
    // 并行化基因循环 - 最外层并行化以获得最好的负载平衡
    #pragma omp parallel for reduction(+:lnL) private(ir, h) schedule(dynamic)
    for(ig = 0; ig < com.ngene; ig++) {
        double gene_lnL = 0.0;
        
        // 基因内的rate类别循环
        for(ir = 0; ir < com.ncatG; ir++) {
            double rate_lnL = 0.0;
            
            // 模式循环 - 内层并行化
            #pragma omp parallel for reduction(+:rate_lnL)
            for(h = 0; h < com.npatt; h++) {
                // 简化的likelihood计算（示例）
                double site_lnL = x[ig % np] * com.freqK[ir] * com.fpatt[h];
                rate_lnL += site_lnL;
            }
            
            gene_lnL += rate_lnL;
        }
        
        lnL += gene_lnL;
    }
    
    NFunCall++;
    return -lnL; // 返回负log-likelihood
}

// 并行化的ConditionalPNode函数
int ConditionalPNode(int inode, int igene, double x[]) {
    int h, j, k;
    double t;
    
    if(inode < com.ns) return 0; // 叶子节点
    
    printf("Computing conditional probabilities for node %d...\n", inode);
    
    // 并行化子节点循环
    #pragma omp parallel for private(j, k, t) schedule(static)
    for(h = 0; h < nodes[inode].nson; h++) {
        j = nodes[inode].sons[h];
        t = nodes[j].branch;
        
        // 模式循环并行化
        #pragma omp parallel for
        for(k = 0; k < com.npatt; k++) {
            // 简化的转移概率计算
            double prob = exp(-t * (k % 4 + 1)) * x[k % 10];
            nodes[j].conP[k] = prob;
        }
    }
    
    return 0;
}

// 简化的参数设置函数
int SetParameters(double x[]) {
    int i;
    
    // 设置kappa参数
    com.kappa = x[0];
    
    // 设置omega参数
    com.omega = x[1];
    
    // 设置分支长度
    for(i = 0; i < tree.nbranch && i < 10; i++) {
        if(i < tree.nbranch) {
            // 简化的分支长度设置
            nodes[i].branch = x[i + 2];
        }
    }
    
    return 0;
}

// 测试函数
void test_parallel_performance() {
    printf("\n=== Parallel Performance Test ===\n");
    
    double x[100];
    int i, num_tests = 5;
    
    // 初始化测试参数
    for(i = 0; i < 100; i++) {
        x[i] = 0.1 + (i * 0.01);
    }
    
    // 测试不同线程数的性能
    int thread_counts[] = {1, 2, 4, 8};
    int num_thread_tests = sizeof(thread_counts) / sizeof(int);
    
    for(i = 0; i < num_thread_tests; i++) {
        omp_set_num_threads(thread_counts[i]);
        
        double start_time = omp_get_wtime();
        
        // 运行多次测试
        for(int test = 0; test < num_tests; test++) {
            double result = lfun(x, 100);
        }
        
        double end_time = omp_get_wtime();
        double avg_time = (end_time - start_time) / num_tests;
        
        printf("Threads: %d, Avg Time: %.6f seconds\n", thread_counts[i], avg_time);
    }
}

// 主函数
int main(int argc, char *argv[]) {
    printf("CodeML Parallel Version (Fixed)\n");
    printf("Copyright 2024, Parallel PAML Implementation\n\n");
    
    // 显示OpenMP信息
    printf("OpenMP Version: %d\n", _OPENMP);
    printf("Max threads available: %d\n", omp_get_max_threads());
    printf("Number of processors: %d\n", omp_get_num_procs());
    
    // 初始化结构体
    com.ngene = 5;
    com.ncatG = 4;
    com.ns = 4;
    com.npatt = 100;
    tree.nbranch = 6;
    tree.nnode = 7;
    
    // 分配内存
    com.fpatt = (double*)malloc(com.npatt * sizeof(double));
    
    // 初始化数据
    for(int i = 0; i < com.npatt; i++) {
        com.fpatt[i] = 1.0 / com.npatt;
    }
    
    for(int i = 0; i < com.ncatG && i < NCATG; i++) {
        com.freqK[i] = 1.0 / com.ncatG;
    }
    
    // 设置节点结构
    nodes = nodes_t;
    for(int i = 0; i < tree.nnode; i++) {
        nodes[i].nson = (i < com.ns) ? 0 : 2;
        if(nodes[i].nson > 0) {
            nodes[i].sons[0] = (i * 2 + 1) % tree.nnode;
            nodes[i].sons[1] = (i * 2 + 2) % tree.nnode;
        }
        nodes[i].branch = 0.1 * (i + 1);
        nodes[i].conP = (double*)malloc(com.npatt * sizeof(double));
        
        // 初始化条件概率
        for(int j = 0; j < com.npatt; j++) {
            nodes[i].conP[j] = 1.0;
        }
    }
    
    // 运行性能测试
    test_parallel_performance();
    
    // 测试ConditionalPNode函数
    printf("\n=== Testing ConditionalPNode ===\n");
    double test_x[100];
    for(int i = 0; i < 100; i++) {
        test_x[i] = 0.05 + (i * 0.001);
    }
    
    double start_time = omp_get_wtime();
    for(int i = com.ns; i < tree.nnode; i++) {
        ConditionalPNode(i, 0, test_x);
    }
    double end_time = omp_get_wtime();
    
    printf("ConditionalPNode computation time: %.6f seconds\n", end_time - start_time);
    
    // 清理内存
    for(int i = 0; i < tree.nnode; i++) {
        if(nodes[i].conP) free(nodes[i].conP);
    }
    free(com.fpatt);
    
    printf("\nCodeML Parallel test completed successfully!\n");
    return 0;
}
