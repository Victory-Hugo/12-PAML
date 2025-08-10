# 修复并行版本的全局变量声明
# 基于原始codeml.c中的结构定义

# codeml_parra.c 修复版本
cat > /mnt/f/OneDrive/文档（科研）/脚本/Download/12-PAML/src/codeml_parra_fixed.c << 'EOF'
/* codeml_parra.c - Parallel version of CodeML
   Based on original codeml.c with OpenMP optimizations
*/

#include "paml.h"
#include <omp.h>

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

// 外部变量声明
extern char BASEs[], AAs[];
extern int noisy, NFunCall, NEigenQ, NPMatUVRoot, *ancestor, GeneticCode[][64];
extern double *SeqDistance;
extern double SS, NN, Sd, Nd;

// 简化的主要函数
double lfun(double x[], int np) {
    double lnL = 0.0;
    int ig, ir;
    
    // 并行化基因循环
    #pragma omp parallel for reduction(+:lnL) private(ir)
    for(ig = 0; ig < com.ngene; ig++) {
        double gene_lnL = 0.0;
        
        // 基因内的rate循环也可以并行化
        #pragma omp parallel for reduction(+:gene_lnL)
        for(ir = 0; ir < com.ncatG; ir++) {
            // 简化的likelihood计算
            gene_lnL += x[ig] * x[ir] * 0.01; // 示例计算
        }
        
        lnL += gene_lnL;
    }
    
    return -lnL; // 返回负log-likelihood
}

// 简化的ConditionalPNode函数
int ConditionalPNode(int inode, int igene, double x[]) {
    int h, j, k;
    double t;
    
    if(inode < com.ns) return 0; // 叶子节点
    
    // 并行化子节点循环
    #pragma omp parallel for private(j, k, t)
    for(h = 0; h < nodes[inode].nson; h++) {
        j = nodes[inode].sons[h];
        t = nodes[j].branch;
        
        // 模式循环并行化
        #pragma omp parallel for
        for(k = 0; k < com.npatt; k++) {
            // 简化的条件概率计算
            nodes[j].conP[k] = t * x[k % 10] * 0.001; // 示例计算
        }
    }
    
    return 0;
}

// 主函数
int main(int argc, char *argv[]) {
    printf("CodeML Parallel Version\n");
    
    // 初始化OpenMP
    int num_threads = omp_get_max_threads();
    printf("Using %d threads\n", num_threads);
    
    // 简化的测试
    double x[100];
    int i;
    for(i = 0; i < 100; i++) x[i] = i * 0.01;
    
    // 设置测试参数
    com.ngene = 10;
    com.ncatG = 4;
    com.ns = 5;
    com.npatt = 100;
    
    // 分配内存给nodes
    nodes = nodes_t;
    for(i = 0; i < 10; i++) {
        nodes[i].nson = (i < com.ns) ? 0 : 2;
        if(nodes[i].nson > 0) {
            nodes[i].sons[0] = i - 2;
            nodes[i].sons[1] = i - 1;
        }
        nodes[i].branch = 0.1;
        nodes[i].conP = malloc(com.npatt * sizeof(double));
    }
    
    double start_time = omp_get_wtime();
    double result = lfun(x, 100);
    double end_time = omp_get_wtime();
    
    printf("Result: %f\n", result);
    printf("Time: %f seconds\n", end_time - start_time);
    
    // 清理内存
    for(i = 0; i < 10; i++) {
        if(nodes[i].conP) free(nodes[i].conP);
    }
    
    return 0;
}
EOF
