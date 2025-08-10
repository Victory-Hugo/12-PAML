/* baseml_parra.ctl - BaseML并行版本配置文件示例
   
   适用于 baseml_parra 程序的配置参数
   包含标准参数和新增的并行控制参数
*/

      seqfile = /mnt/f/OneDrive/文档（科研）/脚本/Download/12-PAML/examples/Example.fasta    * sequence data file name
     treefile = /mnt/f/OneDrive/文档（科研）/脚本/Download/12-PAML/examples/Example.treefile * tree structure file name  
      outfile = /mnt/f/OneDrive/文档（科研）/脚本/Download/12-PAML/output/mlb           * main result file name

        noisy = 3    * 0,1,2,3,9: how much rubbish on the screen
      verbose = 1    * 0: concise; 1: detailed; 2: too much
      runmode = 0    * 0: user tree;  1: semi-automatic;  2: automatic
                     * 3: StepwiseAddition; (4,5):PerturbationNNI

        model = 4    * 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85
                     * 5:T92, 6:TN93, 7:REV, 8:UNREST, 9:REVu; 10:UNRESTu
                     
        Mgene = 0    * 0:rates, 1:separate; 2:pi, 3:kapa, 4:all

    fix_kappa = 0    * 1: kappa fixed, 0: kappa to be estimated
        kappa = 5    * initial or fixed kappa

    fix_alpha = 1    * 0: estimate gamma shape parameter; 1: fix it at alpha
        alpha = 0    * initial or fixed alpha, 0:infinity (constant rate)
       Malpha = 0    * 1: different alpha's for genes, 0: one alpha
        ncatG = 5    * # of categories in the dG, AdG, or nparK models of rates

        clock = 0    * 0:no clock, unrooted tree, 1:clock, rooted tree
        getSE = 0    * 0: don't want them, 1: want S.E.s of estimates
 RateAncestor = 0    * (0,1,2): rates (alpha>0) or ancestral states

   Small_Diff = 7e-6
    cleandata = 1    * remove sites with ambiguity data (1:yes, 0:no)?
*       fix_blength = -1  * 0: ignore, -1: random, 1: initial, 2: fixed, 3: proportional
        method = 0   * 0: simultaneous; 1: one branch at a time

* Molecular clock (clock = 1)
        TipDate = 0 100 * TipDate (1) & time unit
        outgroup = 1 * outgroup root for clock=1

* 新增的并行控制参数 *
     nthreads = 4    * OpenMP线程数: 0=自动检测, >0=指定线程数
                     * 建议设置说明:
                     * - 2-4核CPU: 设置为2-4
                     * - 8核CPU: 设置为4-8  
                     * - 16核CPU: 设置为8-12
                     * - 服务器: 根据负载调整

* Bootstrap and simulation settings
    bootstrap = 0    * 0: no bootstrap; n>0: number of bootstrap replicates
                     * parallel bootstrap analysis supported
        
* Additional parallel optimization notes:
* - Likelihood calculation patterns are parallelized
* - Multiple genes can be processed in parallel
* - Distance matrix calculations use parallel algorithms
* - P-matrix computations are optimized for multi-core
* - Bootstrap replicates run in parallel when bootstrap > 0
* 
* Performance tuning:
* - For large datasets (>100 sequences): increase nthreads
* - For many genes (Mgene > 0): parallel processing very effective  
* - Memory usage scales with nthreads * sequence_length
* - I/O operations remain sequential for data integrity
*
* Environment variable override:
* export OMP_NUM_THREADS=6 && ./baseml_parra baseml.ctl
