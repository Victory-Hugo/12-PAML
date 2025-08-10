/* yn00_parra.ctl - YN00并行版本配置文件示例
   
   适用于 yn00_parra 程序的配置参数
   包含标准参数和新增的并行控制参数
*/

      seqfile = abglobin.nuc  * sequence data filename
      outfile = yn           * main result file
      verbose = 1            * 0: concise; 1: detailed; 2: too much

        icode = 0            * 0:universal code; 1:mammalian mt; 2-10:see below
    weighting = 0            * weighting pathways between codons (0/1)?
   commonkappa = 0           * 0: different kappa for each pair; 1: one kappa
   commonf3x4  = 0           * 0: different pi for each pair; 1: one pi

* 新增的并行控制参数 *
     nthreads = 4            * OpenMP线程数: 0=自动检测, >0=指定线程数
                             * 并行化收益说明:
                             * - 成对比较数量: n*(n-1)/2 (n=序列数)
                             * - 10序列: 45个比较对
                             * - 50序列: 1225个比较对  
                             * - 100序列: 4950个比较对
                             * 序列数越多，并行化收益越显著

* Method comparison settings (when verbose >= 1)
    runmethods = 7           * Bitmask for methods to run in parallel:
                             * 1: YN00 method
                             * 2: NG86 method  
                             * 4: LWL85 method
                             * 7: All three methods (1+2+4)
                             * Methods run simultaneously using parallel sections

* Output control
   matrixout = 1             * 0: no matrices; 1: output distance matrices
                             * Parallel version outputs:
                             * - dS matrix (synonymous distances)
                             * - dN matrix (nonsynonymous distances) 
                             * - dN/dS ratio matrix
                             * - Standard error matrices

* Performance and memory settings
   chunksize = 10            * Number of sequence pairs per parallel chunk
                             * Smaller values: better load balancing
                             * Larger values: less parallel overhead
                             * Auto-adjusted based on nthreads and dataset size

* 遗传密码对照表:
* 0: 通用密码, 1: 脊椎动物线粒体, 2: 酵母线粒体, 3: 霉菌线粒体,
* 4: 无脊椎动物线粒体, 5: 纤毛虫, 6: 棘皮动物线粒体, 7: 真核线粒体,
* 8: 酵母核替代, 9: 海鞘线粒体, 10: 草履虫核

* 并行版本特性说明:
* 1. 成对序列比较并行化 - 主要性能提升点
*    - 每个序列对独立计算dS和dN
*    - 动态负载均衡，自动分配任务
*    - 支持大规模序列集(1000+序列)
*
* 2. 位点计数并行化
*    - 同义和非同义位点统计
*    - 密码子模式循环并行处理
*    - 内存访问优化，减少缓存冲突
*
* 3. 多方法并行比较
*    - YN00, NG86, LWL85方法同时运行
*    - 每个序列对的三种方法并行计算
*    - 结果自动汇总和比较分析
*
* 4. 距离矩阵并行计算
*    - 大规模矩阵运算优化
*    - 分块矩阵计算策略
*    - 内存使用效率提升
*
* 性能调优建议:
* - 小数据集(<20序列): nthreads=2-4
* - 中等数据集(20-100序列): nthreads=4-8
* - 大数据集(>100序列): nthreads=8-16
* - 内存限制: 每线程约需要 sequences²*8字节
*
* 使用示例:
* export OMP_NUM_THREADS=8
* ./yn00_parra yn00.ctl
*
* 或者直接在命令行指定:
* ./yn00_parra sequences.dat  # 直接指定序列文件
