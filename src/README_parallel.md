# PAML 并行版本使用指南

## 📁 文件说明

本次创建的并行化优化文件：

- **`codeml_parra.c`** - CodeML的OpenMP并行版本
- **`baseml_parra.c`** - BaseML的OpenMP并行版本  
- **`yn00_parra.c`** - YN00的OpenMP并行版本
- **`treesub_parra.c`** - 树操作函数的并行版本
- **`Makefile_parra`** - 专用的并行版本Makefile

## 🚀 编译方法

### 1. 检查OpenMP支持
```bash
cd /mnt/f/OneDrive/文档（科研）/脚本/Download/12-PAML/src
make -f Makefile_parra check-openmp
```

### 2. 编译所有并行版本
```bash
make -f Makefile_parra all
```

### 3. 单独编译
```bash
# 编译CodeML并行版本
make -f Makefile_parra codeml_parra

# 编译BaseML并行版本  
make -f Makefile_parra baseml_parra

# 编译YN00并行版本
make -f Makefile_parra yn00_parra
```

## ⚙️ 使用方法

### 1. 设置线程数
```bash
# 方法1：环境变量设置
export OMP_NUM_THREADS=4

# 方法2：在控制文件中添加
echo "nthreads = 4" >> codeml.ctl
```

### 2. 运行程序
```bash
# 运行CodeML并行版本
./codeml_parra ../conf/codeml.ctl

# 运行BaseML并行版本
./baseml_parra ../conf/baseml.ctl

# 运行YN00并行版本
./yn00_parra sequences.dat
```

## 📊 性能基准测试

```bash
# 运行不同线程数的性能测试
make -f Makefile_parra benchmark

# 手动测试不同线程数
for threads in 1 2 4 8; do
    echo "Testing with $threads threads:"
    export OMP_NUM_THREADS=$threads
    time ./codeml_parra test.ctl
done
```

## 🔧 配置文件示例

### CodeML配置文件 (codeml_parra.ctl)
```
seqfile = sequences.dat
treefile = tree.nwk  
outfile = results.out

noisy = 3
verbose = 1
runmode = 0

seqtype = 1
CodonFreq = 2
clock = 0
aaDist = 0

model = 0
NSsites = 0

icode = 0
fix_kappa = 0
kappa = 2

fix_omega = 0
omega = 0.4

fix_alpha = 1
alpha = 0
Malpha = 0
ncatG = 8

getSE = 0
RateAncestor = 0
Small_Diff = .5e-6
cleandata = 1

* 并行参数
nthreads = 4    * OpenMP线程数，0为自动检测
```

### BaseML配置文件 (baseml_parra.ctl)
```
seqfile = dna.dat
treefile = tree.nwk
outfile = baseml.out

noisy = 3
verbose = 1
runmode = 0

model = 4    * HKY85
fix_kappa = 0
kappa = 2

fix_alpha = 1  
alpha = 0.5
Malpha = 0
ncatG = 5

getSE = 0
RateAncestor = 1
Small_Diff = .5e-6
cleandata = 1
fix_blength = 0

* 并行参数  
nthreads = 4
```

### YN00配置文件 (yn00_parra.ctl)
```
seqfile = codons.dat
outfile = yn00.out

verbose = 1
icode = 0
weighting = 0
commonkappa = 0
commonf3x4 = 0

* 并行参数
nthreads = 4
```

## 🎯 并行化优化点

### 1. CodeML并行化
- ✅ **模式循环并行化**：似然函数计算的核心循环
- ✅ **条件概率计算并行化**：节点概率计算
- ✅ **多数据集并行处理**：不同数据集同时分析
- ✅ **NSsites模型并行化**：位点类别参数计算
- ✅ **成对序列比较并行化**：dS/dN计算

### 2. BaseML并行化
- ✅ **似然函数模式循环**：DNA序列似然计算
- ✅ **多基因并行分析**：多个基因同时处理
- ✅ **P矩阵计算并行化**：转移概率矩阵
- ✅ **距离矩阵并行计算**：成对ML距离
- ✅ **Bootstrap并行化**：自举分析加速

### 3. YN00并行化
- ✅ **成对比较并行化**：所有序列对同时计算
- ✅ **位点计数并行化**：同义/非同义位点统计
- ✅ **多方法并行比较**：YN00、NG86、LWL85同时运行
- ✅ **距离矩阵并行化**：大规模矩阵计算

### 4. TreeSub并行化
- ✅ **P矩阵指数运算**：矩阵指数的并行计算
- ✅ **距离矩阵计算**：K80、JC69距离
- ✅ **特征向量计算**：Q矩阵对角化
- ✅ **条件概率节点处理**：树上节点并行处理

## 📈 预期性能提升

| 线程数 | 预期加速比 | 适用场景 |
|--------|------------|----------|
| 2线程  | 1.5-1.8x  | 双核CPU |
| 4线程  | 2.5-3.5x  | 四核CPU |
| 8线程  | 4-6x      | 八核CPU |
| 16线程 | 6-10x     | 多核服务器 |

*注：实际性能提升取决于数据规模和硬件配置*

## 🔍 调试和故障排除

### 1. 编译错误
```bash
# 检查OpenMP支持
gcc -fopenmp --version

# 检查依赖文件
ls -la paml.h tools.c
```

### 2. 运行时错误
```bash
# 启用详细输出
export OMP_DISPLAY_ENV=true
export OMP_PROC_BIND=spread

# 检查线程创建
./codeml_parra test.ctl 2>&1 | grep -i thread
```

### 3. 性能问题
```bash
# 限制线程数
export OMP_NUM_THREADS=4

# 设置线程亲和性
export OMP_PROC_BIND=close

# 监控CPU使用率
htop & ./codeml_parra test.ctl
```

## 🔄 从原版迁移

### 1. 备份原文件
```bash
cp codeml.c codeml_original.c
cp baseml.c baseml_original.c  
cp yn00.c yn00_original.c
```

### 2. 验证结果一致性
```bash
# 运行原版
./codeml_original test.ctl > original.out

# 运行并行版（单线程）
export OMP_NUM_THREADS=1
./codeml_parra test.ctl > parallel.out

# 比较结果
diff original.out parallel.out
```

## 📚 技术细节

### OpenMP指令使用
- `#pragma omp parallel for`：循环并行化
- `reduction(+:variable)`：安全的累加操作
- `schedule(dynamic)`：动态负载均衡
- `critical`：关键区保护
- `sections`：任务并行化

### 内存管理
- 线程私有工作空间分配
- 避免false sharing
- 动态内存分配优化

### 数值稳定性
- 保持浮点运算顺序一致性
- 避免竞争条件影响结果
- 确保随机数生成的可重现性

---

*创建时间：2025年8月10日*
*PAML并行版本 v1.0*
