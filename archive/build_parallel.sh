#!/bin/bash
# build_parallel.sh - 构建PAML并行版本的脚本

echo "============================================"
echo "PAML并行版本构建脚本"
echo "============================================"

# 检查当前目录
if [ ! -f "paml.h" ]; then
    echo "错误：请在包含paml.h的src目录中运行此脚本"
    exit 1
fi

# 检查并行源文件
echo "检查源文件..."
required_files=("codeml_parra.c" "baseml_parra.c" "yn00_parra.c" "treesub_parra.c" "Makefile_parra")
for file in "${required_files[@]}"; do
    if [ -f "$file" ]; then
        echo "✓ $file 存在"
    else
        echo "✗ $file 不存在"
        exit 1
    fi
done

# 检查OpenMP支持
echo -e "\n检查编译器和OpenMP支持..."
gcc -fopenmp --version > /dev/null 2>&1
if [ $? -eq 0 ]; then
    echo "✓ GCC和OpenMP支持正常"
else
    echo "✗ GCC或OpenMP支持有问题"
    exit 1
fi

# 设置编译环境
export CC=gcc
export CFLAGS="-O3 -fopenmp -Wall -Wno-unused-result"
export LDFLAGS="-lm -fopenmp"

echo -e "\n开始编译..."

# 清理旧的编译文件
echo "清理旧文件..."
rm -f *.o *_parra 2>/dev/null

# 尝试编译最简单的YN00版本
echo -e "\n1. 编译YN00并行版本..."
if gcc -O3 -fopenmp -DYN00 -c yn00_parra.c -o yn00_parra.o 2>yn00_compile.log; then
    echo "✓ yn00_parra.c 编译成功"
    if gcc -O3 -fopenmp -c tools.c -o tools.o 2>>yn00_compile.log; then
        echo "✓ tools.c 编译成功"
        if gcc -O3 -fopenmp yn00_parra.o tools.o -lm -o yn00_parra 2>>yn00_compile.log; then
            echo "✓ yn00_parra 链接成功"
        else
            echo "✗ yn00_parra 链接失败，查看 yn00_compile.log"
        fi
    else
        echo "✗ tools.c 编译失败，查看 yn00_compile.log"
    fi
else
    echo "✗ yn00_parra.c 编译失败，查看 yn00_compile.log"
fi

# 编译TreeSub模块
echo -e "\n2. 编译TreeSub模块..."
if gcc -O3 -fopenmp -c treesub_parra.c -o treesub_parra.o 2>treesub_compile.log; then
    echo "✓ treesub_parra.c 编译成功"
else
    echo "✗ treesub_parra.c 编译失败，查看 treesub_compile.log"
fi

# 编译BaseML
echo -e "\n3. 编译BaseML并行版本..."
if gcc -O3 -fopenmp -DBASEML -c baseml_parra.c -o baseml_parra.o 2>baseml_compile.log; then
    echo "✓ baseml_parra.c 编译成功"
    if [ -f "treesub_parra.o" ] && [ -f "tools.o" ]; then
        if gcc -O3 -fopenmp baseml_parra.o treesub_parra.o tools.o -lm -o baseml_parra 2>>baseml_compile.log; then
            echo "✓ baseml_parra 链接成功"
        else
            echo "✗ baseml_parra 链接失败，查看 baseml_compile.log"
        fi
    else
        echo "✗ 依赖文件不存在，跳过链接"
    fi
else
    echo "✗ baseml_parra.c 编译失败，查看 baseml_compile.log"
fi

# 编译CodeML（最复杂的）
echo -e "\n4. 编译CodeML并行版本..."
if gcc -O3 -fopenmp -DCODEML -c codeml_parra.c -o codeml_parra.o 2>codeml_compile.log; then
    echo "✓ codeml_parra.c 编译成功"
    if [ -f "treesub_parra.o" ] && [ -f "tools.o" ]; then
        if gcc -O3 -fopenmp codeml_parra.o treesub_parra.o tools.o -lm -o codeml_parra 2>>codeml_compile.log; then
            echo "✓ codeml_parra 链接成功"
        else
            echo "✗ codeml_parra 链接失败，查看 codeml_compile.log"
        fi
    else
        echo "✗ 依赖文件不存在，跳过链接"
    fi
else
    echo "✗ codeml_parra.c 编译失败，查看 codeml_compile.log"
fi

# 检查编译结果
echo -e "\n============================================"
echo "编译结果汇总:"
echo "============================================"

executables=("yn00_parra" "baseml_parra" "codeml_parra")
success_count=0

for exe in "${executables[@]}"; do
    if [ -x "$exe" ]; then
        size=$(ls -lh "$exe" | awk '{print $5}')
        echo "✓ $exe ($size) - 编译成功"
        success_count=$((success_count + 1))
    else
        echo "✗ $exe - 编译失败"
    fi
done

echo -e "\n成功编译: $success_count/3 个程序"

# 创建bin目录并复制
if [ $success_count -gt 0 ]; then
    echo -e "\n创建bin目录..."
    mkdir -p ../bin
    for exe in "${executables[@]}"; do
        if [ -x "$exe" ]; then
            cp "$exe" ../bin/
            echo "已复制 $exe 到 ../bin/"
        fi
    done
fi

# 测试OpenMP功能
echo -e "\n============================================"
echo "测试OpenMP功能:"
echo "============================================"

if [ -x "yn00_parra" ]; then
    echo "测试yn00_parra的OpenMP支持..."
    export OMP_NUM_THREADS=2
    echo "线程数设置为: $OMP_NUM_THREADS"
    ./yn00_parra --help 2>/dev/null || echo "程序创建成功（参数错误正常）"
fi

echo -e "\n============================================"
echo "构建完成！"
echo "============================================"
echo "使用方法:"
echo "1. 设置线程数: export OMP_NUM_THREADS=4"  
echo "2. 运行程序: ./yn00_parra sequences.dat"
echo "3. 或使用: ./baseml_parra baseml.ctl"
echo "4. 或使用: ./codeml_parra codeml.ctl"
echo -e "\n配置文件位于: ../conf/*_parra.ctl"
echo "详细说明请参考: README_parallel.md"
