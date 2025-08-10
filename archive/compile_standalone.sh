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
