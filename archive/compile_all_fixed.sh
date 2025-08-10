#!/bin/bash
echo "=== 编译所有修复后的并行版本 ==="

CFLAGS="-fopenmp -O3 -Wall"

echo "编译 codeml_parra_fixed_v2..."
gcc $CFLAGS -o codeml_parra_fixed_v2 codeml_parra_fixed_v2.c -lm
if [ $? -eq 0 ]; then
    echo "✓ codeml_parra_fixed_v2 编译成功"
else
    echo "✗ codeml_parra_fixed_v2 编译失败"
fi

echo "编译 baseml_parra_fixed_v2..."
gcc $CFLAGS -o baseml_parra_fixed_v2 baseml_parra_fixed_v2.c -lm  
if [ $? -eq 0 ]; then
    echo "✓ baseml_parra_fixed_v2 编译成功"
else
    echo "✗ baseml_parra_fixed_v2 编译失败"
fi

echo "编译 yn00_parra_fixed_v2..."
gcc $CFLAGS -o yn00_parra_fixed_v2 yn00_parra_fixed_v2.c -lm
if [ $? -eq 0 ]; then
    echo "✓ yn00_parra_fixed_v2 编译成功"
else  
    echo "✗ yn00_parra_fixed_v2 编译失败"
fi

echo "编译 treesub_parra_fixed_v2..."
gcc $CFLAGS -o treesub_parra_fixed_v2 treesub_parra_fixed_v2.c -lm
if [ $? -eq 0 ]; then
    echo "✓ treesub_parra_fixed_v2 编译成功"
else
    echo "✗ treesub_parra_fixed_v2 编译失败"  
fi

echo "=== 编译完成 ==="
ls -la *_parra_fixed_v2
