#!/bin/bash
# build_final_parallel.sh - 编译最终版本的并行可执行文件

echo "=================================================="
echo "       编译PAML并行版本可执行文件"
echo "=================================================="
echo ""

# 编译选项
CFLAGS="-fopenmp -O3 -Wall -march=native"
LIBS="-lm"

# 检查OpenMP支持
echo "🔍 检查编译环境..."
echo "GCC版本: $(gcc --version | head -1)"
echo -n "OpenMP版本: "
echo | gcc -fopenmp -dM -E - | grep _OPENMP | cut -d' ' -f3
echo "CPU核心数: $(nproc)"
echo ""

# 编译函数
compile_program() {
    local source=$1
    local target=$2
    
    echo "编译 $source -> $target..."
    if gcc $CFLAGS -o $target $source $LIBS; then
        echo "✅ $target 编译成功"
        # 检查文件大小
        local size=$(ls -lh $target | awk '{print $5}')
        echo "   文件大小: $size"
    else
        echo "❌ $target 编译失败"
        return 1
    fi
    echo ""
}

echo "=== 开始编译并行版本可执行文件 ==="
echo ""

# 编译所有并行版本
compile_program "codeml_parra.c" "codeml_parra"
compile_program "baseml_parra.c" "baseml_parra" 
compile_program "yn00_parra.c" "yn00_parra"
compile_program "treesub_parra.c" "treesub_parra"

echo "=== 编译完成情况统计 ==="
echo ""

# 检查编译结果
success_count=0
total_count=4

programs=("codeml_parra" "baseml_parra" "yn00_parra" "treesub_parra")

for prog in "${programs[@]}"; do
    if [ -f "$prog" ]; then
        echo "✅ $prog - 编译成功"
        success_count=$((success_count + 1))
    else
        echo "❌ $prog - 编译失败"
    fi
done

echo ""
echo "编译成功: $success_count/$total_count"

if [ $success_count -eq $total_count ]; then
    echo "🎉 所有并行版本编译成功！"
    
    echo ""
    echo "=== 可执行文件信息 ==="
    ls -lh *_parra 2>/dev/null || echo "没有找到编译好的可执行文件"
    
    echo ""
    echo "=== 快速测试 ==="
    echo "运行codeml_parra测试..."
    if [ -f "codeml_parra" ]; then
        timeout 10s ./codeml_parra 2>/dev/null && echo "✅ codeml_parra运行正常" || echo "⚠️  codeml_parra运行测试超时或出错"
    fi
    
    echo ""
    echo "🚀 使用方法:"
    echo "   ./codeml_parra   # 运行CodeML并行版本"
    echo "   ./baseml_parra   # 运行BaseML并行版本"
    echo "   ./yn00_parra     # 运行YN00并行版本" 
    echo "   ./treesub_parra  # 运行TreeSub并行版本"
    
else
    echo "⚠️  部分程序编译失败，请检查错误信息"
fi

echo ""
echo "=================================================="
