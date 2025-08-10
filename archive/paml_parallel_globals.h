/* paml_parallel_globals.h - 并行版本需要的全局变量和函数定义 */
#ifndef PAML_PARALLEL_GLOBALS_H
#define PAML_PARALLEL_GLOBALS_H

#include "paml.h"
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// 全局变量定义（简化但兼容版本）
char BASEs[] = "TCAG";
char AAs[] = "ARNDCQEGHILKMFPSTWYV*";
int noisy = 1;
int NFunCall = 0;
int NEigenQ = 0; 
int NPMatUVRoot = 0;
int *ancestor = NULL;
int GeneticCode[1][64];
double *SeqDistance = NULL;
double SS = 0.0, NN = 0.0, Sd = 0.0, Nd = 0.0;

// 必要的工具函数实现
FILE* gfopen(const char* filename, const char* mode) {
    FILE* fp = fopen(filename, mode);
    if (!fp && noisy) {
        printf("Warning: Could not open file %s\n", filename);
    }
    return fp;
}

void error(const char* msg) {
    printf("Error: %s\n", msg);
    exit(1);
}

void zerror(const char* format, ...) {
    printf("Error: %s\n", format);
    exit(1);
}

double etime(void) {
    return omp_get_wtime();
}

// 简化但兼容的内存管理
void* smalloc(size_t size) {
    void* ptr = malloc(size);
    if (!ptr) error("Memory allocation failed");
    return ptr;
}

void sfree(void* ptr) {
    if (ptr) free(ptr);
}

// 数学工具函数
double square(double x) { return x * x; }
double max2(double a, double b) { return (a > b) ? a : b; }
double min2(double a, double b) { return (a < b) ? a : b; }

#endif /* PAML_PARALLEL_GLOBALS_H */
