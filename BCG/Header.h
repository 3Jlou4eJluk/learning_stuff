#pragma once
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <locale.h>
#include <time.h>
#include <math.h>

double** CreateMatrix(int n, int m);
void DeleteMatrix(double** A, int m);
int ReadMatrix(double** A, int n, char* filename, int my_rank, int size, double* buf);
int ReadMatrixColumns(char* filename, int n, double** columns, int rank, int size);
void PrintMatrix(double** A, int n, int l, int m, int my_rank, int size, double* buf1);
void bicg(double** A, double** A_t, double* b, double* x0, int n, int m, double eps, int max_iter, int my_rank, int size);
int Read_b(double*b, int n, char* filename, int rank, int size);
int FillMatrix(double** matr, double *b, int n, int k, int rank, int size);
int FillMatrixTranspon(double** matr, int n, int k, int rank, int size);