#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include "mpi.h"

void init_matrix(double* matrix, int n, int first_row, int last_row, int max_rows, int my_rank, int p, const char *name){
    int i, j, k;
    double* a = matrix;
    double* buf = new double[max_rows * n];
    memset(buf, 0, max_rows*n*sizeof(double));
    double tmp = 0.;
    int first_row_k = 0, last_row_k = 0, err = 0;
    FILE* fp = nullptr;
    if(my_rank == 0){
        fp = fopen(name, "r");
        if(fp == nullptr){
            err = 1;
        }
    }
    MPI_Bcast(&err, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if(err) return;
    for(k = 0; k < p; ++k){
        if(my_rank == 0){
            if(k == 0){
                for(i = first_row; i <= last_row; ++i){
                    for(j = 0; j < n; ++j){
                        if((fscanf(fp, "%lf", &tmp) != 1)){
                            return;
                        }
                        *(a++) = tmp;
                    }
                }
            }
            else{
                memset(buf, 0, max_rows*n*sizeof(double));
                first_row_k = n * k;
                first_row_k /= p; 
                last_row_k = n * (k + 1);
                last_row_k = last_row_k / p - 1;
                for(i = first_row_k; i <= last_row_k; ++i){
                    for(j = 0; j < n; ++j){
                        if((fscanf(fp, "%lf", &tmp) != 1)){
                            return;
                        }
                        buf[(i - first_row_k)*n + j] = tmp;
                    }
                }
                MPI_Send(buf, max_rows * n, MPI_DOUBLE, k, 0, MPI_COMM_WORLD); 
            }

        }
        else{
            if(my_rank == k){
                MPI_Recv(buf, max_rows*n, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);            
                for(i = 0; i < last_row - first_row + 1; ++i){
                    for(j = 0; j < n; ++j){
                        *(a++) = buf[i*n + j];
                    }
                }
            }

        }
    }

    if(my_rank == 0){
        fclose(fp);
    }
    delete[] buf;
}

void init_vector(double* vector, int n, int first_row, int last_row, int max_rows, int my_rank, int p, const char *name){
    int i, k;
    double* a = vector;
    double* buf = new double[max_rows];
    memset(buf, 0, max_rows*sizeof(double));
    double tmp = 0.;
    int first_row_k = 0, last_row_k = 0, err = 0;
    FILE* fp = nullptr;
    if(my_rank == 0){
        fp = fopen(name, "r");
        if(fp == nullptr){
            err = 1;
        }
    }
    MPI_Bcast(&err, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if(err) return;
    for(k = 0; k < p; ++k){
        if(my_rank == 0){
            if(k == 0){
                for(i = first_row; i <= last_row; ++i){
                    if((fscanf(fp, "%lf", &tmp) != 1)){
                        return;
                    }
                    *(a++) = tmp;
                }
            }
            else{
                memset(buf, 0, max_rows*sizeof(double));
                first_row_k = n * k;
                first_row_k /= p; 
                last_row_k = n * (k + 1);
                last_row_k = last_row_k / p - 1;
                for(i = first_row_k; i <= last_row_k; ++i){
                    if((fscanf(fp, "%lf", &tmp) != 1)){
                        return;
                    }
                    buf[(i - first_row_k)] = tmp;
                }
                MPI_Send(buf, max_rows, MPI_DOUBLE, k, 0, MPI_COMM_WORLD); 
            }

        }
        else{
            if(my_rank == k){
                MPI_Recv(buf, max_rows, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);                
                for(i = 0; i < last_row - first_row + 1; ++i){
                    *(a++) = buf[i];
                }
            }

        }
    }

    if(my_rank == 0){
        fclose(fp);
    }
    delete[] buf;
}

double scalar_product(double *a, double *b, int first_row, int last_row){
    int rows = last_row - first_row + 1;
    int i = 0;
    double sum = 0;
    for(i = 0; i < rows; ++i){
        sum += a[i]*b[i];
    }
    return sum;
}

void matrix_mult_vector(double *a, double *b,double *c, int n, int first_row, int last_row, int my_rank, int p){
    int i, j, k, l;
    int rows = last_row - first_row+ 1;
    int max_rows;
    int first_row_k, last_row_k, rows_k;
    double s;
    int dest, source, tag = 0;
    MPI_Status status;

    max_rows = n / p + n % p;
    for(i = 0; i < rows; ++i){
        c[i] = 0.;
    }

    source = (my_rank + 1) % p;
    if(my_rank == 0){
        dest = p - 1;
    }
    else{
        dest = my_rank - 1;
    }

    for(l = 0; l < p; ++l){
        k = (my_rank + l) % p;
        first_row_k = n * k;
        first_row_k /= p;
        last_row_k = n * (k + 1);
        last_row_k = last_row_k / p - 1;
        rows_k = last_row_k - first_row_k + 1;

        for(i = 0; i < rows; ++i){
            for(s = 0., j = 0; j < rows_k; ++j){
                s += a[i * n + j +first_row_k] * b[j];
            }
            c[i] += s;
        }
        MPI_Sendrecv_replace(b, max_rows, MPI_DOUBLE, dest, tag, source, tag, MPI_COMM_WORLD, &status);
    }
}

void matrix_mult_vector_t(double *a, double *b,double *c, int n, int first_row, int last_row, int my_rank, int p){
    int i, j, k, l;
    int rows = last_row - first_row+ 1;
    int max_rows;
    int first_row_k, last_row_k, rows_k;
    double s;
    int dest, source, tag = 0;
    MPI_Status status;

    max_rows = n / p + n % p;
    for(i = 0; i < rows; ++i){
        c[i] = 0.;
    }

    source = (my_rank + 1) % p;
    if(my_rank == 0){
        dest = p - 1;
    }
    else{
        dest = my_rank - 1;
    }

    for(l = 0; l < p; ++l){
        k = (my_rank + l) % p;
        first_row_k = n * k;
        first_row_k /= p;
        last_row_k = n * (k + 1);
        last_row_k = last_row_k / p - 1;
        rows_k = last_row_k - first_row_k + 1;

        for(i = 0; i < rows; ++i){
            for(s = 0., j = 0; j < rows_k; ++j){
                s = a[i * n + j +first_row_k] * b[i];
                c[j] += s;
            }
        }
        MPI_Sendrecv_replace(c, max_rows, MPI_DOUBLE, dest, tag, source, tag, MPI_COMM_WORLD, &status);
    }
}


void substraction(double* a, double* b, int first_row, int last_row){
    int rows = last_row - first_row + 1;
    int i = 0;
    for(i = 0; i < rows; ++i){
        a[i] = a[i] - b[i];
    }
}

void sgb_method(double* matr, double* b, double* x, int max_iter, int first_row, int last_row, int max_rows, int n, int my_rank, int size, double eps){
    int rows = first_row - last_row + 1;
    rows = max_rows;
    int i = 0, j = 0, err = 0;
    double alpha = 0, beta = 0, sum = 0, res = 0, sum1 = 0, norm = 0, res1 = 0;
    double* r = new double[rows];    
    double* r_current = new double[rows];
    double* p = new double[rows];
    double* p_current = new double[rows];
    double* z = new double[rows];
    double* s = new double[rows];
    double* result = new double[rows];
    double* result1 = new double[rows];
    memset(r, 0, rows * sizeof(double));
    memset(p, 0, rows * sizeof(double));
    memset(z, 0, rows * sizeof(double));
    memset(s, 0, rows * sizeof(double));
    memset(result, 0, rows * sizeof(double));
    memset(r_current, 0, rows * sizeof(double));
    memset(p_current, 0, rows * sizeof(double));

    matrix_mult_vector(matr, x, result, n, first_row, last_row, my_rank, size);
    substraction(b, result, first_row, last_row);
    for(i = 0; i < rows; ++i){
        r[i] = b[i];
        p[i] = r[i];
        z[i] = r[i];
        s[i] = r[i];
    }
    
    for(j = 0; j < max_iter; ++j){
        sum = scalar_product(p, r, first_row, last_row);
        sum1 = scalar_product(r, r, first_row, last_row);
        MPI_Allreduce(&sum, &res, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&sum1, &norm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        if(my_rank == 0){
            norm = std::sqrt(norm);
            if(norm < eps){
                printf("BICG converged. Euclidean norm is %lf\n", norm);
                err = 1;
            }
        }
        MPI_Bcast(&err, 1, MPI_INT, 0, MPI_COMM_WORLD);
        if(err == 1){
            //delete all
            delete[] r;    
            delete[] r_current; 
            delete[] p; 
            delete[] p_current;
            delete[] z;
            delete[] s; 
            delete[] result; 
            delete[] result1; 
            return;
        }
        memset(result, 0, rows*sizeof(double));
        matrix_mult_vector(matr, z, result, n, first_row, last_row, my_rank, size);
        sum = scalar_product(s, result, first_row, last_row);
        MPI_Allreduce(&sum, &res1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        matrix_mult_vector_t(matr, s, result1, n, first_row, last_row, my_rank, size);
        alpha = res / res1;
        for(i = 0; i < rows; ++i){
            x[i] = x[i] + alpha * z[i];
            r_current[i] = r[i] - alpha * result[i];
            p_current[i] = p[i] - alpha * result1[i];
        }
        MPI_Barrier(MPI_COMM_WORLD);
        res = 0;
        res1 = 0;
        sum = scalar_product(p_current, r_current, first_row, last_row);
        MPI_Allreduce(&sum, &res, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        sum = scalar_product(p, r, first_row, last_row);
        MPI_Allreduce(&sum, &res1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        beta = res / res1;
        for(i = 0; i < rows; ++i){
            z[i] = r_current[i] + beta * z[i];
            s[i] = p_current[i] + beta * s[i];
        }
        for(i = 0; i < rows; ++i){
            r[i] = r_current[i];
            p[i] = p_current[i];
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    //delete all
    delete[] r;    
    delete[] r_current; 
    delete[] p; 
    delete[] p_current;
    delete[] z;
    delete[] s; 
    delete[] result; 
    delete[] result1; 
}

void print_matrix_full(double* matr, int first_row, int last_row, int my_rank, int size, int n, int max_rows){
    int i = 0, k = 0, j = 0;
    int rows = last_row - first_row + 1, rows1 = 0;
    double *buf = new double[max_rows * n];
    
    
    if(my_rank == 0){
        for(i = 0; i < rows; ++i){
            for(j = 0; j < n; ++j){
                printf(" %12.6lf", matr[i*n + j]);
            }
            printf("\n");
        }
    }

    for(i = 1; i < size; ++i){
        
        if(my_rank == 0){
            MPI_Recv(&rows1, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(buf, max_rows * n, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for(k = 0; k < rows1; ++k){
                for(j = 0; j < n; ++j){
                    printf(" %12.6lf", buf[k*n + j]);
                }
                printf("\n");
            }
        }
        else{
            if(my_rank == i){                
                MPI_Send(&rows, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
                MPI_Send(matr, max_rows * n, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            }
        }
    }
}

int main(int argc, char** argv){
    int my_rank = 0; // ранг текущего процесса
    int p = 0; // общее число процессов
    int n = 0; // размер матрицы
    int max_iter = 0;
    double t = 0; //время работы всей задачи
    double eps = 1e-8;
    int first_row = 0, last_row = 0;//, rows;
    int max_rows = 0;
    double *matrix = nullptr;
    double *vector = nullptr;
    double *result = nullptr; // результирующий вектор
    char *name = 0;
    char* name_vector = 0;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank); // номер текущего процесса
    MPI_Comm_size(MPI_COMM_WORLD, &p); // общее число процессов

    if(!((argc == 6) 
    && sscanf(argv[1], "%d", &n) == 1
    && sscanf(argv[2], "%d", &max_iter) == 1
    && sscanf(argv[3], "%lf", &eps) == 1)){
        if(my_rank == 0){
            printf("Usage %s n max_iter eps [file] [file]\n", argv[0]);
        }
        MPI_Finalize();
        return 0;
    }
    name = argv[4];
    name_vector = argv[5];

    first_row = n * my_rank;
    first_row /= p;

    last_row = n * (my_rank + 1);
    last_row = last_row / p - 1;

    max_rows = n / p + n % p;

    if(!(matrix = (double*)malloc(max_rows * n * sizeof(double)))){
        printf("Not enough memory!\n");
        MPI_Finalize();
        return 0;
    }

    if(!(vector = (double*)malloc(max_rows * sizeof(double)))){
        printf("Not enough memory!\n");
        MPI_Finalize();
        return 0;
    }

    if(!(result = (double*)malloc(max_rows * sizeof(double)))){
        printf("Not enough memory!\n");
        MPI_Finalize();
        return 0;
    }

    init_matrix(matrix, n, first_row, last_row, max_rows, my_rank, p, name);
    init_vector(vector, n, first_row, last_row, max_rows, my_rank, p, name_vector);
    if(my_rank == 0){
        printf("Matrix:\n");
    }
    if(n <= 20) print_matrix_full(matrix, first_row, last_row, my_rank, p, n, max_rows);
    if(my_rank == 0){
        printf("Vector:\n");
    }
    if(n <= 20) print_matrix_full(vector, first_row, last_row, my_rank, p, 1, max_rows); //vector size n = 1
    memset(result, 0, max_rows*sizeof(double));
    MPI_Barrier(MPI_COMM_WORLD);
    t = MPI_Wtime();
    sgb_method(matrix, vector, result, max_iter, first_row, last_row, max_rows, n, my_rank, p, eps);
    MPI_Barrier(MPI_COMM_WORLD);
    t = MPI_Wtime() - t;
    if(my_rank == 0){
        printf("Result:\n");
    }
    print_matrix_full(result, first_row, last_row, my_rank, p, 1, max_rows);
    if(my_rank == 0){
        printf("Total time = %le\n", t);
    }
    free(matrix);
    free(vector);
    free(result);

    MPI_Finalize();
    return 0;
}