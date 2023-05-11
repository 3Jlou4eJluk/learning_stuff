#include "Header.h"

double** CreateMatrix(int n1, int n2){
    int i;
    double **m;
    m = (double**)malloc(n1*sizeof(double*));
    for (i = 0; i < n1; i++)
    {
        m[i] = (double*)malloc(n2*sizeof(double));
    }
    return m;
}


void DeleteMatrix(double** matr, int n1){
    int i;
    for (i = 0; i < n1; i++) {
        free(matr[i]);
    }
    free(matr);
}

int Read_b(double*b, int n, char* filename, int rank, int size){
    int d, i, j, c = 0;
    MPI_Status status;
    FILE* fin;
    if (rank == 0){
        fin = fopen(filename, "r");
        if (!fin){
            printf("Error opening file %s\n", filename);
            c = -1;
        }
    }
    MPI_Bcast(&c, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if(c) return -1;
    /*for (i = 0; i < n; ++i){
        buf1[i] = 0;
    }*/
    double a = 0;
    for (i = 0; i < n; i++){
        d = i % size;
        if (rank == 0){
            //for (j = 0; j < n; j++){
                if (fscanf(fin, "%lf", &a) != 1){
                    printf("Error reading matr[i][j]!");
                    c = -1;
                }
            //}
            if (d != 0){
                MPI_Send(&a, 1, MPI_DOUBLE, d, 0, MPI_COMM_WORLD);
            }
        }
        else if (rank == d){
            MPI_Recv(&a, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
        }
        if (rank == d){
            //b[i/size] = 0;
            for(j = 0; j < n; j++){
                if (j % 2 == 0){
                    b[i/size] = a;
                }
            }
        }
    }
    if (rank == 0) fclose (fin);
    MPI_Bcast(&c, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (c) return -1;
    return 0;
}

int ReadMatrix(double** matr, double *b, int n, char* filename, int rank, int size, double *buf1){
    double t;
    int d, i, j, c = 0;
    MPI_Status status;
    FILE* fin;
    if (rank == 0){
        fin = fopen(filename, "r");
        if (!fin){
            printf("Error opening file %s\n", filename);
            c = -1;
        }
    }
    MPI_Bcast(&c, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if(c) return -1;
    for (i = 0; i < n; ++i){
        buf1[i] = 0;
    }
    for (i = 0; i < n; i++){
        d = i % size;
        if (rank == 0){
            for (j = 0; j < n; j++){
                if (fscanf(fin, "%lf", &buf1[j]) != 1){
                    printf("Error reading matr[i][j]!");
                    c = -1;
                }
            }
            if (d != 0){
                MPI_Send(buf1, n, MPI_DOUBLE, d, 0, MPI_COMM_WORLD);
            }
        }
        else if (rank == d){
            MPI_Recv(buf1, n, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
        }
        if (rank == d){
            b[i/size] = 0;
            for(j = 0; j < n; j++){
                t = buf1[j];
                matr[i/size][j] = t;
                if (j % 2 == 0){
                    b[i/size] = b[i/size] + t;
                }
            }
        }
    }
    if (rank == 0) fclose (fin);
    MPI_Bcast(&c, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (c) return -1;
    return 0;
}

int ReadMatrixColumns(char* filename, int n, double* columns, int rank, int size){
    double* buf;
    FILE *f;
    int l;
    if (rank == 0){
        // Проверка файла
        f = fopen(filename,"r");
        if(!f){
            printf("Can't open file\n");
            return -1;
        }
    
    }

    // Если всё хорошо, то делаем буфер для строки
    buf = (double*)malloc(n*sizeof(double));
    for (int i=0;i<n;i++){
        if (rank==0){
            for (int j=0;j<n;j++)
                fscanf(f, "%lf", &buf[j]);
        }
        // Отправляем строку
        MPI_Bcast(buf, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        // Каждый достает нужное из строки
        l=0;
        for (int k=rank;k<n;k+=size){
            columns[n*l+i] = buf[k];
            l++;
        }
    }

    MPI_Bcast(buf, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    if (rank == 0)
        fclose(f);
    free(buf);
    return 0;
}

void PrintMatrix(double **matr, int n, int l, int m, int rank, int size, double *buf1){
    int i, j;
    double t;
    MPI_Status status;
    if (m == -1) {
        m = n;
    }
    if ((m > 0) || (m == -1)){
        for (i = 0; i < fmin (m, l); i++){
            if ((rank == i%size) && (rank != 0)) {
                for (j = 0; j < fmin (n, m); j++){
                    t = matr[i/size][j];
                    buf1[j] = t;
                }
                MPI_Send(buf1, m, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            }
            if ((rank == i%size) && (rank == 0)){
                for (j = 0; j < fmin (n, m); j++){
                    printf(" %10.3e", matr[i/size][j]);
                }
                fputc('\n', stdout);
            }
            if ((rank != i%size) && (rank == 0)){
                MPI_Recv(buf1, m, MPI_DOUBLE, i%size, 0, MPI_COMM_WORLD, &status);
                for (j = 0; j < fmin (n, m); j++){
                    printf(" %10.3e", buf1[j]);
                }
                fputc('\n', stdout);
            }
        }
    }
}


void mat_vec(double** A, double* x, double* tmp_vec1, int n, int my_rank, int size){
    int i = 0, j = 0;
    double sum = 0;
    for(i = 0; i < n; ++i){
        if(my_rank == i%size){
            sum = 0;
            for(j = 0; j < n; ++j){
                sum = sum + A[i/size][j] * x[j];
            }
            tmp_vec1[i] = sum;
        }
    }
}

void mat_vec_T(double* A, double* x, double* tmp_vec1, int n, int my_rank, int size){
    int m = (n % size >= my_rank + 1) ? (n/size + 1) : (n/size); //количество строк у процесса
    int i = 0, j = 0;
    double sum = 0;
    for(i = 0; i < m; ++i){
        sum = 0;
        for(j = 0; j < n; ++j){
            sum = sum + A[i*n + j] * x[j];
        }
        tmp_vec1[i] = sum;
        
    }
}


void bicg(double** A, double** A_t, double* b, double* x0, int n, int /*m*/, double eps, int max_iter, int my_rank, int size){
   double* tmp_vec1 = nullptr, *tmp_vec2 = nullptr, *x = nullptr;
    double* r = nullptr, *p = nullptr, *z = nullptr, *s = nullptr, *z1 = nullptr, *s1 = nullptr;
    double* x_current = nullptr, *r_current = nullptr, *p_current = nullptr, *z_current = nullptr, *s_current = nullptr;
    double alpha_current = 0, beta_current = 0, r_norm = 0;
    double alpha_current_1 = 0, beta_current_1 = 0, r_norm_1 = 0;
    double tmp = 0, tmp1 = 0;
    int i = 0, j = 0, err = 0;

    tmp_vec1 = new double[n];
    tmp_vec2 = new double[n];
    x = new double[n];
    r = new double[n];
    p = new double[n];
    z = new double[n];
    s = new double[n];
    z1 = new double[n];
    s1 = new double[n];
    x_current = new double[n];
    r_current = new double[n];
    p_current = new double[n];
    z_current = new double[n];
    s_current = new double[n];
    
    memset(r, 0, n*sizeof(double));
    memset(p, 0, n*sizeof(double));
    memset(z, 0, n*sizeof(double));
    memset(s, 0, n*sizeof(double));
    memset(z1, 0, n*sizeof(double));
    memset(s1, 0, n*sizeof(double));
    
    memset(x_current, 0, n*sizeof(double));
    memset(r_current, 0, n*sizeof(double));
    memset(p_current, 0, n*sizeof(double));
    memset(z_current, 0, n*sizeof(double));
    memset(s_current, 0, n*sizeof(double));
    if(my_rank == 0){
        if (x0 != NULL) {
            for (i = 0; i < n; i++) {
                x[i] = x0[i];
            }
        }
    }
    MPI_Bcast(x, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if(my_rank == 0){
        tmp_vec1[0] = 0;
    }    
    MPI_Bcast(tmp_vec1, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    //computing tmp_vec2 =  A*x0
    mat_vec(A, x, tmp_vec1, n, my_rank, size);
    MPI_Allreduce(tmp_vec1, tmp_vec2, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    //computing r  = b - tmp_vec2
    //fill p, z, s = r
    for(i = 0; i < n; ++i){
        if(i%size == my_rank){
            r[i] = b[i/size] - tmp_vec2[i];
            p[i] = r[i];
            z1[i] = r[i];
            s1[i] = r[i];
        }
    }
    MPI_Allreduce(z1, z, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(s1, s, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    for(j = 0; j < max_iter; ++j){
        alpha_current = 0;
        alpha_current_1 = 0;
        r_norm = 0;
        r_norm_1 = 0;
        MPI_Barrier(MPI_COMM_WORLD);
        for (i = 0; i < n; i++){
            if(i%size == my_rank){
                alpha_current_1 += p[i] * r[i];
                r_norm_1 += r[i] * r[i];
            }
        }
        //MPI_Barrier(MPI_COMM_WORLD);
        MPI_Allreduce(&alpha_current_1, &alpha_current, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        //MPI_Barrier(MPI_COMM_WORLD);
        MPI_Allreduce(&r_norm_1, &r_norm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        //MPI_Barrier(MPI_COMM_WORLD);
        if(my_rank == 0){
            r_norm = std::sqrt(r_norm);
            if(r_norm < eps){
                printf("BICG converged. Euclidean norm is %lf\n", r_norm);
                err = 1;
            }
            
        }
        //MPI_Barrier(MPI_COMM_WORLD);
        MPI_Bcast(&err, 1, MPI_INT, 0, MPI_COMM_WORLD);
        //MPI_Barrier(MPI_COMM_WORLD);                             //podumat
        if(err == 1){
            memset(x, 0, n*sizeof(double));
            for (i = 0; i < n; i++) {
                if(i%size == my_rank){
                    x[i] = x_current[i];
                }
            }
            MPI_Allreduce(x, x0, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            //ВСЕ УДАЛИТЬ!!!!!!!!!
            delete[] tmp_vec1; 
            delete[] tmp_vec2; 
            delete[] r; 
            delete[] p;
            delete[] z; 
            delete[] s; 
            delete[] z1; 
            delete[] s1; 
            delete[] x_current; 
            delete[] r_current;
            delete[] p_current;
            delete[] z_current;
            delete[] s_current;    
            MPI_Barrier(MPI_COMM_WORLD);
			return;
		}
        
        //computing A*z(k-1)
        memset(tmp_vec1, 0, n*sizeof(double));
        memset(tmp_vec2, 0, n*sizeof(double));
        mat_vec(A, z, tmp_vec2, n, my_rank, size);
        MPI_Allreduce(tmp_vec2, tmp_vec1, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        //MPI_Barrier(MPI_COMM_WORLD);
        tmp = 0, tmp1 = 0;
        for (i = 0; i < n; i++){
            if(i%size == my_rank){
                tmp1 += s[i] * tmp_vec1[i];
            }
        }
        //MPI_Barrier(MPI_COMM_WORLD);
        MPI_Allreduce(&tmp1, &tmp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        alpha_current /= tmp;

        MPI_Barrier(MPI_COMM_WORLD); //&нужно ли?

         // computing x_k
        for (i = 0; i < n; i++) {
            if(i%size == my_rank){
                x_current[i] = x[i] + z[i] * alpha_current;
            }
        }
        //MPI_Barrier(MPI_COMM_WORLD);
         // computing r_k
        for (i = 0; i < n; i++) {
            if(i%size == my_rank){
                r_current[i] = r[i] + tmp_vec1[i] * alpha_current * (-1);
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
        // computing p_k
        memset(tmp_vec1, 0, n*sizeof(double));
        memset(tmp_vec2, 0, n*sizeof(double));

        mat_vec(A_t, s, tmp_vec2, n, my_rank, size);//A_t        
        MPI_Allreduce(tmp_vec2, tmp_vec1, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        for (i = 0; i < n; i++) {
            if(i%size == my_rank){
                p_current[i] = p[i] + tmp_vec1[i] * alpha_current * (-1);
            }
        }

        //MPI_Barrier(MPI_COMM_WORLD);
        // computing beta_k
        beta_current = 0, beta_current_1 = 0;
        tmp = 0, tmp1 = 0;
        for (i = 0; i < n; i++) {
            if(i%size == my_rank){
                beta_current_1 += p_current[i] * r_current[i];
                tmp1 += p[i] * r[i];
            }
        }
        //MPI_Barrier(MPI_COMM_WORLD);
        MPI_Allreduce(&tmp1, &tmp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        //MPI_Barrier(MPI_COMM_WORLD);
        MPI_Allreduce(&beta_current_1, &beta_current, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        //MPI_Barrier(MPI_COMM_WORLD);
        beta_current /= tmp;

        MPI_Barrier(MPI_COMM_WORLD);
         // computing z_k
        for (i = 0; i < n; i++) {
            if(i % size == my_rank){
                z_current[i] = r_current[i] + z[i] * beta_current;            
            }
        }
        //MPI_Barrier(MPI_COMM_WORLD);
        for (i = 0; i < n; i++) {
            if(i%size == my_rank){
                s_current[i] = p_current[i] + s[i] * beta_current;
            }
        }
        //MPI_Barrier(MPI_COMM_WORLD);
        memset(z1, 0, n*sizeof(double));
        memset(s1, 0, n*sizeof(double));
        for (i = 0; i < n; i++) {
            if(i%size == my_rank){
                x[i] = x_current[i];
                r[i] = r_current[i];
                p[i] = p_current[i];
                z1[i] = z_current[i];
                s1[i] = s_current[i];
            }
        }
        MPI_Allreduce(s1, s, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(z1, z, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);  
    
    for (i = 0; i < n; i++) {
        if(i%size == my_rank){
            x[i] = x_current[i];
        }
    }
    MPI_Allreduce(x, x0, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    if(my_rank == 0) printf("norm = %10.3e\n", r_norm);
    MPI_Barrier(MPI_COMM_WORLD);  
    delete[] tmp_vec1; 
    delete[] tmp_vec2; 
    delete[] x;
    delete[] r; 
    delete[] p;
    delete[] z; 
    delete[] s; 
    delete[] z1; 
    delete[] s1; 
    delete[] x_current; 
    delete[] r_current;
    delete[] p_current;
    delete[] z_current;
    delete[] s_current;
}