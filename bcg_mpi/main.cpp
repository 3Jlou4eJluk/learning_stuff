#include<mpi.h>
#include "Header.h"

int main(int argc, char** argv){
    int my_rank, size, n, max_iter, m = 0, err = 0;
    double t1 =0, eps = 1e-3;
    char* filename, *file, *file3;

    double **A = nullptr, *b = nullptr, *x = nullptr, *buf1 = nullptr, **A_t = nullptr, *columns = nullptr, *b1 = nullptr, *b2 = nullptr;
    //MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank); //номер текущего процесса
    MPI_Comm_size(MPI_COMM_WORLD, &size); //всего процессов

    if (!((argc == 6)
		&& sscanf(argv[1], "%d", &n) == 1
		&& sscanf(argv[2], "%d", &max_iter) == 1)) {
			if(my_rank == 0){
				printf("Usage %s n max_iter [file]\n", argv[0]);
			}
		MPI_Finalize();
		return 0;
	}
	filename = argv[3];
    file = argv[4];
    file3 = argv[5];
    m = (n % size >= my_rank + 1) ? (n/size + 1) : (n/size); //количество строк у процесса
    A = CreateMatrix(m , n);
    A_t = CreateMatrix(m, n);
    b = new double[m];
    x = new double[m];
    b1 = new double[m];
    b2 = new double[m];
    buf1 = new double[n];
    columns = new double[m*n];
    memset(x, 0, m*sizeof(double));
    memset(columns, 0, m*n*sizeof(double));
    err = ReadMatrix(A, b2, n, filename, my_rank, size, buf1);
    
    if(err < 0){
        printf("Error read matrix!\n");
        DeleteMatrix(A, m);
        delete[] b;
        delete[] x;
        delete[] buf1;
        MPI_Finalize();
        return 0;
    }

    MPI_Barrier(MPI_COMM_WORLD);
    err = ReadMatrix(A_t, b1, n, file, my_rank, size, buf1);    
    if(err < 0){
        printf("Error read matrix!\n");
        DeleteMatrix(A, m);
        delete[] b;
        delete[] x;
        delete[] buf1;
        MPI_Finalize();
        return 0;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    Read_b(b, n, file3, my_rank, size);
    MPI_Barrier(MPI_COMM_WORLD);
    if(my_rank == 0) t1 = MPI_Wtime();
    bicg(A, A_t, b, x, n, m, eps, max_iter, my_rank, size); //A_t
    MPI_Barrier(MPI_COMM_WORLD);
    if(my_rank == 0) t1 = MPI_Wtime() - t1; 


    MPI_Barrier(MPI_COMM_WORLD);
    PrintMatrix(A, n, n, n, my_rank, size, buf1);
    MPI_Barrier(MPI_COMM_WORLD);
    
    if(my_rank == 0) printf("T  = %10.3e\n", t1);

    if(my_rank == 0){
        for (int i = 0; i < n; i++) {
            printf(" %10.3e", x[i]);
        }
    }

    if(my_rank == 0) printf("\n");

    
    DeleteMatrix(A, m);
    DeleteMatrix(A_t, m);
    delete[] b;
    delete[] x;
    delete[] buf1;
    MPI_Finalize();
    return 0;
}
