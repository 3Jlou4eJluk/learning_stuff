
#include<stdio.h>

#include<cblas.h>  // заголовочный файл C-интерфейса библиотеки BLAS

#define M 300

#define N 400

#define K 500

int main()

{ int i,j;

  float A[M*K],   // массив расположен в памяти одним непрерывным блоком

        B[K][N],  // этот тоже

        *C;       // и этот тоже будет непрерывным

  C=(float*)malloc(M*N*sizeof(float)); // ну вот, непрерывный ;)

  for (i=0;i<M;i++)     // инициализируем массив A

    for (j=0;j<K;j++)

      A[i*K+j]=3*j+2*i;

  for (i=0;i<K;i++)     // инициализируем массив B

    for (j=0;j<N;j++)

      B[i][j]=5*i+j;

  for (i=0;i<M*N;i++) C[i]=5; // инициализируем массив C (зачем?)

  // перемножаем!!!

  // нам нужно C=AB, поэтому сделаем beta=0

  cblas_sgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,

              M,N,K,1.0,A,K,&B[0][0],N,0.0,C,N);

  for (i=0;i<M;i++) // выводим результат на экран (наверно не влезет: 300*400)

  { for (j=0;j<N;j++) printf("%.2f ",C[i*N+j]);

    printf("\n");

  }

  free(C);

  return 0;

}