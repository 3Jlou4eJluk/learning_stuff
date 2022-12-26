#include <iostream>
#include <iomanip>

extern "C" {
    extern int dgesdd_(char*,int*,int*, double*, int*, double*, double*, int*, double*, int*, double*, int*, int*, int*);
}

int main() {
    int m, n;
    std::cin >> m >> n;

    double* mat = static_cast<double*>(malloc(m * n * sizeof(double)));

    double tmp = 0;
    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < n; j++) {
            std::cin >> tmp;
            mat[i * n + j] = tmp;
        }
    }

    int sigma_size = (m < n ? m : n);

    double* sigma = static_cast<double*>(malloc(sigma_size * sizeof(double)));

    double* u = static_cast<double*>(malloc(m * m * sizeof(double)));

    double* vs = static_cast<double*>(malloc(n * n * sizeof(double)));

    char q[] = "A";

    int* iwork = static_cast<int*>(malloc(8 * sigma_size * sizeof(int)));


    int info = 0;

    double work_size;
    double * work = &work_size;
    int lwork = -1;

    int lda = m;

    dgesdd_(q, &m, &n, mat, &lda, sigma, u, &m, vs, &n, work, &lwork, iwork, &info);

    lwork = work_size;

    work = static_cast<double*>(malloc(lwork * sizeof(double)));

    dgesdd_(q, &m, &n, mat, &lda, sigma, u, &m, vs, &n, work, &lwork, iwork, &info);

    std::cout << "U is " << std::endl;

    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < m; j++) {
            std::cout << std::setw(10) << std::setprecision(3) <<  u[i * m + j] <<  " ";
        }
        std::cout << std::endl;
    }

    std::cout << "Sigma is " << std::endl;
    for (size_t i = 0; i < sigma_size; i++) {
        std::cout << std::setw(10) << std::setprecision(3) << sigma[i] << " ";
    }
    std::cout << std::endl;

    std::cout << "V is " << std::endl;
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            std::cout << std::setw(10) << std::setprecision(3) << vs[i * n + j] << " ";
        }
        std::cout << std::endl;
    }
    return 0;
}