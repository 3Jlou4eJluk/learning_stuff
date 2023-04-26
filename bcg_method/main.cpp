#include <iostream>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <chrono>

// Matrix is expected to implement the matvec function
// matvec(vec_to_mul, vec_to_save)
template <class T, class Matrix>
void bicg(const Matrix &A, const T *b, T threshold, T *x0=NULL, size_t n=3,
                    size_t num_threads=4, size_t max_iter=10) {
    T* tmp_vec1 = (T*)malloc(sizeof(T) * n);
    T* tmp_vec2 = (T*)malloc(sizeof(T) * n);
    T* tmp_vec3 = (T*)malloc(sizeof(T) * n);

    T* r = (T*)malloc(sizeof(T) * n);
    T* p = (T*)malloc(sizeof(T) * n);
    T* z = (T*)malloc(sizeof(T) * n);
    T* s = (T*)malloc(sizeof(T) * n);
    T* x = (T*)malloc(sizeof(T) * n);
    if (x0 != NULL) {
        #pragma omp parallel for num_threads(num_threads)
        for (size_t i = 0; i < n; i++) {
            x[i] = x0[i];
        }
    }

    // computing r0 = b - Ax0
    A.matvec(x, tmp_vec1);

    #pragma omp parallel for num_threads(num_threads)
    for (size_t i = 0; i < n; i++) {
        r[i] = b[i] - tmp_vec1[i];
    }
    // now tmp_vec2 stores r0

    #pragma omp parallel for num_threads(num_threads)
    for (size_t i = 0; i < n; i++) {
        p[i] = r[i];
        z[i] = r[i];
        s[i] = r[i];
    }

    T alpha_current;
    T beta_current;
    T tmp;

    T* x_current = (T*)malloc(sizeof(T) * n);
    T* r_current = (T*)malloc(sizeof(T) * n);
    T* p_current = (T*)malloc(sizeof(T) * n);
    T* z_current = (T*)malloc(sizeof(T) * n);
    T* s_current = (T*)malloc(sizeof(T) * n);

    T r_norm = 0;

    for (size_t i = 0; i < max_iter; i++) {
        // computing alpha_k
        alpha_current = 0;

        r_norm = 0;
        #pragma omp parallel for num_threads(num_threads)
        for (size_t i = 0; i < n; i++) {
            alpha_current += p[i] * r[i];
            r_norm += r[i] * r[i];
        }
        r_norm = std::sqrt(r_norm);
        if (r_norm < threshold) {
            std::cout << "BICG converged. Euclidean norm is " << r_norm << std::endl;
            return;
        }

        // computing Az_{k-1}
        A.matvec(z, tmp_vec1);

        tmp = 0;
        #pragma omp parallel for num_threads(num_threads)
        for (size_t i = 0; i < n; i++) {
            tmp += s[i] * tmp_vec1[i];
        }
        alpha_current /= tmp;

        // computing x_k
        #pragma omp parallel for num_threads(num_threads)
        for (size_t i = 0; i < n; i++) {
            x_current[i] = x[i] + z[i] * alpha_current;
        }

        // computing r_k
        #pragma omp parallel for num_threads(num_threads)
        for (size_t i = 0; i < n; i++) {
            r_current[i] = r[i] + tmp_vec1[i] * alpha_current * (-1);
        }

        // computing p_k
        A.tmatvec(s, tmp_vec1);
        #pragma omp parallel for num_threads(num_threads)
        for (size_t i = 0; i < n; i++) {
            p_current[i] = p[i] + tmp_vec1[i] * alpha_current * (-1);
        }

        // computing beta_k
        beta_current = 0;
        tmp = 0;
        #pragma omp parallel for num_threads(num_threads)
        for (size_t i = 0; i < n; i++) {
            beta_current += p_current[i] * r_current[i];
            tmp += p[i] * r[i];
        }
        beta_current /= tmp;

        // computing z_k
        #pragma omp parallel for num_threads(num_threads)
        for (size_t i = 0; i < n; i++) {
            z_current[i] = r_current[i] + z[i] * beta_current;
        }

        #pragma omp parallel for num_threads(num_threads)
        for (size_t i = 0; i < n; i++) {
            s_current[i] = p_current[i] + s[i] * beta_current;
        }

        #pragma omp parallel for num_threads(num_threads)
        for (size_t i = 0; i < n; i++) {
            x[i] = x_current[i];
            r[i] = r_current[i];
            p[i] = p_current[i];
            z[i] = z_current[i];
            s[i] = s_current[i];
        }
    }
    #pragma omp parallel for num_threads(num_threads)
    for (size_t i = 0; i < n; i++) {
        x0[i] = x_current[i];
    }

    std::cout << r_norm << std::endl;

    free(tmp_vec1);
    free(tmp_vec2);
    free(tmp_vec3);

    free(r);
    free(p);
    free(z);
    free(s);
    free(x);

    free(x_current);
    free(r_current);
    free(p_current);
    free(z_current);
    free(s_current);
    return;

}

template <typename T>
class Matrix {
public:
    size_t rows_;
    size_t cols_;
    size_t num_threads_;
    T* data_;

    Matrix(size_t rows, size_t cols, size_t num_threads, std::ifstream &in);
    ~Matrix();

    void print();
    void matvec(T* vec, T* save_vec) const;
    void tmatvec(T* vec, T* save_vec) const;
};


template <typename T>
Matrix<T>::Matrix(size_t rows, size_t cols, size_t num_threads, std::ifstream &in)
    : rows_(rows)
    , cols_(cols)
    , data_(nullptr)
    , num_threads_(num_threads)

{
    data_ = (T*)malloc(sizeof(T) * rows * cols);
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            in >> data_[i * cols + j];
        }
    }
}

template <typename T>
Matrix<T>::~Matrix() {
    free(data_);
}

template <typename T>
void Matrix<T>::print() {
    for (size_t i = 0; i < rows_; i++) {
        for (size_t j = 0; j < cols_; j++) {
            std::cout << data_[i][j] << " ";
        }
        std::cout << std::endl;
    }
}


template <typename T>
void Matrix<T>::matvec(T* vec, T* save_vec) const {
    #pragma omp for parallel num_threads(num_threads_)
    for (size_t i = 0; i < rows_; i++) {
        T mul_sum = 0;
        for (size_t j = 0; j < cols_; j++) {
            mul_sum += data_[i* cols_ + j] * vec[j];
        }
        save_vec[i] = mul_sum;
    }
}


template <typename T>
void Matrix<T>::tmatvec(T* vec, T* save_vec) const {
    #pragma omp for parallel num_threads(num_threads_)
    for (size_t i = 0; i < rows_; i++) {
        T mul_sum = 0;
        for (size_t j = 0; j < cols_; j++) {
            mul_sum += data_[j * rows_ + i] * vec[j];
        }
        save_vec[i] = mul_sum;
    }
}


int main(int argc, char *argv[]) {
    if (argc < 5) {
        std::cout << "Not enough parameters" << std::endl;
        return 1;
    }

    std::ifstream mat_in;
    std::ifstream vec_in;

    mat_in.open(argv[1]);
    vec_in.open(argv[2]);
    size_t n = atoi(argv[3]);
    size_t max_iter = static_cast<size_t>(atoi(argv[4]));

    auto b = static_cast<double*>(malloc(n * sizeof(double)));
    auto x0 = static_cast<double*>(malloc(n * sizeof(double)));

    Matrix<double> A(n, n, 4, mat_in);

    for (size_t i = 0; i < n; i++) {
        vec_in >> b[i];
        x0[i] = 0;
    }

    auto start = std::chrono::steady_clock::now();
    bicg(A, b, 1e-3, x0, n, 6, max_iter);
    auto end = std::chrono::steady_clock::now();

    std::chrono::duration<double> elapsed_seconds = end - start;

    std::cout << elapsed_seconds.count() << std::endl;

    for (size_t i = 0; i < n; i++) {
        std::cout << x0[i] << std::endl;
    }

    return 0;
}