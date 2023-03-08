#include <iostream>
#include <omp.h>
#include <stdlib.h>
#include <cmath>
#include "cblas.h"

#define size_type int32_t
#define THREAD_COUNT 3


double* create_vector(size_type vec_size) {
    return static_cast<double*>(std::malloc(vec_size * sizeof(double)));
}

void vector_add(const double* vec1, const double* vec2, double* res_vec, size_type vec_size) {
    #pragma omp parallel for num_threads(THREAD_COUNT)
    for (size_type i = 0; i < vec_size; i++) {
        res_vec[i] = vec1[i] + vec2[i];
    }
}

void vector_num_mul(const double* vec, double alpha, size_type n, double* res_vec) {
    #pragma omp parallel for num_threads(THREAD_COUNT)
    for (size_type i = 0; i < n; i++) {
        res_vec[i] = vec[i] * alpha;
    }
}

void copy_vector(double* to_vec, const double* from_vec, size_type vec_size) {
    #pragma omp parallel for num_threads(THREAD_COUNT)
    for (size_type i = 0; i < vec_size; i++) {
        to_vec[i] = from_vec[i];
    }

}


void mat_vec_mul(const double* matrix, const double* vec, double* res_vec, size_type n) {
    for (size_type i = 0; i < n; i++) {
        res_vec[i] = 0;
        #pragma omp parallel for num_threads(THREAD_COUNT)
        for (size_type j = 0; j < n; j++) {
            res_vec[i] += matrix[i * n + j] * vec[j];
        }
    }
}

void tmat_vec_mul(const double* matrix, const double* vec, double* res_vec, size_type n) {
    for (size_type i = 0; i < n; i++) {
        res_vec[i] = 0;
        #pragma omp parallel for num_threads(THREAD_COUNT)
        for (size_type j = 0; j < n; j++) {
            res_vec[i] += matrix[j * n + i] * vec[j];
        }
    }
}

double ddot(const double* vec1, const double* vec2, size_type n) {
    double res = 0;

    #pragma omp parallel for num_threads(THREAD_COUNT)
    for (size_type i = 0; i < n; i++) {
        res += vec1[i] * vec2[i];
    }

    return res;
}

/*
 * https://en.wikipedia.org/wiki/Biconjugate_gradient_method
 */
void bcg_method(double* matrix, double* b, double* x0, size_type n, size_type max_iter) {
    double* tmp_vec1 = create_vector(n);
    double* tmp_vec2 = create_vector(n);
    double* tmp_vec3 = create_vector(n);

    double* r = create_vector(n);
    double* p = create_vector(n);
    double* z = create_vector(n);
    double* s = create_vector(n);
    double* x = create_vector(n);
    copy_vector(x, x0, n);

    // computing r0 = b - Ax0
    mat_vec_mul(matrix, x0, tmp_vec1, n);
    vector_num_mul(tmp_vec1, -1, n, tmp_vec2);
    vector_add(b, tmp_vec2, r, n);

    // initializing p0, z0, s0 with r0
    copy_vector(p, r, n);
    copy_vector(z, r, n);
    copy_vector(s, r, n);

    // start iteration
    // cycle temporaries
    double alpha_current;
    double beta_current;

    double* x_current = create_vector(n);
    double* r_current = create_vector(n);
    double* p_current = create_vector(n);
    double* z_current = create_vector(n);
    double* s_current = create_vector(n);

    double tmp_double;

    for (size_type i = 0; i < max_iter; i++) {

        // computing alpha_k
        alpha_current = ddot(p, r, n);
        // computing Az_{k-1}
        mat_vec_mul(matrix, z, tmp_vec1, n);
        // now tmp_vec1 stores Az_{k-1}
        tmp_double = ddot(s, tmp_vec1, n);
        alpha_current /= tmp_double;

        // computing x_k
        vector_num_mul(z, alpha_current, n, tmp_vec2);
        vector_add(x, tmp_vec2, x_current, n);

        // computing r_k
        vector_num_mul(tmp_vec1, (-1) * alpha_current, n, tmp_vec3);
        vector_add(r, tmp_vec3, r_current, n);

        // computing p_k
        tmat_vec_mul(matrix, s, tmp_vec1, n);
        vector_num_mul(tmp_vec1, (-1) * alpha_current, n, tmp_vec2);
        vector_add(p, tmp_vec2, p_current, n);

        // computing beta_k
        beta_current = ddot(p_current, r_current, n);
        beta_current /= ddot(p, r, n);

        // computing z_k
        vector_num_mul(z, beta_current, n, tmp_vec1);
        vector_add(r_current, tmp_vec1, z_current, n);

        // computing s_k
        vector_num_mul(s, beta_current, n, tmp_vec1);
        vector_add(p_current, tmp_vec1, s_current, n);

        copy_vector(x, x_current, n);
        copy_vector(r, r_current, n);
        copy_vector(p, p_current, n);
        copy_vector(z, z_current, n);
        copy_vector(s, s_current, n);
    }

    copy_vector(x0, x_current, n);

    free(r);
    free(p); free(z);
    free(s); free(x);
    free(x_current); free(r_current);
    free(p_current); free(z_current);
    free(s_current); free(tmp_vec1);
    free(tmp_vec2); free(tmp_vec3);

}


int main() {

    auto A = static_cast<double*>(malloc(9 * sizeof(double)));
    auto b = static_cast<double*>(malloc(3 * sizeof(double)));
    auto x0 = static_cast<double*>(malloc(3 * sizeof(double)));

    for (size_type i = 0; i < 3; i++) {
        for (size_type j = 0; j < 3; j++) {
            std::cin >> A[i * 3 + j];
        }
    }

    for (size_type i = 0; i < 3; i++) {
        std::cin >> b[i];
        x0[i] = 0;
    }

    bcg_method(A, b, x0, 3, 25);

    for (size_type i = 0; i < 3; i++) {
        std::cout << x0[i] << std::endl;
    }

    return 0;
}
