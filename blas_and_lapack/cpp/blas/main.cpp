#include <iostream>
#include <cblas.h>
#include <iomanip>



int main() {
	int m, n, k;
	std::cin >> m >> n >> k;

	float* mat1 = static_cast<float*>(malloc(m * n * sizeof(float)));
	float* mat2 = static_cast<float*>(malloc(n * k * sizeof(float)));

	float* mat3 = static_cast<float*>(malloc(m * k * sizeof(float)));

	float tmp = 0;
	for (size_t i = 0; i < m; i++) {
		for (size_t j = 0; j < n; j++) {
			std::cin >> tmp;
			mat1[i * n + j] = tmp;
		}
	}

	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < k; j++) {
			std::cin >> tmp;
			mat2[i * k + j] = tmp;
		}
	}

	for (size_t i = 0; i < m; i++) {
		for (size_t j = 0; j < k; j++) {
			mat3[i * k + j] = 0;
		}
	}

	cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1.0, mat1, m, mat2, n, 0.0, mat3, m);
	for (size_t i = 0; i < m; i++) {
		for (size_t j = 0; j < k; j++) {
			std::cout << std::setw(5) << std::setprecision(3) << mat3[i * k + j];
		}
		std::cout << std::endl;
	}

	free(mat1);
	free(mat2);
	free(mat3);
	return 0;
}