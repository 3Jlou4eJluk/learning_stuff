#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

// Подключение библиотеки INMOST, предполагая что она уже сконфигурирована и установлена
#include "inmost.h"

using namespace INMOST;

const int grid_step_x = 1;
const int grid_step_y = 1;

int main(int, char** argv) {
    Solver::Initialize(&argc, &argv);
    const int grid_points = 1000, total_nodes = grid_points * grid_points;
    const double mesh_size = 1.0 / grid_points;

    std::vector<double> reference_solution(total_nodes);
    for(int node_index = 0; node_index < total_nodes; node_index++) {
        int i = node_index / grid_points, j = node_index % grid_points;
        reference_solution[node_index] = std::sin(M_PI * i * mesh_size) * std::sin(M_PI * j * mesh_size) / (2 * M_PI * M_PI);
    }

    Sparse::Matrix system_matrix;
    Sparse::Vector right_hand_side;
    Sparse::Vector solution_vector;

    system_matrix.SetInterval(0, total_nodes);
    right_hand_side.SetInterval(0, total_nodes);
    solution_vector.SetInterval(0, total_nodes);

    for(int node_index = 0; node_index < total_nodes; ++node_index) {
        int i = node_index / grid_points, j = node_index % grid_points;
        double x = mesh_size * i, y = mesh_size * j;
        if(i * j == 0 || i == grid_points - 1 || j == grid_points - 1) {
            system_matrix[node_index][node_index] = 1.0;
            right_hand_side[node_index] = 0.0;
        } else {
            system_matrix[node_index][node_index] = 2.0 * (grid_step_x + grid_step_y);
            system_matrix[node_index][node_index + grid_points] = -grid_step_x;
            system_matrix[node_index][node_index - grid_points] = -grid_step_x;
            system_matrix[node_index][node_index - 1] = -grid_step_y;
            system_matrix[node_index][node_index + 1] = -grid_step_y;
            right_hand_side[node_index] = std::sin(M_PI * x) * std::sin(M_PI * y) * mesh_size * mesh_size;
        }
    }

    Solver approximate_solution(Solver::INNER_ILU2);
    approximate_solution.SetParameter("absolute_tolerance", "1e-12");
    approximate_solution.SetParameter("verbosity", "2");
    approximate_solution.SetParameter("relative_tolerance", "1e-7");
    approximate_solution.SetMatrix(system_matrix);
    bool solution_converged = approximate_solution.Solve(right_hand_side, solution_vector);

    // Вычисление нормы разности решений
    double norm_l2 = 0, norm_max = 0;
    for(int i = 0; i < total_nodes; ++i) {
        double diff = reference_solution[i] - solution_vector[i];
        norm_l2 += diff * diff;
        norm_max = std::max(norm_max, std::abs(diff));
    }
    norm_l2 = std::sqrt(norm_l2) / total_nodes;

    Solver::Finalize();

    std::cout << std::endl <<  "L2 Norm: " << norm_l2 << "\nMaximum Norm: " << norm_max << std::endl;

    return solution_converged ? 0 : 1; // Возвращаем статус завершения решения
}
