#include "eval.h"

std::vector<double> dense_matrix(const std::vector<double>& locs, Kernel kernel)
{
    int n_locs = locs.size();
    std::vector<double> matrix(n_locs * n_locs);
    for (int i = 0; i < n_locs; i++) {
        for (int j = 0; j < n_locs; j++) {
            int entry = i * n_locs + j;
            if (i == j) {
                matrix[entry] = 0;
            } else {
                matrix[entry] = kernel(locs[i], locs[j]); 
            }
        }
    }
    return matrix;
}

std::vector<double> solve_system(const std::vector<double>& A,
                                 const std::vector<double>& b) {

}
