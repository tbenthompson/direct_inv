#include <iostream>
#include "eval.h"
#include "tree.h"

using Eigen::MatrixXd; 
using Eigen::VectorXd; 
using Eigen::SparseMatrix;

MatrixXd dense_matrix(const VectorXd& locs, Kernel kernel)
{
    int n_locs = locs.size();
    MatrixXd mat(n_locs, n_locs);
    for (int i = 0; i < n_locs; i++) {
        for (int j = 0; j < n_locs; j++) {
            if (i == j) {
                mat(i,j) = 0;
            } else {
                mat(i,j) = kernel(locs[i], locs[j]); 
            }
        }
    }
    return mat;
}

VectorXd solve_system(const Eigen::MatrixXd& A, const Eigen::VectorXd& b) {
    VectorXd res = A.partialPivLu().solve(b);
    return res;
}

SparseMatrix<double> build_P2P(const Tree& t) {
    /* SparseMatrix<double> p2p( */
}
