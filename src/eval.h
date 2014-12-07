#ifndef __QKHQWLHJLHKA_EVAL_H
#define __QKHQWLHJLHKA_EVAL_H
#include <functional>
#include <vector>
#include "Eigen/Dense"
#include "Eigen/SparseCore"
#include "tree.h"

typedef std::function<double(double,double)> Kernel;

Eigen::MatrixXd dense_matrix(const Eigen::VectorXd& locs, Kernel kernel);
Eigen::VectorXd solve_system(const Eigen::MatrixXd& mat, const Eigen::VectorXd& b);

typedef Eigen::SparseMatrix<double,Eigen::RowMajor> SpMat;
typedef Eigen::Triplet<double> Triplet;

std::pair<SpMat,Eigen::VectorXd> banded_lu(const SpMat& A_in);
Eigen::VectorXd solve_lu(const std::pair<SpMat,Eigen::VectorXd>& lu,
                         const Eigen::VectorXd& b);

struct InvTreecode {
    // InvTreecode(SpMat matrix, int n_pts);

    // int n_pts;
    // int n_cells;
    // SpMat P2P_block;
    // SpMat M2P_block;
    // SpMat P2M_block;
    // SpMat M2M_block;
};

struct TreecodeMatrix {
    TreecodeMatrix(Kernel k, const Tree& t, double mac);

    const Kernel kernel;
    const Tree tree;
    const double mac;

    std::vector<int> cell_to_dof;
    std::vector<std::vector<int>> cells_by_level;
    std::vector<int> level_start_dof;

    SpMat matrix;

    std::vector<Triplet> triplets;

    void P2P_M2P();
    void P2P_M2P_cell(int pt_idx, int cell_idx);
    void P2M_M2M();
    void build_matrix();
    Eigen::VectorXd eval(const Eigen::VectorXd& vals);
    void invert();
};

#endif
