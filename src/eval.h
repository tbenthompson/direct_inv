#ifndef __QKHQWLHJLHKA_EVAL_H
#define __QKHQWLHJLHKA_EVAL_H
#include <functional>
#include <vector>
#include "Eigen/Dense"
#include "Eigen/SparseCore"

typedef std::function<double(double,double)> Kernel;

Eigen::MatrixXd dense_matrix(const Eigen::VectorXd& locs, Kernel kernel);
Eigen::VectorXd solve_system(const Eigen::MatrixXd& mat, const Eigen::VectorXd& b);

class Tree;

Eigen::SparseMatrix<double> build_P2P(const Tree& t);
Eigen::SparseMatrix<double> build_M2P(const Tree& t);
Eigen::SparseMatrix<double> build_P2M(const Tree& t);
Eigen::SparseMatrix<double> build_M2M(const Tree& t);


#endif
