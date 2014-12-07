#include <iostream>
#include "eval.h"
#include "tree.h"
#include "Eigen/SparseLU"

using Eigen::MatrixXd; 
using Eigen::VectorXd; 
using Eigen::SparseMatrix;

MatrixXd dense_matrix(const VectorXd& locs, Kernel kernel)
{
    int n_locs = locs.size();
    MatrixXd mat(n_locs, n_locs);
    for (int i = 0; i < n_locs; i++) {
        for (int j = 0; j < n_locs; j++) {
            mat(i,j) = kernel(locs[i], locs[j]); 
        }
    }
    return mat;
}

VectorXd solve_system(const Eigen::MatrixXd& A, const Eigen::VectorXd& b) {
    VectorXd res = A.partialPivLu().solve(b);
    return res;
}

TreecodeMatrix::TreecodeMatrix(Kernel k, const Tree& t, double mac):
    kernel(k), tree(t), mac(mac),
    cell_to_dof(t.cells.size()),
    cells_by_level(t.max_depth + 1),
    level_start_dof(t.max_depth + 1)
{
    int n_pts = tree.elements.size();
    int n_cells = tree.cells.size();
    int n_dofs = n_pts + n_cells;
    matrix.resize(n_dofs, n_dofs);

    // Global ordering of treecode multipole dofs.
    for (std::size_t i = 0; i < tree.cells.size(); i++) {
        const auto c = tree.cells[i];
        cells_by_level[c.depth].push_back(i);
    }

    int next_dof = n_pts;
    for (int level = tree.max_depth; level >= 0; level--) {
        level_start_dof[level] = next_dof;
        for (int c_idx: cells_by_level[level]) {
            cell_to_dof[c_idx] = next_dof;
            next_dof++;
        }
    }
}

void TreecodeMatrix::P2P_M2P() {
    int n_pts = tree.elements.size();
    for (int i = 0; i < n_pts; i++) {
        P2P_M2P_cell(i, tree.root_index());
    }
}

void TreecodeMatrix::P2P_M2P_cell(int pt_idx, int cell_idx) {
    auto cell = tree.cells[cell_idx];
    auto pt = tree.elements[pt_idx];
    // Compare distance to MAC -- multipole acceptance criteria
    if (std::fabs((cell.center - pt) / cell.width) > mac) {
        // farfield, well-separated, so use M2P
        int cell_dof = cell_to_dof[cell_idx];
        triplets.push_back(Triplet(pt_idx, cell_dof, kernel(cell.center, pt)));
    } else if (cell.depth == tree.max_depth) {
        // nearfield, but cell is leaf so eval P2P
        for (int src_idx = cell.begin_idx; src_idx < cell.end_idx; src_idx++) {
            auto src = tree.elements[src_idx];
            triplets.push_back(Triplet(pt_idx, src_idx, kernel(src, pt)));
        }
    } else {
        // too close, but not leaf, so recurse
        P2P_M2P_cell(pt_idx, cell.children[0]);
        P2P_M2P_cell(pt_idx, cell.children[1]);
    }
}

void TreecodeMatrix::P2M_M2M() {
    int max_level = tree.max_depth;
    //P2M
    for (std::size_t c = 0; c < cells_by_level[max_level].size(); c++) {
        int cell_idx = cells_by_level[max_level][c];
        int cell_dof = cell_to_dof[cell_idx];
        auto cell = tree.cells[cell_idx];
        triplets.push_back(Triplet(cell_dof, cell_dof, -1.0));
        for (int src_idx = cell.begin_idx; src_idx < cell.end_idx; src_idx++) {
            triplets.push_back(Triplet(cell_dof, src_idx, 1.0));
        }
    }

    for (int level = max_level - 1; level >= 0; level--) {
        for (std::size_t c = 0; c < cells_by_level[level].size(); c++) {
            int cell_idx = cells_by_level[level][c];
            int cell_dof = cell_to_dof[cell_idx];
            auto cell = tree.cells[cell_idx];
            triplets.push_back(Triplet(cell_dof, cell_dof, -1.0));
            int child_dof0 = cell_to_dof[cell.children[0]];
            int child_dof1 = cell_to_dof[cell.children[1]];
            triplets.push_back(Triplet(cell_dof, child_dof0, 1.0));
            triplets.push_back(Triplet(cell_dof, child_dof1, 1.0));
        }
    }
}

void TreecodeMatrix::build_matrix() {
    matrix.setFromTriplets(triplets.begin(), triplets.end());
}

VectorXd TreecodeMatrix::eval(const VectorXd& vals) {
    int n_pts = tree.elements.size();
    int n_cells = tree.cells.size();
    SpMat P2P_block = matrix.topLeftCorner(n_pts, n_pts);
    SpMat M2P_block = matrix.topRightCorner(n_pts, n_cells);
    SpMat P2M_block = matrix.bottomLeftCorner(n_cells, n_pts);
    SpMat M2M_block = matrix.bottomRightCorner(n_cells, n_cells);

    VectorXd leaf_multipoles = P2M_block * vals;
    VectorXd all_multipoles = M2M_block.triangularView<Eigen::Lower>()
                                   .solve(-leaf_multipoles);

    VectorXd out = M2P_block * all_multipoles + P2P_block * vals;
    return out;
}

void TreecodeMatrix::invert() {
    int n_pts = tree.elements.size();
    int n_cells = tree.cells.size();
    SpMat P2P_block = matrix.topLeftCorner(n_pts, n_pts);
    SpMat M2P_block = matrix.topRightCorner(n_pts, n_cells);
    SpMat P2M_block = matrix.bottomLeftCorner(n_cells, n_pts);
    SpMat M2M_block = matrix.bottomRightCorner(n_cells, n_cells);

    Eigen::SparseLU<SpMat> solver;
    solver.compute(P2P_block);
    SparseMatrix<double> I(n_pts, n_pts);
    I.setIdentity();
    // auto P2P_inv = solver.solve(I);
    // auto P2P_inv2 = banded_lu(P2P_block); 
}

std::pair<SpMat,VectorXd> banded_lu(const SpMat& A_in) {
    //in place LU decomposition, but copies input parameter
    assert(A_in.rows() == A_in.cols());
    SpMat A(A_in);

    VectorXd diag(A.rows());
    std::vector<Triplet> out_entries;
    for (int i = 0; i < A.rows() - 1; i++) {
        // Break up matrix like:
        // [a_11 A_12] = [0    0][u_11 U_12]
        // [A_21 A_22]   [L_21 0][0    U_22]

        // u_11 = a_11 -- nothing happens
        diag(i) = A.coeff(i, i);
        std::cout << diag(i) << std::endl;
        std::cout << A.coeff(i + 1, i + 1) << std::endl;

        // U_12 = A_12 -- nothing happens

        // L_21 = (1/a_11) * A_21
        double inv_a11 = 1.0 / diag(i);
        for (int row = i + 1; row < A.rows(); row++) {
            A.coeffRef(row, i) *= inv_a11;
        }

        // A_22 - L_21 * U_12 = L_22*U_22
        // Subtract schur complement
        for (SparseMatrix<double>::InnerIterator it_col(A, i);
             it_col; ++it_col) {
            if (it_col.col() > i) {
                for (int row = i + 1; row < A.rows(); row++) {
                    for (SparseMatrix<double>::InnerIterator it_row(A, row); 
                         it_row; ++it_row) {
                        if (it_row.col() == i) {
                            A.coeffRef(it_row.row(), it_col.col()) -= 
                                it_row.value() * it_col.value(); 
                        }
                    }
                }
            }
        }
    }

    return {A, diag};
}

VectorXd solve_lu(const std::pair<SpMat,VectorXd>& lu, const VectorXd& b) {
    VectorXd y = lu.first.triangularView<Eigen::Lower>().solve(b);
    y = y.cwiseQuotient(lu.second);
    VectorXd soln = lu.first.triangularView<Eigen::Upper>().solve(y);
    return soln;
}

// InvTreecode::InvTreecode(SpMat matrix, int n_pts):
//     n_pts(n_pts),
//     n_cells(matrix.rows() - n_pts),
//     P2P_block(matrix.topLeftCorner(n_pts, n_cells)),
//     M2P_block(matrix.topRightCorner(n_pts, n_cells)),
//     P2M_block(matrix.bottomLeftCorner(n_cells, n_pts)),
//     M2M_block(matrix.bottomRightCorner(n_cells, n_cells))
// {
//     
// }
