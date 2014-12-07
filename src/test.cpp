#include "tree.h"
#include "eval.h"
#include "util.h"
#include <random>
#include <iostream>
#include <cassert>
#include "Eigen/SparseLU"

double dist_kernel(double x, double y) {
    double dist = std::fabs(x - y);
    if (dist < 1e-10) {
        return 1;
    }
    return 1.0 / std::fabs(x - y);
}

void test_dense_matrix() {
    auto vec = to_eigen(random_list(100));
    auto mat = dense_matrix(vec, dist_kernel);
    auto vec2 = to_eigen(random_list(100));
    auto res = solve_system(mat, vec2);
    assert(res.size() == 100);
}

void check_begin_end(const Tree& t, const Cell& c) {
    if (c.children[0] == -1) {
        return;
    }

    auto left = t.cells[c.children[0]];
    auto right = t.cells[c.children[1]];
    assert(c.begin_idx == left.begin_idx);
    assert(left.end_idx == right.begin_idx);
    assert(c.end_idx == right.end_idx);

    check_begin_end(t, left);
    check_begin_end(t, right);
}

void test_treebuild() {
    auto pts = random_list(1000);
    auto t = Tree(pts, 10, 300);
    check_begin_end(t, t.root());
}

void test_treebuild_maxdepth() {
    auto pts = random_list(1000);
    auto t = Tree(pts, 1, 3);
    int max_depth = 0;
    for (auto c: t.cells) {
        max_depth = std::max(c.depth, max_depth);
    }
    check_begin_end(t, t.root());
    assert(max_depth == 3);
}

void test_treebuild_empty_cells() {
    auto pts = random_list(6);
    int max_depth = 5;
    auto t = Tree(pts, -1, max_depth);
    int last_end = 0;
    for (auto c: t.cells) {
        if(c.depth < max_depth) {
            continue;
        }
        assert(c.begin_idx == last_end);
        last_end = c.end_idx;
    }
}

void test_treecode_matrix_init() {
    auto pts = random_list(1000);
    auto t = Tree(pts, -1, 9);
    auto tmat = TreecodeMatrix(dist_kernel, t, 3.0);
    for (int level = 0; level <= t.max_depth; level++) {
        assert(tmat.cells_by_level[level].size() == pow(2, level));
    }
    //Better ordering test.
}

void test_treecode_eval_P2P() {
    // Compare the treecode Matrix-vector product with just P2P with an exact
    // matrix vector product
    auto pts = random_list(1000);
    auto t = Tree(pts, -1, 9);
    auto tmat = TreecodeMatrix(dist_kernel, t, 10000.0);
    tmat.P2P_M2P(); 
    tmat.build_matrix();

    auto vec = to_eigen(random_list(1000));
    Eigen::VectorXd est = tmat.matrix.block(0, 0, 1000, 1000) * vec;

    auto densemat = dense_matrix(to_eigen(t.elements), dist_kernel);
    Eigen::VectorXd exact = densemat * vec;
    for (int i = 0; i < 1000; i++) {
        assert(std::fabs((exact(i) - est(i)) / exact(i)) < 1e-2);
    }
}

void test_treecode_eval() {
    int n_pts = 1000;
    auto pts = random_list(n_pts);
    auto t = Tree(pts, -1, 9);
    auto tmat = TreecodeMatrix(dist_kernel, t, 5.0);
    tmat.P2P_M2P(); 
    tmat.P2M_M2M();
    tmat.build_matrix();

    auto vec = to_eigen(random_list(n_pts));
    auto est = tmat.eval(vec);

    auto densemat = dense_matrix(to_eigen(t.elements), dist_kernel);
    Eigen::VectorXd exact = densemat * vec;

    for (int i = 0; i < n_pts; i++) {
        assert(std::fabs((exact(i) - est(i)) / exact(i)) < 1e-2);
    }
}

void test_treecode_invert() {
    int n_pts = 200;
    auto pts = random_list(n_pts);
    auto t = Tree(pts, -1, 6);
    auto tmat = TreecodeMatrix(dist_kernel, t, 20.0);
    tmat.P2P_M2P(); 
    tmat.P2M_M2M();
    tmat.build_matrix();

    auto lu = banded_lu(tmat.matrix);
    auto mat = Eigen::MatrixXd(tmat.matrix);

    Eigen::VectorXd vec2 = to_eigen(random_list(n_pts));
    auto densemat = dense_matrix(to_eigen(t.elements), dist_kernel);
    Eigen::VectorXd exact = solve_system(densemat, vec2);

    Eigen::VectorXd extended_vec2(vec2.size() + t.cells.size());
    Eigen::VectorXd some_zeros = Eigen::VectorXd::Zero(t.cells.size());
    extended_vec2 << vec2, some_zeros;
    Eigen::VectorXd est = mat.partialPivLu().solve(extended_vec2);
    auto est2 = solve_lu(lu, extended_vec2);

    double error = (exact - est.segment(0, n_pts)).norm() / exact.norm();
    double error2 = (exact - est2.segment(0, n_pts)).norm() / exact.norm();
    std::cout << error << std::endl;
    std::cout << error2 << std::endl;
    assert(error < 1e-1);
}


int main() {
    test_dense_matrix();
    test_treebuild();
    test_treebuild_maxdepth();
    test_treebuild_empty_cells();
    test_treecode_matrix_init();
    test_treecode_eval_P2P();
    test_treecode_eval();
    test_treecode_invert();
    std::cout << "All tests passed." << std::endl;
}
