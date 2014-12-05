#include "tree.h"
#include "eval.h"
#include "util.h"
#include <random>
#include <iostream>
#include <cassert>

void test_dense_matrix() {
    auto vec = to_eigen(random_list(100));
    auto mat = dense_matrix(vec, [](double x, double y) {return 1 / std::fabs(x - y);});
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

int main() {
    test_dense_matrix();
    test_treebuild();
    test_treebuild_maxdepth();
    std::cout << "All tests passed." << std::endl;
}
