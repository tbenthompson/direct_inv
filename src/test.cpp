#include "tree.h"
#include "eval.h"
#include <random>
#include <iostream>
#include <cassert>

std::vector<double> random_list(int N) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0, 1);
    std::vector<double> es(N);
    for (int i = 0; i < N; ++i) {
        es[i] = dis(gen);
    }
    return es;
}

void test_dense_matrix() {
    auto vec = random_list(10);
    auto mat = dense_matrix(vec, [](double x, double y) {return 1 / std::fabs(x - y);});
}

void test_remove_if() {
    auto vec = random_list(10000); 
    auto out = include_if(vec, [](double x) {return x < 0.5;});
    // Very low probability of failing
    assert(out.size() < 7500);
}

void test_treebuild() {
    auto pts = random_list(1000);
    auto t = Tree(pts, 1);
    assert(t.left->right->elements.size() < 350);
}

int main() {
    test_dense_matrix();
    test_remove_if();
    test_treebuild();
}
