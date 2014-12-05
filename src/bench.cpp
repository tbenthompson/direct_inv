#include <chrono>
#include "eval.h"
#include "tree.h"
#include "util.h"

/* Some simple timing functions. */
#define TIC\
    std::chrono::high_resolution_clock::time_point start =\
        std::chrono::high_resolution_clock::now();\
    int time_ms;
#define TIC2\
    start = std::chrono::high_resolution_clock::now();
#define TOC(name)\
    time_ms = std::chrono::duration_cast<std::chrono::milliseconds>\
                (std::chrono::high_resolution_clock::now() - start).count();\
    std::cout << name << " took "\
              << time_ms\
              << "ms.\n";


void bench_eigenlu() {
    int n = 2500;
    TIC
    auto vec = to_eigen(random_list(n));
    auto mat = dense_matrix(vec, [](double x, double y) {return 1 / std::fabs(x - y);});
    auto vec2 = to_eigen(random_list(n));
    auto res = solve_system(mat, vec2);
    TOC("Dense inverse on " + std::to_string(n) + " points")
}

void bench_treeconstruct() {
    TIC
    int n = (int)1e6;
    auto pts = random_list(n);
    auto t = Tree(pts, 10, 1000);
    TOC("Tree construct on " + std::to_string(n) + " points")
}

int main() {
    bench_eigenlu();
    bench_treeconstruct();
}
