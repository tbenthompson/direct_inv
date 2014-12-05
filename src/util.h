#ifndef __JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ_UTIL_H
#define __JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ_UTIL_H
#include "Eigen/Dense"

inline std::vector<double> random_list(int N) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0, 1);
    std::vector<double> es(N);
    for (int i = 0; i < N; ++i) {
        es[i] = dis(gen);
    }
    return es;
}

inline Eigen::VectorXd to_eigen(const std::vector<double>& a) {
    return Eigen::VectorXd::Map(&a[0], a.size());
}

#endif
