#ifndef __QKHQWLHJLHKA_EVAL_H
#define __QKHQWLHJLHKA_EVAL_H
#include <functional>
#include <vector>

typedef std::function<double(double,double)> Kernel;

std::vector<double> dense_matrix(const std::vector<double>& locs, Kernel kernel);

#endif
