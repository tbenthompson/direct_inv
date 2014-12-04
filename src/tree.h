#ifndef __QWELJLAHHSLHA_TREE_H
#define __QWELJLAHHSLHA_TREE_H
#include <vector>
#include <algorithm>
#include <iostream>

/* Return a vector including only those elements of the input vector
 * that pass the predicate.
 */
inline std::vector<double> include_if(const std::vector<double>& in,
                                      const std::function<bool(double)>& pred) {
    std::vector<double> out;
    for (auto val: in) {
        if (pred(val)) {
            out.push_back(val); 
        }
    }
    return out;
}

/* Inhomogeneous tree where all the elements are placed in the leaves
 * and the branches simply separate space at the interval level.
 */
struct Tree {
    Tree(const std::vector<double>& pts, std::size_t max_per_cell):
        min_pos(*std::min_element(pts.begin(), pts.end())),
        max_pos(*std::max_element(pts.begin(), pts.end())),
        split((min_pos + max_pos) / 2),
        elements(pts),
        left(subtree(pts, max_per_cell, [&](double x){return x <= split;})),
        right(subtree(pts, max_per_cell, [&](double x){return x > split;}))
    {}

    ~Tree() {
        delete left; 
        delete right;
    }

    static Tree* subtree(const std::vector<double>& pts,
                         std::size_t max_per_cell,
                         const std::function<bool(double)>& pred) {
        if (pts.size() <= max_per_cell) {
            return nullptr;
        }
        auto subtree_pts = include_if(pts, pred);
        return new Tree(subtree_pts, max_per_cell);
    }
    
    const double min_pos;
    const double max_pos;
    const double split;
    const std::vector<double> elements;
    const Tree* left;
    const Tree* right;
};

#endif
