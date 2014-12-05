#ifndef __QWELJLAHHSLHA_TREE_H
#define __QWELJLAHHSLHA_TREE_H
#include <vector>
#include <algorithm>
#include <iostream>
#include <array>

struct Cell {
    int depth;

    double min_pos;
    double max_pos;
    double width;

    int begin_idx;
    int end_idx;

    std::array<int,2> children;
};

struct Tree {
    Tree(const std::vector<double>& p_elements, 
         int max_per_cell,
         int max_depth):
        elements(p_elements),
        max_per_cell(max_per_cell),
        max_depth(max_depth)
    {
        std::sort(elements.begin(), elements.end());
        double min_pos = elements[0];
        double max_pos = elements[elements.size() - 1];
        Cell root = {0, min_pos, max_pos, max_pos - min_pos,
                     0, (int)elements.size(), {-1, -1}};
        if ((int)elements.size() <= max_per_cell || max_depth > 0) {
            root.children = build_children(root);
        }
        cells.push_back(root);
    }

    std::array<int,2> build_children(const Cell& cell) {
        return { build_child(cell, 0), build_child(cell, 1) };
    }

    int build_child(const Cell& cell, int i) {
        double min_pos = cell.min_pos + i * (cell.width / 2);
        double max_pos = cell.min_pos + (i + 1) * (cell.width / 2);

        int begin_idx = cell.begin_idx;
        int end_idx = cell.end_idx;
        if (i == 0) {
            for (int j = begin_idx; j < end_idx; j++) {
                if (elements[j] > max_pos) {
                    end_idx = j;
                    break;
                }
            }
        } else {
            for (int j = begin_idx; j < end_idx; j++) {
                if (elements[j] > min_pos) {
                    begin_idx = j;
                    break;
                }
            }
        }

        Cell child = {cell.depth + 1, min_pos, max_pos, (max_pos - min_pos), 
                      begin_idx, end_idx, {-1, -1}}; 
        if (end_idx - begin_idx > max_per_cell && max_depth > child.depth) {
            child.children = build_children(child);
        }
        cells.push_back(child);
        return cells.size() - 1;
    }

    Cell root() {
        return cells[cells.size() - 1];
    }

    std::vector<double> elements;
    std::vector<Cell> cells;
    int max_per_cell;
    int max_depth; 
};

#endif
