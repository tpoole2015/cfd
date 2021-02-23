#pragma once
#include <vector>
#include "grid.h"

struct Solution
{
    Solution(const Grid &grid);

    double operator()(const Grid::Index &idx) const;
    double& operator()(const Grid::Index &idx);
private:
    std::vector<double> values_;
};


