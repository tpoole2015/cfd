#pragma once
#include <vector>

struct Grid;
struct Solution
{
    Solution(const Grid &grid);

    const Grid &GetGridReference() const;

    double operator()(const Grid::Index &idx) const;
    double& operator()(const Grid::Index &idx);
private:
    const Grid &grid_;
    std::vector<double> values_;
};


