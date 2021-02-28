#pragma once
#include <vector>
#include <string>
#include "grid.h"

struct Solution
{
    Solution(const Grid &grid);

    double operator()(const Grid::Index &idx) const;
    double& operator()(const Grid::Index &idx);

    void WriteHeatMapToFile(const std::string &file) const;
private:
    const Grid* gridPtr_;
    std::vector<double> values_;
};


