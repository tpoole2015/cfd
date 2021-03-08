#pragma once
#include <vector>
#include <string>
#include "grid.h"

struct Solution
{
    Solution(const Grid *gridPtr);

    void Init(const std::vector<double> &values);

    double operator()(const Grid::Point &pt) const;
    double& operator()(const Grid::Point &pt);

    void WriteHeatMapToFile(const std::string &file) const;
private:
    const Grid* gridPtr_;
    std::vector<double> values_;
};


