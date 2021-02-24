#include "solution.h"

Solution::Solution(const Grid &grid)
    : values_(grid.GetInitialValues())
{}

double Solution::operator()(const Grid::Index &idx) const
{
    return values_[idx.ToArrayIndex()];
}

double& Solution::operator()(const Grid::Index &idx)
{
    return values_[idx.ToArrayIndex()];
}


