#include "solution.h"

Solution::Solution(const Grid &grid)
    : grid_(grid), values_(grid.NumXValues*grid.NumYValues)
{}

const Grid &Solution::GetGridReference() const
{
    return grid_;
}

double Solution::operator()(const Grid::Index &idx) const
{
    return values_[idx.YIndex*grid_.NumXValues + idx.XIndex];
}

double& Solution::operator()(const Grid::Index &idx)
{
    return values_[idx.YIndex*grid_.NumXValues + idx.XIndex];
}


