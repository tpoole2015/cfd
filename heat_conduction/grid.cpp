#include <cassert>
#include "grid.h"
#include "solution.h"

Grid::Grid(const std::vector<double> &controlVolXValues, const std::vector<double> &controlVolYValues) :
    controlVolXValues_(controlVolXValues),
    controlVolYValues_(controlVolYValues)
{
    // Take the grid points to be the centers of the control volumes
    assert(controlVolXValues.size() > 1);
    xValues_.reserve(controlVolXValues.size() - 1);
    std::sort(controlVolXValues_.begin(), controlVolXValues_.end());  // sort left to right
    for (int i = 0; i < static_cast<int>(controlVolXValues.size()) - 1; ++i)
    {
        xValues_[i] = .5*(controlVolXValues[i] + controlVolXValues[i+1]);
    }

    assert(controlVolYValues.size() > 1);
    yValues_.reserve(controlVolYValues.size() - 1);
    std::sort(controlVolYValues_.begin(), controlVolYValues_.end());  // sort bottom to top 
    for (int i = 0; i < static_cast<int>(controlVolYValues.size()) - 1; ++i)
    {
        yValues_[i] = .5*(controlVolYValues[i] + controlVolYValues[i+1]);
    }
}

Grid::Point Grid::IndexToPoint(const Index &i) const
{
    return {xValues_[i.XIndex], yValues_[i.YIndex]};
}

Grid::Coefficients Grid::GetDiscretizationCoeffs(const Index &idx, const Solution &soln) const
{
    (void) idx;
    (void) soln;
    Coefficients c;
    return c;
}

Grid::Index Grid::GetBottomLeft() const
{
    return {0, NumYValues() - 1};
}

bool Grid::HasLeft(const Index &idx) const
{
    return idx.XIndex > 0;
}

bool Grid::HasRight(const Index &idx) const
{
    return idx.XIndex < NumXValues() - 1;
}

bool Grid::HasUp(const Index &idx) const
{
    return idx.YIndex > 0;
}

bool Grid::HasDown(const Index &idx) const
{
    return idx.YIndex < NumYValues() - 1;
}

Grid::Index Grid::MoveLeft(const Index &idx) const
{
    return {idx.XIndex-1, idx.YIndex};
}

Grid::Index Grid::MoveRight(const Index &idx) const
{
    return {idx.XIndex+1, idx.YIndex};
}

Grid::Index Grid::MoveUp(const Index &idx) const
{
    return {idx.XIndex, idx.YIndex-1};
}

Grid::Index Grid::MoveDown(const Index &idx) const
{
    return {idx.XIndex, idx.YIndex+1};
}

int Grid::NumXValues() const
{
    return static_cast<int>(xValues_.size());
}

int Grid::NumYValues() const
{
    return static_cast<int>(yValues_.size());
}


