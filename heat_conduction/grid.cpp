#include <cassert>
#include "grid.h"
#include "solution.h"

Grid::Grid(std::vector<double> controlVolXValues, 
           std::vector<double> controlVolYValues, 
           std::vector<std::vector<double>> kVertical,
           std::vector<std::vector<double>> kHorizontal,
           std::vector<std::vector<double>> sources) :
    xValues_(controlVolXValues.size() - 1),
    yValues_(controlVolYValues.size() - 1)
{
    // we only do this once so it doesn't have to be fast
 
    // take the grid points to be the centers of the control volumes
    assert(controlVolXValues.size() > 1);
    std::sort(controlVolXValues.begin(), controlVolXValues.end());  // sort left to right
    for (int i = 0; i < static_cast<int>(controlVolXValues.size()) - 1; ++i)
    {
        xValues_[i] = .5*(controlVolXValues[i] + controlVolXValues[i+1]);
    }

    assert(controlVolYValues.size() > 1);
    std::sort(controlVolYValues.begin(), controlVolYValues.end());  // sort bottom to top 
    for (int i = 0; i < static_cast<int>(controlVolYValues.size()) - 1; ++i)
    {
        yValues_[i] = .5*(controlVolYValues[i] + controlVolYValues[i+1]);
    }

    // construct what we can of the discretization coefficients
    assert(kVertical.size() == controlVolYValues.size() - 1);
    assert(kVertical[0].size() == controlVolXValues.size());
    assert(kHorizontal.size() == controlVolYValues.size());
    assert(kHorizontal[0].size() == controlVolXValues.size() - 1);
    assert(sources.size() == yValues_.size());
    assert(sources[0].size() == xValues_.size());
    for (int i = 0; i < static_cast<int>(yValues_.size()); ++i) // bottom to top
    {
        const std::vector<double> &kUp = kHorizontal[i]; 
        const std::vector<double> &kLeftRight = kVertical[i]; 
        const std::vector<double> &kDown = kHorizontal[i+1]; 

        const double deltaY_d = yValues_[i] - controlVolYValues[i];
        const double deltaY_u = controlVolYValues[i+1] - yValues_[i];
        const double deltaY = deltaY_d + deltaY_u;
        for (int j = 0; j < static_cast<int>(xValues_.size()); ++j) // left to right
        {
            const double ku = kUp[j];
            const double kd = kDown[j];
            const double kl = kLeftRight[j];
            const double kr = kLeftRight[j+1];

            const double deltaX_l = xValues_[j] - controlVolXValues[j];
            const double deltaX_r = controlVolXValues[j+1] - xValues_[j];
            const double deltaX = deltaX_l + deltaX_r;

            leftCoeffs_[i*xValues_.size() + j] = kl*deltaY / deltaX_l;
            rightCoeffs_[i*xValues_.size() + j] = kr*deltaY / deltaX_r;
            upCoeffs_[i*xValues_.size() + j] = ku*deltaX / deltaY_u;
            downCoeffs_[i*xValues_.size() + j] = kd*deltaX / deltaY_d;
            deltaXdeltaY_[i*xValues_.size() + j] = deltaX * deltaY;
            sources_[i*xValues_.size() + j] = sources[i][j];
        }
    }
}

Grid::Point Grid::IndexToPoint(const Index &i) const
{
    return {xValues_[i.XIndex], yValues_[i.YIndex]};
}

Grid::Coefficients Grid::GetDiscretizationCoeffs(double alpha, double dt, const Index &idx, const Solution &soln) const
{
    const int n = idx.YIndex * NumXValues() + idx.XIndex;
    const double a = alpha*deltaXdeltaY_[n]/dt;
    Coefficients c;
    c.Left = leftCoeffs_[n];
    c.Right = rightCoeffs_[n];
    c.Up = upCoeffs_[n];
    c.Down = downCoeffs_[n];
    c.Constant = a*soln(idx);
    c.Center = c.Left + c.Right + c.Up + c.Down + a - sources_[n]*deltaXdeltaY_[n]; 
    return c;
}

Grid::Index Grid::GetOrigin() const
{
    return {0, 0};
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
    return idx.YIndex < NumYValues() - 1;
}

bool Grid::HasDown(const Index &idx) const
{
    return idx.YIndex > 0;
}

Grid::Index Grid::MoveLeft(const Index &idx) const
{
    return {idx.XIndex-1, idx.YIndex};
}

Grid::Index Grid::MoveRight(const Index &idx) const
{
    return {idx.XIndex+1, idx.YIndex};
}

Grid::Index Grid::MoveDown(const Index &idx) const
{
    return {idx.XIndex, idx.YIndex-1};
}

Grid::Index Grid::MoveUp(const Index &idx) const
{
    return {idx.XIndex, idx.YIndex+1};
}

int Grid::NumXValues() const
{
    return xValues_.size();
}

int Grid::NumYValues() const
{
    return yValues_.size();
}


