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
    Index idx{0, 0, *this};
    for ( ; idx.YIndex < static_cast<int>(yValues_.size()); ++idx.YIndex) // bottom to top
    {
        const std::vector<double> &kUp = kHorizontal[idx.YIndex]; 
        const std::vector<double> &kLeftRight = kVertical[idx.YIndex]; 
        const std::vector<double> &kDown = kHorizontal[idx.YIndex+1]; 

        const double deltaY_d = yValues_[idx.YIndex] - controlVolYValues[idx.YIndex];
        const double deltaY_u = controlVolYValues[idx.YIndex+1] - yValues_[idx.YIndex];
        const double deltaY = deltaY_d + deltaY_u;
        for ( ; idx.XIndex < static_cast<int>(xValues_.size()); ++idx.XIndex) // left to right
        {
            const double ku = kUp[idx.XIndex];
            const double kd = kDown[idx.XIndex];
            const double kl = kLeftRight[idx.XIndex];
            const double kr = kLeftRight[idx.XIndex+1];

            const double deltaX_l = xValues_[idx.XIndex] - controlVolXValues[idx.XIndex];
            const double deltaX_r = controlVolXValues[idx.XIndex+1] - xValues_[idx.XIndex];
            const double deltaX = deltaX_l + deltaX_r;

            const auto arrIdx = idx.ToArrayIndex();
            leftCoeffs_[arrIdx] = kl*deltaY / deltaX_l;
            rightCoeffs_[arrIdx] = kr*deltaY / deltaX_r;
            upCoeffs_[arrIdx] = ku*deltaX / deltaY_u;
            downCoeffs_[arrIdx] = kd*deltaX / deltaY_d;
            deltaXdeltaY_[arrIdx] = deltaX * deltaY;
            sources_[arrIdx] = sources[idx.YIndex][idx.XIndex];
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
    return {0, 0, *this};
}

int Grid::NumXValues() const
{
    return xValues_.size();
}

int Grid::NumYValues() const
{
    return yValues_.size();
}


