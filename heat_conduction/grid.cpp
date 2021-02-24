#include <cassert>
#include "grid.h"
#include "solution.h"

namespace
{
    void ConstructCenter(const InputVariables::ControlVolume &vol,
                         double &xValue,
                         double &yValue,
                         double &leftCoeff,
                         double &rightCoeff,
                         double &upCoeff,
                         double &downCoeff,
                         double &deltaXdeltaY)
    {
        xValue = 0.5*(vol.BottomRight.X + vol.TopLeft.X);
        yValue = 0.5*(vol.BottomRight.Y + vol.TopLeft.Y);

        const double deltaY_d = yValue - vol.BottomRight.Y;
        const double deltaY_u = vol.TopLeft.Y - yValue;
        const double deltaY = deltaY_d + deltaY_u;
        const double deltaX_l = xValue - vol.TopLeft.X;
        const double deltaX_r = vol.BottomRight.X - xValue;
        const double deltaX = deltaX_l + deltaX_r;

        const double kl = vol.K[0];
        const double kr = vol.K[1];
        const double ku = vol.K[2];
        const double kd = vol.K[3];

        leftCoeff = kl*deltaY / deltaX_l;
        rightCoeff = kr*deltaY / deltaX_r;
        upCoeff = ku*deltaX / deltaY_u;
        downCoeff = kd*deltaX / deltaY_d;
        deltaXdeltaY = deltaX * deltaY;
    }

    void ConstructBoundary(const InputVariables::ControlVolume &vol,
                           double &xValue,
                           double &yValue,
                           double &constantCoeff,
                           double &deltaXdeltaY)
    {
        // TODO: handle flux boundary conditions
        assert(vol.BVType == InputVariables::ControlVolume::Scalar && vol.BoundaryValues.size() == 1);
        xValue = 0.5*(vol.BottomRight.X + vol.TopLeft.X);
        yValue = 0.5*(vol.BottomRight.Y + vol.TopLeft.Y);

        const double deltaX = vol.BottomRight.X - vol.TopLeft.X;
        const double deltaY = vol.BottomRight.Y - vol.TopLeft.Y;

        constantCoeff = vol.BoundaryValues[0];
        deltaXdeltaY = deltaX * deltaY;
    }
}

Grid::Grid(InputVariables input) :
    NumXValues(input.NumXVols),
    NumYValues(input.NumYVols),
    xValues_(input.NumXVols*input.NumYVols, 0),
    yValues_(input.NumXVols*input.NumYVols, 0),
    leftCoeffs_(input.NumXVols*input.NumYVols, 0),
    rightCoeffs_(input.NumXVols*input.NumYVols, 0),
    upCoeffs_(input.NumXVols*input.NumYVols, 0),
    downCoeffs_(input.NumXVols*input.NumYVols, 0),
    constantCoeffs_(input.NumXVols*input.NumYVols, 0),
    deltaXdeltaY_(input.NumXVols*input.NumYVols, 0),
    sources_(input.NumXVols*input.NumYVols, 0),
    initialValues_(input.NumXVols*input.NumYVols, 0)
{
    // we only do this once so it doesn't have to be fast
    
    assert(static_cast<int>(input.Volumes.size()) == input.NumXVols*input.NumYVols);
    auto cmp = [](const InputVariables::ControlVolume &lhs, const InputVariables::ControlVolume &rhs) -> bool {
        // sort bottom to top then left to right
        if (lhs.TopLeft.Y == rhs.TopLeft.Y)
        {
            return lhs.TopLeft.X < rhs.TopLeft.X;
        } 
        return lhs.TopLeft.Y < rhs.TopLeft.Y;
    };
    std::sort(input.Volumes.begin(), input.Volumes.end(), cmp);

    double xValue = 0;
    double yValue = 0;
    double leftCoeff = 0;
    double rightCoeff = 0;
    double upCoeff = 0;
    double downCoeff = 0;
    double constantCoeff = 0;
    double deltaXdeltaY = 0;

    Index idx{0, 0, *this};
    for ( ; idx.YIndex < input.NumYVols; ++idx.YIndex)
    {
        for ( ; idx.XIndex < input.NumXVols; ++idx.XIndex)
        {
            const auto arrIdx = idx.ToArrayIndex();
            const InputVariables::ControlVolume &vol = input.Volumes[arrIdx];
            if (idx.XIndex > 0 && idx.XIndex < NumXValues - 1 &&
                idx.YIndex > 0 && idx.YIndex < NumYValues - 1)
            {
                ConstructCenter(vol, xValue, yValue, leftCoeff, rightCoeff, upCoeff, downCoeff, deltaXdeltaY);
            }
            else
            {
                ConstructBoundary(vol, xValue, yValue, constantCoeff, deltaXdeltaY);
            }

            xValues_[arrIdx] = xValue;
            yValues_[arrIdx] = yValue;
            leftCoeffs_[arrIdx] = leftCoeff;
            rightCoeffs_[arrIdx] = rightCoeff;
            upCoeffs_[arrIdx] = upCoeff;
            downCoeffs_[arrIdx] = downCoeff;
            constantCoeffs_[arrIdx] = constantCoeff;
            deltaXdeltaY_[arrIdx] = deltaXdeltaY;
            sources_[arrIdx] = vol.SourceValue;
            initialValues_[arrIdx] = vol.InitialValue;
        }
    }
}

Point Grid::IndexToPoint(const Index &i) const
{
    return {xValues_[i.XIndex], yValues_[i.YIndex]};
}

Grid::Coefficients Grid::GetDiscretizationCoeffs(double alpha, double dt, const Index &idx, const Solution &soln) const
{
    const auto n = idx.ToArrayIndex();
    const double a = alpha*deltaXdeltaY_[n]/dt;
    Coefficients c;
    if (idx.XIndex > 0 && idx.XIndex < NumXValues - 1 &&
        idx.YIndex > 0 && idx.YIndex < NumYValues - 1)
    {
        c.Left = leftCoeffs_[n];
        c.Right = rightCoeffs_[n];
        c.Up = upCoeffs_[n];
        c.Down = downCoeffs_[n];
        c.Constant = a*soln(idx);
        c.Center = c.Left + c.Right + c.Up + c.Down + a - sources_[n]*deltaXdeltaY_[n]; 
    }
    else
    {
        // boundary value
        c.Left = 0;
        c.Right = 0;
        c.Up = 0;
        c.Down = 0;
        c.Center = 1;
        c.Constant = constantCoeffs_[n];
    }
    return c;
}

Grid::Index Grid::GetOrigin() const
{
    return {0, 0, *this};
}

const std::vector<double> &Grid::GetInitialValues() const
{
    return initialValues_;
}

