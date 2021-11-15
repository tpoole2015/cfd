#include <cassert>
#include <iostream>
#include "heat_equation.h"
#include "solver.h"

namespace {
void ConstructCenter(const InputVariables::ControlVolume &vol,
                     double &midPointXValue,
                     double &midPointYValue,
                     double &leftCoeff,
                     double &rightCoeff,
                     double &upCoeff,
                     double &downCoeff,
                     double &deltaXdeltaY)
{
    midPointXValue = 0.5*(vol.BottomRight.X + vol.TopLeft.X);
    midPointYValue = 0.5*(vol.BottomRight.Y + vol.TopLeft.Y);

    const double deltaY_d = midPointYValue - vol.BottomRight.Y;
    const double deltaY_u = vol.TopLeft.Y - midPointYValue;
    const double deltaY = deltaY_d + deltaY_u;
    const double deltaX_l = midPointXValue - vol.TopLeft.X;
    const double deltaX_r = vol.BottomRight.X - midPointXValue;
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
                       double &midPointXValue,
                       double &midPointYValue,
                       double &constantCoeff,
                       double &deltaXdeltaY)
{
    midPointXValue = 0.5*(vol.BottomRight.X + vol.TopLeft.X);
    midPointYValue = 0.5*(vol.BottomRight.Y + vol.TopLeft.Y);

    const double deltaX = vol.BottomRight.X - vol.TopLeft.X;
    const double deltaY = vol.TopLeft.Y - vol.BottomRight.Y;

    constantCoeff = vol.BoundaryValue;
    deltaXdeltaY = deltaX * deltaY;
}
}

HeatEquation::HeatEquation(const InputVariables &input, double alpha, double dt) :
    grid_(input.NumXVols, input.NumYVols),
    alpha_(alpha),
    dt_(dt),
    prev_(&grid_),
    next_(&grid_),
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
   
    double midPointXValue = 0;
    double midPointYValue = 0;
    double leftCoeff = 0;
    double rightCoeff = 0;
    double upCoeff = 0;
    double downCoeff = 0;
    double constantCoeff = 0;
    double deltaXdeltaY = 0;

    for (auto topBottom = grid_.GetTopLeft(); topBottom.InGrid(); --topBottom.Y)
    {
        for (auto pt = topBottom; pt.InGrid(); ++pt.X)
        {
            const auto arrIdx = pt.ToArrayIndex();
            const InputVariables::ControlVolume &vol = input.Volumes[arrIdx];
            if (pt.IsBoundary())
            {
                ConstructBoundary(vol, midPointXValue, midPointYValue, constantCoeff, deltaXdeltaY);
            }
            else
            {
                ConstructCenter(vol, midPointXValue, midPointYValue, leftCoeff, rightCoeff, upCoeff, downCoeff, deltaXdeltaY);
            }

            grid_.xValues_[arrIdx] = midPointXValue;
            grid_.yValues_[arrIdx] = midPointYValue;
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

    prev_.Init(initialValues_);
    next_.Init(initialValues_);
}

Coefficients HeatEquation::operator()(const Grid::Point& pt) const
{
    const auto n = pt.ToArrayIndex();
    const double a = alpha_*deltaXdeltaY_[n]/dt_;
    Coefficients c;
    if (pt.IsBoundary())
    {
        // boundary value
        c.Left = 0;
        c.Right = 0;
        c.Up = 0;
        c.Down = 0;
        c.Center = 1;
        c.Constant = constantCoeffs_[n]; // ignoring boundary sources
    }
    else
    {
        c.Left = leftCoeffs_[n];
        c.Right = rightCoeffs_[n];
        c.Up = upCoeffs_[n];
        c.Down = downCoeffs_[n];
        c.Constant = a*prev_(pt) + sources_[n]*deltaXdeltaY_[n];
        c.Center = c.Left + c.Right + c.Up + c.Down + a; 
    }

    return c;
}

int HeatEquation::SolveOneStep()
{
    const int numIterations = Solver::Solve(*this, this->grid_, this->next_);
    this->prev_ = this->next_;
    return numIterations;
}


