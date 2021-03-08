// Chapter 4 of Numerical Heat Transfer and Fluid Flow by Patankar
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <cmath>
#include <chrono>
#include <cstdlib>
#include <cassert>
#include "matrix.h"
#include "grid.h"
#include "solution.h"
#include "solver.h"

struct InputVariables
{
    struct ControlVolume
    {
        enum BoundaryValueType { Scalar, Flux };

        Point TopLeft;
        Point BottomRight;

        double InitialValue;
        double SourceValue;
        std::array<double, 4> K; // 0 = Left face, 1 = Right face, 2 = Up face, 3 = Down face

        BoundaryValueType BVType; // Only used for boundary volumes

        // If type is Scalar: BoundaryValues = { value at grid point }
        // If type is Flux:
        //  If position is TopLeft: BoundaryValues = { flux top, flux left }
        //  If position is TopRight: BoundaryValues = { flux top, flux right }
        //  If position is BottomLeft: BoundaryValues = { flux bottom, flux left}
        //  If position is BottomRight: BoundaryValues = { flux bottom, flux right}
        //  If position is Left: BoundaryValues = { flux left }
        //  If position is Right: BoundaryValues = { flux right }
        //  If position is Bottom: BoundaryValues = { flux bottom }
        //  If position is Top: BoundaryValues = { flux top }
        // Here position is determined by location of our ControlVolume in the sorted array Volumes
        std::vector<double> BoundaryValues;
    };
    int NumXVols; 
    int NumYVols;
    std::vector<ControlVolume> Volumes; // [Y*NumXVols + X]
};

struct HeatEquation 
{
    HeatEquation(InputVariables input, double alpha, double dt);
    Coefficients operator()(const Grid::Point& pt) const;

    Grid grid_;
    double alpha_;
    double dt_;
    Solution prev_;
    Solution next_;

    std::vector<double> leftCoeffs_;
    std::vector<double> rightCoeffs_;
    std::vector<double> upCoeffs_;
    std::vector<double> downCoeffs_;
    std::vector<double> constantCoeffs_;
    std::vector<double> deltaXdeltaY_;
    std::vector<double> sources_;
    std::vector<double> initialValues_;
};

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

HeatEquation::HeatEquation(InputVariables input, double alpha, double dt) :
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
   
    double xValue = 0;
    double yValue = 0;
    double leftCoeff = 0;
    double rightCoeff = 0;
    double upCoeff = 0;
    double downCoeff = 0;
    double constantCoeff = 0;
    double deltaXdeltaY = 0;

    Grid::Point pt{0, 0, &grid_};
    for (pt.J = 0 ; pt.J < input.NumYVols; ++pt.J)
    {
        for (pt.I = 0 ; pt.I < input.NumXVols; ++pt.I)
        {
            const auto arrIdx = pt.ToArrayIndex();
            const InputVariables::ControlVolume &vol = input.Volumes[arrIdx];
            if (pt.I > 0 && pt.I < grid_.XDim - 1 && pt.J > 0 && pt.J < grid_.YDim - 1)
            {
                ConstructCenter(vol, xValue, yValue, leftCoeff, rightCoeff, upCoeff, downCoeff, deltaXdeltaY);
            }
            else
            {
                ConstructBoundary(vol, xValue, yValue, constantCoeff, deltaXdeltaY);
            }

            grid_.xValues_[arrIdx] = xValue;
            grid_.yValues_[arrIdx] = yValue;
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
    if (pt.I > 0 && pt.I < grid_.XDim - 1 && pt.J > 0 && pt.J < grid_.YDim - 1)
    {
        c.Left = leftCoeffs_[n];
        c.Right = rightCoeffs_[n];
        c.Up = upCoeffs_[n];
        c.Down = downCoeffs_[n];
        c.Constant = a*prev_(pt) + sources_[n]*deltaXdeltaY_[n];
        c.Center = c.Left + c.Right + c.Up + c.Down + a; 
    }
    else
    {
        // boundary value
        c.Left = 0;
        c.Right = 0;
        c.Up = 0;
        c.Down = 0;
        c.Center = 1;
        c.Constant = constantCoeffs_[n]; // ignoring boundary sources
    }
    return c;
}

InputVariables BuildInputVariables()
{
    const double kL = 0, kR = 0, kU = 10, kD = 10;
    const double I = 0;
    const double S = 100;
    const double bvTop = 0, bvBottom = 0, bvLeft = 0, bvRight = 0; 
    InputVariables iv;
    iv.NumXVols = 1000;
    iv.NumYVols = 100;
    iv.Volumes.resize(iv.NumXVols*iv.NumYVols);

    const double volumeWidth = 1 / static_cast<double>(iv.NumXVols);
    const double volumeHeight = 1 / static_cast<double>(iv.NumYVols);

    // Top & Bottom
    for (int X = 0; X < iv.NumXVols; ++X)
    {
        // Top 
        const int topY = iv.NumYVols-1;
        auto &tvol = iv.Volumes[topY*iv.NumXVols + X];
        tvol.TopLeft = {X*volumeWidth, (topY+1)*volumeHeight};
        tvol.BottomRight = {(X+1)*volumeWidth, topY*volumeHeight};
        tvol.InitialValue = 0;
        tvol.K = {{kL, kR, kU, kD}};
        tvol.BVType = InputVariables::ControlVolume::Scalar; 
        tvol.BoundaryValues = {bvTop};
 
        // Bottom
        const int bottomY = 0;
        auto &bvol = iv.Volumes[bottomY*iv.NumXVols + X];
        bvol.TopLeft = {X*volumeWidth, (bottomY + 1)*volumeHeight};
        bvol.BottomRight = {(X+1)*volumeWidth, bottomY*volumeHeight};
        bvol.InitialValue = 0;
        bvol.K = {{kL, kR, kU, kD}};
        bvol.BVType = InputVariables::ControlVolume::Scalar; 
        bvol.BoundaryValues = {bvBottom};
    }

    // Left & Right
    for (int Y = 1; Y < iv.NumYVols - 1; ++Y) // skip Top & Bottom as we did those above
    {
        // Left 
        const int leftX = 0;
        auto &lvol = iv.Volumes[Y*iv.NumXVols + leftX];
        lvol.TopLeft = {leftX*volumeWidth, (Y+1)*volumeHeight};
        lvol.BottomRight = {(leftX+1)*volumeWidth, Y*volumeHeight};
        lvol.InitialValue = 0;
        lvol.K = {{kL, kR, kU, kD}};
        lvol.BVType = InputVariables::ControlVolume::Scalar; 
        lvol.BoundaryValues = {bvLeft};

        // Right 
        const int rightX = iv.NumXVols-1;
        auto &rvol = iv.Volumes[Y*iv.NumXVols + rightX];
        rvol.TopLeft = {(rightX-1)*volumeWidth, (Y+1)*volumeHeight};
        rvol.BottomRight = {rightX*volumeWidth, Y*volumeHeight};
        rvol.InitialValue = 0;
        rvol.K = {{kL, kR, kU, kD}};
        rvol.BVType = InputVariables::ControlVolume::Scalar; 
        rvol.BoundaryValues = {bvRight};
    }

    // Now fill in the interaior
    for (int X = 1; X < iv.NumXVols - 1; ++X)
    {
        // heat up the bottom half
        for (int Y = 1; Y < iv.NumYVols/2; ++Y)
        {
            auto &vol = iv.Volumes[Y*iv.NumXVols + X];
            vol.TopLeft = {X*volumeWidth, (Y+1)*volumeHeight};
            vol.BottomRight = {(X+1)*volumeWidth, Y*volumeHeight};
            vol.InitialValue = I;
            vol.SourceValue = S;
            vol.K = {{kL, kR, kU, kD}};
            vol.BVType = InputVariables::ControlVolume::Scalar; 
        }
        // don't heat up the top half
        for (int Y = iv.NumYVols/2; Y < iv.NumYVols - 1; ++Y)
        {
            auto &vol = iv.Volumes[Y*iv.NumXVols + X];
            vol.TopLeft = {X*volumeWidth, (Y+1)*volumeHeight};
            vol.BottomRight = {(X+1)*volumeWidth, Y*volumeHeight};
            vol.InitialValue = I;
            vol.SourceValue = 0;
            vol.K = {{kL, kR, kU, kD}};
            vol.BVType = InputVariables::ControlVolume::Scalar; 
        }
    }
    return iv;
}

int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        std::cout << "Usage: <output file> <time>\n";
        return 1;
    }

    const std::string fileName(argv[1]);
    const double alpha = 1;
    const double dt = 0.1;
    const double t = std::strtod(argv[2], nullptr);

    HeatEquation he(BuildInputVariables(), alpha, dt);

    int totalNumIterations = 0; 
    const int numSteps = static_cast<int>(t/dt); 
    const auto start = std::chrono::steady_clock::now();
    for (int i = 0; i < numSteps; ++i)
    {
        totalNumIterations += Solver::Solve(he, he.grid_, he.next_);
        he.prev_ = he.next_;
    }
    const auto end = std::chrono::steady_clock::now();
    std::cout << "took " << 1000*std::chrono::duration<double>(end - start).count() << "ms\n";

    he.next_.WriteHeatMapToFile(fileName);
    std::cout << "Avg num iterations for convergence: " << totalNumIterations/numSteps << "\n";
}

