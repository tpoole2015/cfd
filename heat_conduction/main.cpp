// Chaidxer 4 of Numerical Heat Transfer and Fluid Flow by Patankar
#include <vector>
#include <algorithm>
#include <cassert>
#include <iostream>
#include <sstream>
#include <cmath>
#include "matrix.h"
#include "grid.h"
#include "solution.h"

#define TOLERANCE      10e-5
#define MAX_ITERATIONS 100

std::pair<TridiagonalMatrix, std::vector<double>> BuildTriDiagonalEquations(Grid::Index idx, const Solution &prev, const Solution &next, double alpha, double dt) 
{
    const Grid &grid = *(idx.GridPtr);
    TridiagonalMatrix m(grid.NumYValues);
    std::vector<double> rhs(m.Order);

    int i = 0;
    for ( ; idx.InGrid(); idx.MoveUp())
    {
        const Grid::Coefficients eqnCoefs = grid.GetDiscretizationCoeffs(alpha, dt, idx, prev);
        const double sum = std::fabs(eqnCoefs.Up) + std::fabs(eqnCoefs.Down) + std::fabs(eqnCoefs.Left) + std::fabs(eqnCoefs.Right);
        assert(sum < std::fabs(eqnCoefs.Center)); // numerical requirement

        m.GetDiagonal()[i] = eqnCoefs.Center;
        if (i > 0)
        {
            m.GetSubDiagonal()[i-1] = -1*eqnCoefs.Down; // -1 is IMPORTANT, we've moved the coefficient to the LHS!
        }
        if (i < m.Order)
        {
            m.GetSuperDiagonal()[i] = -1*eqnCoefs.Up;
        }
        double c = eqnCoefs.Constant;
        auto copyIdx = idx;
        copyIdx.MoveLeft();
        if (copyIdx.InGrid())
        {
            c += eqnCoefs.Left*next(copyIdx);
        }
        copyIdx = idx;
        copyIdx.MoveRight();
        if (copyIdx.InGrid())
        {
            c += eqnCoefs.Right*next(copyIdx);
        }
        rhs[i] = c;

        ++i;
    }

    return {m, rhs}; // Is m copied here? (need to check)
}

InputVariables BuildInputVariables()
{
    const double k = 1;
    const double I = 150;
    const double bvTop = 0, bvBottom = 300, bvLeft = 0, bvRight = 300; 
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
        tvol.SourceValue = 0;
        tvol.K = {{k, k, k, k}};
        tvol.BVType = InputVariables::ControlVolume::Scalar; 
        tvol.BoundaryValues = {bvTop};
 
        // Bottom
        const int bottomY = 0;
        auto &bvol = iv.Volumes[bottomY*iv.NumXVols + X];
        bvol.TopLeft = {X*volumeWidth, (bottomY + 1)*volumeHeight};
        bvol.BottomRight = {(X+1)*volumeWidth, bottomY*volumeHeight};
        bvol.InitialValue = 0;
        bvol.SourceValue = 0;
        bvol.K = {{k, k, k, k}};
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
        lvol.SourceValue = 0;
        lvol.K = {{k, k, k, k}};
        lvol.BVType = InputVariables::ControlVolume::Scalar; 
        lvol.BoundaryValues = {bvLeft};

        // Right 
        const int rightX = iv.NumXVols-1;
        auto &rvol = iv.Volumes[Y*iv.NumXVols + rightX];
        rvol.TopLeft = {(rightX-1)*volumeWidth, (Y+1)*volumeHeight};
        rvol.BottomRight = {rightX*volumeWidth, Y*volumeHeight};
        rvol.InitialValue = 0;
        rvol.SourceValue = 0;
        rvol.K = {{k, k, k, k}};
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
            vol.SourceValue = 0;
            vol.K = {{k, k, k, k}};
            vol.BVType = InputVariables::ControlVolume::Scalar; 
        }
        // don't heat up the top half
        for (int Y = iv.NumYVols/2; Y < iv.NumYVols - 1; ++Y)
        {
            auto &vol = iv.Volumes[Y*iv.NumXVols + X];
            vol.TopLeft = {X*volumeWidth, (Y+1)*volumeHeight};
            vol.BottomRight = {(X+1)*volumeWidth, Y*volumeHeight};
            vol.InitialValue = 0;
            vol.SourceValue = 0;
            vol.K = {{k, k, k, k}};
            vol.BVType = InputVariables::ControlVolume::Scalar; 
        }
    }
    return iv;
}

int main(int argc, char *argv[])
{
    if (argc < 1)
    {
        std::cout << "Usage: <output file>\n";
        return 1;
    }

    const std::string fn(argv[1]);
    const double alpha = 1;
    const double dt = 0.1;
    const double t = 10;

    Grid grid(BuildInputVariables());

    Solution next(grid);
    Solution prev(grid);
 
    int totalNumIterations = 0; 
    const int numSteps = static_cast<int>(t/dt); 
    for (int i = 0; i < numSteps; ++i)
    {
        int numIterations = 0;
        while (numIterations++ < MAX_ITERATIONS)
        {
            double sumResiduals = 0;
            for (auto leftRightIdx = grid.GetBottomLeft(); leftRightIdx.InGrid(); leftRightIdx.MoveRight())
            {
                auto triDiag = BuildTriDiagonalEquations(leftRightIdx, prev, next, alpha, dt);
                triDiag.first.SolveLinear(&triDiag.second);
               
                // update soln
                int i = 0;
                for (auto downUpIdx = leftRightIdx; downUpIdx.InGrid(); downUpIdx.MoveUp())
                {
                    const double newValue = triDiag.second[i++];
                    sumResiduals += (newValue - next(downUpIdx))*(newValue - next(downUpIdx));
                    next(downUpIdx) = newValue;
                }
            }

            const double err = sumResiduals / static_cast<double>(grid.NumXValues*grid.NumYValues);
            if (err < TOLERANCE)
            {
                break;
            }
        }
        totalNumIterations += numIterations;
        prev = next;

        std::stringstream fileName;
        fileName << fn << "_" << i;
        next.WriteHeatMapToFile(fileName.str());
    }    
    std::cout << "Avg num iterations for convergence: " << totalNumIterations/numSteps << "\n";
}

