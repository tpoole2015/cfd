// Chaidxer 4 of Numerical Heat Transfer and Fluid Flow by Patankar
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <cmath>
#include <chrono>
#include "matrix.h"
#include "grid.h"
#include "solution.h"

#define TOLERANCE      10e-6
#define MAX_ITERATIONS 100

std::pair<TridiagonalMatrix, std::vector<double>> BuildTriDiagonalEquationsLeftRight(Grid::Index idx, const Solution &prev, const Solution &next, double alpha, double dt) 
{
    const Grid &grid = *(idx.GridPtr);
    std::pair<TridiagonalMatrix, std::vector<double>> returnPair = 
        std::make_pair(TridiagonalMatrix(grid.NumYValues), std::vector<double>(grid.NumYValues, 0));

    auto &m = returnPair.first;
    auto &rhs = returnPair.second;
    int i = 0;
    for ( ; idx.InGrid(); idx.MoveUp())
    {
        const Grid::Coefficients eqnCoefs = grid.GetDiscretizationCoeffs(alpha, dt, idx, prev);

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

    return returnPair;
}

std::pair<TridiagonalMatrix, std::vector<double>> BuildTriDiagonalEquationsDownUp(Grid::Index idx, const Solution &prev, const Solution &next, double alpha, double dt) 
{
    const Grid &grid = *(idx.GridPtr);
    std::pair<TridiagonalMatrix, std::vector<double>> returnPair = 
        std::make_pair(TridiagonalMatrix(grid.NumXValues), std::vector<double>(grid.NumXValues, 0));

    auto &m = returnPair.first;
    auto &rhs = returnPair.second;
    int i = 0;
    for ( ; idx.InGrid(); idx.MoveRight())
    {
        const Grid::Coefficients eqnCoefs = grid.GetDiscretizationCoeffs(alpha, dt, idx, prev);

        m.GetDiagonal()[i] = eqnCoefs.Center;
        if (i > 0)
        {
            m.GetSubDiagonal()[i-1] = -1*eqnCoefs.Left; // -1 is IMPORTANT, we've moved the coefficient to the LHS!
        }
        if (i < m.Order)
        {
            m.GetSuperDiagonal()[i] = -1*eqnCoefs.Right;
        }
        double c = eqnCoefs.Constant;
        auto copyIdx = idx;
        copyIdx.MoveDown();
        if (copyIdx.InGrid())
        {
            c += eqnCoefs.Down*next(copyIdx);
        }
        copyIdx = idx;
        copyIdx.MoveUp();
        if (copyIdx.InGrid())
        {
            c += eqnCoefs.Up*next(copyIdx);
        }
        rhs[i] = c;

        ++i;
    }

    return returnPair;
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

    const std::string fileName(argv[1]);
    const double alpha = 1;
    const double dt = 0.1;
    const double t = 10;

    Grid grid(BuildInputVariables());

    Solution next(grid);
    Solution prev(grid);
 
    int totalNumIterations = 0; 
    const int numSteps = static_cast<int>(t/dt); 
    const auto start = std::chrono::steady_clock::now();
    for (int i = 0; i < numSteps; ++i)
    {
        // Gauss-Seidel (alternate solving left/right and down/up)
        int numIterations = 0;
        bool workLeftRight = true;
        while (numIterations++ < MAX_ITERATIONS)
        {
            double sumResiduals = 0;
            if (workLeftRight)
            {
                for (auto leftRightIdx = grid.GetBottomLeft(); leftRightIdx.InGrid(); leftRightIdx.MoveRight())
                {
                    auto triDiag = BuildTriDiagonalEquationsLeftRight(leftRightIdx, prev, next, alpha, dt);
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
            } 
            else 
            {
                for (auto downUpIdx = grid.GetBottomLeft(); downUpIdx.InGrid(); downUpIdx.MoveUp())
                {
                    auto triDiag = BuildTriDiagonalEquationsDownUp(downUpIdx, prev, next, alpha, dt);
                    triDiag.first.SolveLinear(&triDiag.second);
                   
                    // update soln
                    int i = 0;
                    for (auto leftRightIdx = downUpIdx; leftRightIdx.InGrid(); leftRightIdx.MoveRight())
                    {
                        const double newValue = triDiag.second[i++];
                        sumResiduals += (newValue - next(leftRightIdx))*(newValue - next(leftRightIdx));
                        next(leftRightIdx) = newValue;
                    }
                }
            }
            workLeftRight = !workLeftRight;
            const double err = sumResiduals / static_cast<double>(grid.NumXValues*grid.NumYValues);
            if (err < TOLERANCE)
            {
                break;
            }
        }
        totalNumIterations += numIterations;
        prev = next;
    }    
    const auto end = std::chrono::steady_clock::now();
    std::cout << "took " << 1000*std::chrono::duration<double>(end - start).count() << "ms\n";

    next.WriteHeatMapToFile(fileName);
    std::cout << "Avg num iterations for convergence: " << totalNumIterations/numSteps << "\n";
}

