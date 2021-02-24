// Chaidxer 4 of Numerical Heat Transfer and Fluid Flow by Patankar
#include <vector>
#include <algorithm>
#include <cassert>
#include "matrix.h"
#include "grid.h"
#include "solution.h"

#define TOLERANCE      10e-6
#define MAX_ITERATIONS 10

std::pair<Matrix, std::vector<double>> BuildTriDiagonalEquations(Grid::Index idx, const Solution &prev, const Solution &next, double alpha, double dt) 
{
    const Grid &grid = idx.GridReference;
    TridiagonalMatrix m(grid.NumYValues);
    std::vector<double> rhs(m.Order);

    int i = 0;
    for ( ; idx.HasUp(); idx.MoveUp())
    {
        const Grid::Coefficients eqnCoefs = grid.GetDiscretizationCoeffs(alpha, dt, idx, prev);

        m.GetDiagonal()[i] = eqnCoefs.Center;
        if (i > 0)
        {
            m.GetSubDiagonal()[i-1] = eqnCoefs.Down;
        }
        if (i < m.Order)
        {
            m.GetSuperDiagonal()[i] = eqnCoefs.Up;
        }
        double c = eqnCoefs.Constant;
        if (idx.HasLeft())
        {
            auto copyIdx = idx;
            copyIdx.MoveLeft();
            c += eqnCoefs.Left*next(copyIdx);
        }
        if (idx.HasRight())
        {
            auto copyIdx = idx;
            copyIdx.MoveRight();
            c += eqnCoefs.Right*next(copyIdx);
        }
        rhs[i] = c;

        ++i;
    }

    return {m, rhs}; // Is m copied here? (need to check)
}

int main(int argc, char *argv[])
{
    const double k = 2;
    const double S = 3; 
    const double alpha = 1;
    const double dt = 0.1;
    const double t = 5;

    InputVariables input;
    Grid grid(input);

    Solution next(grid);
    Solution prev(grid);
    for (int i = 0; i < static_cast<int>(t / dt); ++i)
    {
        int numIterations = 0;
        while (numIterations < MAX_ITERATIONS)
        {
            double sumResiduals = 0;
            for (auto leftRightIdx = grid.GetOrigin(); leftRightIdx.HasRight(); leftRightIdx.MoveRight())
            {
                auto triDiag = BuildTriDiagonalEquations(leftRightIdx, prev, next, alpha, dt);
                triDiag.first.SolveLinear(&triDiag.second);
                
                // update soln
                int i = 0;
                for (auto downUpIdx = leftRightIdx; downUpIdx.HasUp(); downUpIdx.MoveUp())
                {
                    const double newValue = triDiag.second[i++];
                    sumResiduals += (newValue - next(downUpIdx))*(newValue - next(downUpIdx));
                    next(downUpIdx) = newValue;
                }
            }

            const double err = sumResiduals / (grid.NumXValues*grid.NumYValues);
            if (err < TOLERANCE)
            {
                break;
            }
            ++numIterations;
        }
        prev = next;
    }    
}

