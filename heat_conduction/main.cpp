// Chaidxer 4 of Numerical Heat Transfer and Fluid Flow by Patankar
#include <vector>
#include <algorithm>
#include <cassert>
#include "matrix.h"
#include "grid.h"
#include "solution.h"

std::pair<Matrix, std::vector<double>> BuildTriDiagonalEquations(Grid::Index idx, const Solution &soln) 
{
    const Grid &grid = soln.GetGridReference();
    TridiagonalMatrix m(grid.NumYValues());
    std::vector<double> rhs(m.Order);

    int i = 0;
    for ( ; grid.HasUp(idx); idx = grid.MoveUp(idx))
    {
        const Grid::Coefficients eqnCoefs = grid.GetDiscretizationCoeffs(idx, soln);

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
        if (grid.HasLeft(idx))
        {
            c += eqnCoefs.Left*soln(grid.MoveLeft(idx));
        }
        if (grid.HasRight(idx))
        {
            c += eqnCoefs.Right*soln(grid.MoveRight(idx));
        }
        rhs[i] = c;

        ++i;
    }

    return {m, rhs}; // Is m copied here? (need to check)
}

int main(int argc, char *argv[])
{
    Grid grid({0, 2, 4, 5}, {0, 1, 8, 9});
    const int numTimeSteps = 10;

    Solution soln(grid);
    for (int t = 0; t < numTimeSteps; ++t)
    {
        for (auto leftRightIdx = grid.GetBottomLeft(); grid.HasRight(leftRightIdx); leftRightIdx = grid.MoveRight(leftRightIdx))
        {
            auto triDiag = BuildTriDiagonalEquations(leftRightIdx, soln);
            triDiag.first.SolveLinear(&triDiag.second);
            
            // update soln
            int i = 0;
            for (auto downUpIdx = leftRightIdx; grid.HasUp(downUpIdx); downUpIdx = grid.MoveUp(downUpIdx))
            {
                soln(downUpIdx) = triDiag.second[i++];
            }
        }
    }    
}

