#pragma once
#include <utility>
#include <vector>
#include "grid.h"
#include "solution.h"
#include "matrix.h"

#define TOLERANCE      10e-4
#define MAX_ITERATIONS 10

namespace
{
template <class T>
std::pair<TridiagonalMatrix, std::vector<double>> BuildTriDiagonalEquationsLeftRight(Grid::Point pt, const Solution &soln, const T &coeff) 
{
    const Grid &grid = pt.GetGrid();
    std::pair<TridiagonalMatrix, std::vector<double>> returnPair = std::make_pair(TridiagonalMatrix(grid.YDim), std::vector<double>(grid.YDim, 0));

    auto &m = returnPair.first;
    auto &rhs = returnPair.second;
    int i = 0;
    for ( ; pt.InGrid(); ++pt.Y)
    {
        const Coefficients eqnCoefs = coeff(pt);
        m.GetDiagonal()[i] = eqnCoefs.Center;
        if (i > 0)
        {
            m.GetSubDiagonal()[i-1] = -1*eqnCoefs.Down; // -1 is IMPORTANT, we've moved the coefficient to the LHS!
        }
        if (i < m.Order)
        {
            m.GetSuperDiagonal()[i] = -1*eqnCoefs.Up;
        }
        rhs[i] = eqnCoefs.Constant;
        auto copyPt = pt;
        --copyPt.X;
        if (copyPt.InGrid())
        {
            rhs[i] += eqnCoefs.Left*soln(copyPt);
        }
        copyPt = pt;
        ++copyPt.X;
        if (copyPt.InGrid())
        {
            rhs[i] += eqnCoefs.Right*soln(copyPt);
        }
        ++i;
    }
    return returnPair;
}

template <class T>
std::pair<TridiagonalMatrix, std::vector<double>> BuildTriDiagonalEquationsDownUp(Grid::Point pt, const Solution &soln, const T &coeff) 
{
    const Grid &grid = pt.GetGrid();
    std::pair<TridiagonalMatrix, std::vector<double>> returnPair = std::make_pair(TridiagonalMatrix(grid.XDim), std::vector<double>(grid.XDim, 0));

    auto &m = returnPair.first;
    auto &rhs = returnPair.second;
    int i = 0;
    for ( ; pt.InGrid(); ++pt.X)
    {
        const Coefficients eqnCoefs = coeff(pt);
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
        auto copyPt = pt;
        --copyPt.Y;
        if (copyPt.InGrid())
        {
            c += eqnCoefs.Down*soln(copyPt);
        }
        copyPt = pt;
        ++copyPt.Y;
        if (copyPt.InGrid())
        {
            c += eqnCoefs.Up*soln(copyPt);
        }
        rhs[i] = c;

        ++i;
    }
    return returnPair;
}
}

namespace Solver
{
    template <class T>
    int Solve(const T &coeffs, const Grid &grid, Solution &soln)
    {
        // Gauss-Seidel (alternate solving left/right and down/up)
        int numIterations = 0;
        bool workLeftRight = true;
        while (numIterations++ < MAX_ITERATIONS)
        {
            double sumResiduals = 0;
            if (workLeftRight)
            {
                for (auto leftRightPt = grid.GetBottomLeft(); leftRightPt.InGrid(); ++leftRightPt.X)
                {
                    auto triDiag = BuildTriDiagonalEquationsLeftRight(leftRightPt, soln, coeffs);
                    triDiag.first.SolveLinear(&triDiag.second);
                   
                    // update soln
                    int i = 0;
                    for (auto downUpPt = leftRightPt; downUpPt.InGrid(); ++downUpPt.Y)
                    {
                        const double newValue = triDiag.second[i++];
                        sumResiduals += (newValue - soln(downUpPt))*(newValue - soln(downUpPt));
                        soln(downUpPt) = newValue;
                    }
                }
            } 
            else 
            {
                for (auto downUpPt = grid.GetBottomLeft(); downUpPt.InGrid(); ++downUpPt.Y)
                {
                    auto triDiag = BuildTriDiagonalEquationsDownUp(downUpPt, soln, coeffs);
                    triDiag.first.SolveLinear(&triDiag.second);
                   
                    // update soln
                    int i = 0;
                    for (auto leftRightPt = downUpPt; leftRightPt.InGrid(); ++leftRightPt.X)
                    {
                        const double newValue = triDiag.second[i++];
                        sumResiduals += (newValue - soln(leftRightPt))*(newValue - soln(leftRightPt));
                        soln(leftRightPt) = newValue;
                    }
                }
            }
            workLeftRight = !workLeftRight;
            const double err = sumResiduals / static_cast<double>(grid.XDim*grid.YDim);
            if (err < TOLERANCE)
            {
                break;
            }
        }
        return numIterations;
    }
}
