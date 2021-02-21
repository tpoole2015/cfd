// Chaidxer 4 of Numerical Heat Transfer and Fluid Flow by Patankar
#include <vector>
#include <algorithm>
#include <cassert>
#include "matrix.h"

struct Solution;
struct Grid
{
    struct Index
    {
        int XIndex;
        int YIndex;
    };

    struct Point
    {
        double X;
        double Y;
    };

    struct Coefficients
    {
        double Center;
        double Left;
        double Right;
        double Up;
        double Down;
        double Constant;
    };

    Grid(const std::vector<double> &xValues, const std::vector<double> &yValues);
    /*
     x----x-x--x
     |    | |  |
     |    | |  |
     x----x-x--x
     |    | |  |
     |    | |  |
     |    | |  |
     x----x-x--x
   (0,0)

     Here xValues = {0, 5, 7, 10}
          yValues = {0, 3, 7}
    */

    Point IndexToPoint(const Index &i) const;

    Coefficients GetDiscretizationCoeffs(const Index &idx, const Solution &soln) const;

    Index GetBottomLeft() const;

    bool HasLeft(const Index &idx) const;
    bool HasRight(const Index &idx) const;
    bool HasUp(const Index &idx) const;
    bool HasDown(const Index &idx) const;

    Index MoveLeft(const Index &idx) const;
    Index MoveRight(const Index &idx) const;
    Index MoveUp(const Index &idx) const;
    Index MoveDown(const Index &idx) const;

    const int NumXValues;
    const int NumYValues;

private:
    std::vector<double> xValues_;
    std::vector<double> yValues_;
};

Grid::Grid(const std::vector<double> &xValues, const std::vector<double> &yValues) :
    NumXValues(xValues.size()), 
    NumYValues(yValues.size()), 
    xValues_(xValues),
    yValues_(yValues)
{
    assert(xValues.size() && yValues.size());
    std::sort(xValues_.begin(), xValues_.end());  // sort left to right
    std::sort(yValues_.begin(), yValues_.end());  // sort bottom to top 
}

Grid::Point Grid::IndexToPoint(const Index &i) const
{
    return {xValues_[i.XIndex], yValues_[i.YIndex]};
}

Grid::Coefficients Grid::GetDiscretizationCoeffs(const Index &idx, const Solution &soln) const
{
    (void) idx;
    (void) soln;
    Coefficients c;
    return c;
}

Grid::Index Grid::GetBottomLeft() const
{
    return {0, NumYValues - 1};
}

bool Grid::HasLeft(const Index &idx) const
{
    return idx.XIndex > 0;
}

bool Grid::HasRight(const Index &idx) const
{
    return idx.XIndex < NumXValues - 1;
}

bool Grid::HasUp(const Index &idx) const
{
    return idx.YIndex > 0;
}

bool Grid::HasDown(const Index &idx) const
{
    return idx.YIndex < NumYValues - 1;
}

Grid::Index Grid::MoveLeft(const Index &idx) const
{
    return {idx.XIndex-1, idx.YIndex};
}

Grid::Index Grid::MoveRight(const Index &idx) const
{
    return {idx.XIndex+1, idx.YIndex};
}

Grid::Index Grid::MoveUp(const Index &idx) const
{
    return {idx.XIndex, idx.YIndex-1};
}

Grid::Index Grid::MoveDown(const Index &idx) const
{
    return {idx.XIndex, idx.YIndex+1};
}

struct Solution
{
    Solution(const Grid &grid);

    const Grid &GetGridReference() const;

    double operator()(const Grid::Index &idx) const;
    double& operator()(const Grid::Index &idx);
private:
    const Grid &grid_;
    std::vector<double> values_;
};

Solution::Solution(const Grid &grid)
    : grid_(grid), values_(grid.NumXValues*grid.NumYValues)
{}

const Grid &Solution::GetGridReference() const
{
    return grid_;
}

double Solution::operator()(const Grid::Index &idx) const
{
    return values_[idx.YIndex*grid_.NumXValues + idx.XIndex];
}

double& Solution::operator()(const Grid::Index &idx)
{
    return values_[idx.YIndex*grid_.NumXValues + idx.XIndex];
}

std::pair<Matrix, std::vector<double>> BuildTriDiagonalEquations(Grid::Index idx, const Solution &soln) 
{
    const Grid &grid = soln.GetGridReference();
    TridiagonalMatrix m(grid.NumYValues);
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

