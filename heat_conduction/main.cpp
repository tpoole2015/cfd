// Chaidxer 4 of Numerical Heat Transfer and Fluid Flow by Patankar
#include <vector>
#include <algorithm>
#include "matrix.h"

struct Grid
{
    struct Index
    {
        // TODO: understand brace initilization
        int XIndex{0};
        int YIndex{0};
    };

    struct Point
    {
        double X{0.0};
        double Y{0.0};
    }

    struct Coefficients
    {
        double Center{0.0};
        double Left{0.0};
        double Right{0.0};
        double Up{0.0};
        double Down{0.0};
        double Constant{0.0};
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
    vector<double> xValues_;
    vector<double> yValues_;
};

Grid::Grid(const std::vector<double> &xValues, const std::vector<double> &yValues) :
    NumXValues(xValues.size()), 
    NumYValues(yValues.size()), 
    xValues_(xValues),
    yValues_(yValues)
{
    assert(xValues.size() && yValues.size());
    std::sort(xValues.begin(), xValues.end());  // sort left to right
    std::sort(yValues.begin(), yValues.end());  // sort bottom to top 
}

Grid::Point Grid::IndexToPoint(const Index &i) const
{
    return {xValues_[i.XIndex], yValues_[i.YIndex];
}

Coefficients GetDiscretizationCoeffs(const Index &idx, const Solution &soln)
{
    (void) idx;
    (void) soln;
    Coefficients c;
    return c;
}

Index GetBottomLeft()
{
    return {0, NumYValues - 1};
}

bool HasLeft(const Index &idx)
{
    return idx.XIndex > 0;
}

bool HasRight(const Index &idx)
{
    return idx.XIndex < NumXValues - 1;
}

bool HasUp(const Index &idx)
{
    return idx.YIndex > 0;
}

bool HasDown(const Index &idx);
{
    return idx.YIndex < NumYValues - 1;
}

Index MoveLeft(const Index &idx)
{
    return {idx.XIndex-1, idx.YIndex};
}

Index MoveRight(const Index &idx)
    return {idx.XIndex+1, idx.YIndex};
}

Index MoveUp(const Index &idx)
{
    return {idx.XIndex, idx.YIndex-1};
}

Index MoveDown(const Index &idx)
    return {idx.XIndex, idx.YIndex+1};
}

struct Solution
{
    Solution(const Grid &grid);

    const Grid &GetGridReference();

    double operator()(const Grid::Index &idx) const;
    double& operator()(const Grid::Index &idx);
private:
    const Grid &grid_;
    std::vector<double> values_;
};

Solution::Solution(const Grid &grid)
    : grid_(grid), values_[grid.NumXValues*grid.NumYValues]
{}

const Grid &GetGridReference()
{
    return grid_;
}

double operator()(const Grid::Index &idx) const
{
    return (*this)(idx);
}

double& operator()(const Grid::Index &idx)
{
    return values_[idx.YIndex*grid_.NumXValues + idx.XIndex];
}

std::pair<Matrix, std::vector<double>> BuildTriDiagonalEquations(Grid::Index idx, const Solution &soln) 
{
    TridiagonalMatrix m(grid.NumYIndexs);
    std::vector<double> rhs(m.Order);
    const Grid &grid = soln.GetGridReference();

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
        if (g.HasRight(idx))
        {
            c += eqnCoefs.Right*soln(grid.MoveRight(idx));
        }
        rhs[i] = c;

        ++i;
    }

    return {m, c}; // Is m copied here? (need to check)
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
            const auto triDiag = BuildTriDiagonalMatrix(leftRightIdx, soln);
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

