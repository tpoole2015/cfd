#pragma once
#include <vector>

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

    Grid(const std::vector<double> &controlVolXValues, const std::vector<double> &controlVolYValues);
    /*
     If this is what our control volumes look like, 

     x----x-x--x
     |    | |  |
     |    | |  |
     x----x-x--x
     |    | |  |
     |    | |  |
     |    | |  |
     x----x-x--x
   (0,0)

     then controlVolXValues = {0, 5, 7, 10}
          controlVolYValues = {0, 3, 7}
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

    int NumXValues() const;
    int NumYValues() const;

private:
    std::vector<double> controlVolXValues_;
    std::vector<double> controlVolYValues_;
    std::vector<double> xValues_;
    std::vector<double> yValues_;
};

