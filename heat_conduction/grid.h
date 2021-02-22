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

    Grid(std::vector<double> controlVolXValues, 
         std::vector<double> controlVolYValues,
         std::vector<std::vector<double>> kVertical,
         std::vector<std::vector<double>> kHorizontal,
         std::vector<std::vector<double>> sources);
    /*
     If this is what our control volumes look like (along with the thermal conductivity k and source S), 

     x----x-----x----x   horizontal faces have k = (k1,k2,k3)
     |    |     |    |
     | S1 | S2  | S3 |   vertical faces have k = (k4,k5,k6,k7)  
     |    |     |    |
     x----x-----x----x   horizontal faces have k = (k8,k9,k10)
     |    |     |    |
     |    |     |    |
     | S4 | S5  | S6 |   vertical faces have k = (k11,k12,k13,k14)
     |    |     |    |
     |    |     |    |
     x----x-----x----x   horizontal faces have k = (k15,k16,k17)
   (0,0)

     then controlVolXValues = {0, 5, 11, 16}
          controlVolYValues = {0, 4, 10}
          kVertical = {{k4,k5,k6,k7},{k11,k12,k13,k14}}
          kHorizontal = {{k1,k2,k3},{k8,k9,k10},{k15,k16,k17}}
          sources = {{S1,S2,S3},{S4,S5,S6}}
    */

    Point IndexToPoint(const Index &i) const;

    Coefficients GetDiscretizationCoeffs(double alpha, double dt, const Index &idx, const Solution &soln) const;

    Index GetOrigin() const;

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
    // we can cache at construction time some of the parameters used in the discretization coefficients
    std::vector<double> xValues_;
    std::vector<double> yValues_;
    std::vector<double> leftCoeffs_;
    std::vector<double> rightCoeffs_;
    std::vector<double> upCoeffs_;
    std::vector<double> downCoeffs_;
    std::vector<double> deltaXdeltaY_;
    std::vector<double> sources_;
};

