#pragma once
#include <array>
#include <vector>
#include "grid.h"

struct InputVariables
{
    struct ControlVolume
    {
        Point2d TopLeft;
        Point2d BottomRight;
        double InitialValue;
        double BoundaryValue;
        double SourceValue;
        std::array<double, 4> K; // 0 = Left face, 1 = Right face, 2 = Up face, 3 = Down face
    };
    int NumXVols; 
    int NumYVols;
    std::vector<ControlVolume> Volumes; // [Y*NumXVols + X]
};


