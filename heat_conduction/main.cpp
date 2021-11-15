// Chapter 4 of Numerical Heat Transfer and Fluid Flow by Patankar
#include <iostream>
#include <cmath>
#include <chrono>
#include <cstdlib>
#include "input_variables.h"
#include "heat_equation.h"

InputVariables BuildInputVariables()
{
    const double kL = 1, kR = 1, kU = 10, kD = 10;
    const double I = 0;
    const double bvTop = 1000, bvBottom = 0, bvLeft = 0, bvRight = 0; 
    InputVariables iv;
    iv.NumXVols = 1000;
    iv.NumYVols = 1000;
    iv.Volumes.resize(iv.NumXVols*iv.NumYVols);
    const double dX = 1 / static_cast<double>(iv.NumXVols);
    const double dY = 1 / static_cast<double>(iv.NumYVols);
    const double S = 1 / (dX*dY);

    // Top & Bottom
    for (int X = 0; X < iv.NumXVols; ++X)
    {
        // Top 
        const int topY = iv.NumYVols-1;
        auto &tvol = iv.Volumes[topY*iv.NumXVols + X];
        tvol.TopLeft = {X*dX, (topY+1)*dY};
        tvol.BottomRight = {(X+1)*dX, topY*dY};
        tvol.InitialValue = 0;
        tvol.SourceValue = 0;
        tvol.K = {{kL, kR, kU, kD}};
        tvol.BoundaryValue = bvTop;
 
        // Bottom
        const int bottomY = 0;
        auto &bvol = iv.Volumes[bottomY*iv.NumXVols + X];
        bvol.TopLeft = {X*dX, (bottomY + 1)*dY};
        bvol.BottomRight = {(X+1)*dX, bottomY*dY};
        bvol.InitialValue = 0;
        bvol.SourceValue = 0;
        bvol.K = {{kL, kR, kU, kD}};
        bvol.BoundaryValue = bvBottom;
    }

    // Left & Right
    for (int Y = 1; Y < iv.NumYVols - 1; ++Y) // skip Top & Bottom as we did those above
    {
        // Left 
        const int leftX = 0;
        auto &lvol = iv.Volumes[Y*iv.NumXVols + leftX];
        lvol.TopLeft = {leftX*dX, (Y+1)*dY};
        lvol.BottomRight = {(leftX+1)*dX, Y*dY};
        lvol.InitialValue = 0;
        lvol.SourceValue = 0;
        lvol.K = {{kL, kR, kU, kD}};
        lvol.BoundaryValue = bvLeft;

        // Right 
        const int rightX = iv.NumXVols-1;
        auto &rvol = iv.Volumes[Y*iv.NumXVols + rightX];
        rvol.TopLeft = {(rightX-1)*dX, (Y+1)*dY};
        rvol.BottomRight = {rightX*dX, Y*dY};
        rvol.InitialValue = 0;
        rvol.SourceValue = 0;
        rvol.K = {{kL, kR, kU, kD}};
        rvol.BoundaryValue = bvRight;
    }

    // Now fill in the interior
    for (int X = 1; X < iv.NumXVols - 1; ++X)
    {
        // heat up the bottom half
        for (int Y = 1; Y < iv.NumYVols/2; ++Y)
        {
            auto &vol = iv.Volumes[Y*iv.NumXVols + X];
            vol.TopLeft = {X*dX, (Y+1)*dY};
            vol.BottomRight = {(X+1)*dX, Y*dY};
            vol.InitialValue = I;
            vol.SourceValue = S;
            vol.K = {{kL, kR, kU, kD}};
        }
        // don't heat up the top half
        for (int Y = iv.NumYVols/2; Y < iv.NumYVols - 1; ++Y)
        {
            auto &vol = iv.Volumes[Y*iv.NumXVols + X];
            vol.TopLeft = {X*dX, (Y+1)*dY};
            vol.BottomRight = {(X+1)*dX, Y*dY};
            vol.InitialValue = I;
            vol.SourceValue = 0;
            vol.K = {{kL, kR, kU, kD}};
        }
    }
    return iv;
}

int main(int argc, char *argv[])
{
    if (argc < 4)
    {
        std::cout << "Usage: <output file> <dt> <t>\n";
        return 1;
    }

    const std::string fileName(argv[1]);
    const double alpha = 1;
    const double dt = std::strtod(argv[2], nullptr);
    const double t = std::strtod(argv[3], nullptr);

    HeatEquation he(BuildInputVariables(), alpha, dt);

    int totalNumIterations = 0; 
    const int numSteps = static_cast<int>(t/dt); 
    std::cout << "Num steps=" << numSteps << "\n";

    const auto start = std::chrono::steady_clock::now();
    for (int i = 0; i < numSteps; ++i)
    {
        totalNumIterations += he.SolveOneStep();
    }
    const auto end = std::chrono::steady_clock::now();
    std::cout << "took " << 1000*std::chrono::duration<double>(end - start).count() << "ms\n";

    he.GetSolution().WriteHeatMapToFile(fileName);
    std::cout << "Avg num iterations for convergence: " << static_cast<double>(totalNumIterations)/numSteps << "\n";
}

