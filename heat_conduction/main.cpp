// Chapter 4 of Numerical Heat Transfer and Fluid Flow by Patankar
#include <iostream>
#include <cmath>
#include <chrono>
#include <cstdlib>
#include "input_variables.h"
#include "heat_equation.h"

InputVariables BuildInputVariables()
{
    const double kL = 0, kR = 0, kU = 10, kD = 10;
    const double I = 0;
    const double S = 100;
    const double bvTop = 0, bvBottom = 0, bvLeft = 0, bvRight = 0; 
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
        tvol.K = {{kL, kR, kU, kD}};
        tvol.BoundaryValue = bvTop;
 
        // Bottom
        const int bottomY = 0;
        auto &bvol = iv.Volumes[bottomY*iv.NumXVols + X];
        bvol.TopLeft = {X*volumeWidth, (bottomY + 1)*volumeHeight};
        bvol.BottomRight = {(X+1)*volumeWidth, bottomY*volumeHeight};
        bvol.InitialValue = 0;
        bvol.K = {{kL, kR, kU, kD}};
        bvol.BoundaryValue = bvBottom;
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
        lvol.K = {{kL, kR, kU, kD}};
        lvol.BoundaryValue = bvLeft;

        // Right 
        const int rightX = iv.NumXVols-1;
        auto &rvol = iv.Volumes[Y*iv.NumXVols + rightX];
        rvol.TopLeft = {(rightX-1)*volumeWidth, (Y+1)*volumeHeight};
        rvol.BottomRight = {rightX*volumeWidth, Y*volumeHeight};
        rvol.InitialValue = 0;
        rvol.K = {{kL, kR, kU, kD}};
        rvol.BoundaryValue = bvRight;
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
            vol.SourceValue = S;
            vol.K = {{kL, kR, kU, kD}};
        }
        // don't heat up the top half
        for (int Y = iv.NumYVols/2; Y < iv.NumYVols - 1; ++Y)
        {
            auto &vol = iv.Volumes[Y*iv.NumXVols + X];
            vol.TopLeft = {X*volumeWidth, (Y+1)*volumeHeight};
            vol.BottomRight = {(X+1)*volumeWidth, Y*volumeHeight};
            vol.InitialValue = I;
            vol.SourceValue = 0;
            vol.K = {{kL, kR, kU, kD}};
        }
    }
    return iv;
}

int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        std::cout << "Usage: <output file> <time>\n";
        return 1;
    }

    const std::string fileName(argv[1]);
    const double alpha = 1;
    const double dt = 0.1;
    const double t = std::strtod(argv[2], nullptr);

    HeatEquation he(BuildInputVariables(), alpha, dt);

    int totalNumIterations = 0; 
    const int numSteps = static_cast<int>(t/dt); 
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

