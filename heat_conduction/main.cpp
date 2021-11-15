// Chapter 4 of Numerical Heat Transfer and Fluid Flow by Patankar
#include <iostream>
#include <cmath>
#include <chrono>
#include <cstdlib>
#include "input_variables.h"
#include "heat_equation.h"

InputVariables BuildInputVariables(double dx, double dy, double I, double S, double bvTop, double bvBottom, double bvLeft, double bvRight)
{
    const double kL = 1, kR = 1, kU = 1, kD = 1;
    InputVariables iv;
    iv.NumXVols = 1.0 / dx;
    iv.NumYVols = 1.0 / dy;
    iv.Volumes.resize(iv.NumXVols*iv.NumYVols);

    // Top & Bottom
    for (int X = 0; X < iv.NumXVols; ++X)
    {
        // Top 
        const int topY = iv.NumYVols-1;
        auto &tvol = iv.Volumes[topY*iv.NumXVols + X];
        tvol.TopLeft = {X*dx, (topY+1)*dy};
        tvol.BottomRight = {(X+1)*dx, topY*dy};
        tvol.InitialValue = 0;
        tvol.SourceValue = 0;
        tvol.K = {{kL, kR, kU, kD}};
        tvol.BoundaryValue = bvTop;
 
        // Bottom
        const int bottomY = 0;
        auto &bvol = iv.Volumes[bottomY*iv.NumXVols + X];
        bvol.TopLeft = {X*dx, (bottomY + 1)*dy};
        bvol.BottomRight = {(X+1)*dx, bottomY*dy};
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
        lvol.TopLeft = {leftX*dx, (Y+1)*dy};
        lvol.BottomRight = {(leftX+1)*dx, Y*dy};
        lvol.InitialValue = 0;
        lvol.SourceValue = 0;
        lvol.K = {{kL, kR, kU, kD}};
        lvol.BoundaryValue = bvLeft;

        // Right 
        const int rightX = iv.NumXVols-1;
        auto &rvol = iv.Volumes[Y*iv.NumXVols + rightX];
        rvol.TopLeft = {(rightX-1)*dx, (Y+1)*dy};
        rvol.BottomRight = {rightX*dx, Y*dy};
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
            vol.TopLeft = {X*dx, (Y+1)*dy};
            vol.BottomRight = {(X+1)*dx, Y*dy};
            vol.InitialValue = I;
            vol.SourceValue = S;
            vol.K = {{kL, kR, kU, kD}};
        }
        // don't heat up the top half
        for (int Y = iv.NumYVols/2; Y < iv.NumYVols - 1; ++Y)
        {
            auto &vol = iv.Volumes[Y*iv.NumXVols + X];
            vol.TopLeft = {X*dx, (Y+1)*dy};
            vol.BottomRight = {(X+1)*dx, Y*dy};
            vol.InitialValue = I;
            vol.SourceValue = 0;
            vol.K = {{kL, kR, kU, kD}};
        }
    }
    return iv;
}

int main(int argc, char *argv[])
{
    if (argc < 12)
    {
        std::cout << "Usage: <output file> <t> <dx> <dy> <dt> <S> <I> <bv top> <bv bottom> <bv left> <bv right>\n";
        return 1;
    }

    const std::string fileName(argv[1]);
    const double alpha = 1;
    const double t = std::strtod(argv[2], nullptr);
    const double dx = std::strtod(argv[3], nullptr);
    const double dy = std::strtod(argv[4], nullptr);
    const double dt = std::strtod(argv[5], nullptr);
    const double S = std::strtod(argv[6], nullptr) / (dx*dy);
    const double I = std::strtod(argv[7], nullptr);
    const double bvTop = std::strtod(argv[8], nullptr);
    const double bvBottom = std::strtod(argv[9], nullptr);
    const double bvLeft = std::strtod(argv[10], nullptr);
    const double bvRight = std::strtod(argv[11], nullptr);

    HeatEquation he(BuildInputVariables(dx, dy, I, S, bvTop, bvBottom, bvLeft, bvRight), alpha, dt);

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

