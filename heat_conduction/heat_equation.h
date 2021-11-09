#pragma once
#include <vector>
#include "solution.h"
#include "input_variables.h"
#include "grid.h"

struct HeatEquation 
{
    HeatEquation(const InputVariables &input, double alpha, double dt);
    Coefficients operator()(const Grid::Point& pt) const;

    int SolveOneStep();
    const Solution &GetSolution() const { return next_; }

private:
    Grid grid_;
    double alpha_;
    double dt_;
    Solution prev_;
    Solution next_;

    std::vector<double> leftCoeffs_;
    std::vector<double> rightCoeffs_;
    std::vector<double> upCoeffs_;
    std::vector<double> downCoeffs_;
    std::vector<double> constantCoeffs_;
    std::vector<double> deltaXdeltaY_;
    std::vector<double> sources_;
    std::vector<double> initialValues_;
};


