#include <fstream>
#include <sstream>
#include "solution.h"
#include "grid.h"

Solution::Solution(const Grid *gridPtr)
    : gridPtr_(gridPtr)
{}

void Solution::Init(const std::vector<double> &values)
{
    values_ = values;
}

double Solution::operator()(const Grid::Point &pt) const
{
    return values_[pt.ToArrayIndex()];
}

double& Solution::operator()(const Grid::Point &pt)
{
    return values_[pt.ToArrayIndex()];
}

void Solution::WriteHeatMapToFile(const std::string &file) const
{
    std::ofstream ofs(file);
    if (!ofs.is_open())
    {
        std::stringstream errMsg;
        errMsg << "error opening " << file;
        throw std::runtime_error(errMsg.str());
    }

    for (auto topBottom = gridPtr_->GetTopLeft(); topBottom.InGrid(); topBottom.MoveDown())
    {
        for (auto leftRight = topBottom; leftRight.InGrid(); leftRight.MoveRight())
        {
            const auto pt = leftRight.GetCoordinates();
            const double val = this->operator()(leftRight);
            ofs << pt.first << " " << pt.second << " " << val << "\n";
        }
        ofs << "\n";
    }

}


