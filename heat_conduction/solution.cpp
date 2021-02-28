#include <fstream>
#include <sstream>
#include "solution.h"
#include "grid.h"

Solution::Solution(const Grid &grid)
    : gridPtr_(&grid), values_(grid.GetInitialValues())
{}

double Solution::operator()(const Grid::Index &idx) const
{
    return values_[idx.ToArrayIndex()];
}

double& Solution::operator()(const Grid::Index &idx)
{
    return values_[idx.ToArrayIndex()];
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
            const Point pt = gridPtr_->IndexToPoint(leftRight);

            const double val = this->operator()(leftRight);
            ofs << pt.X << " " << pt.Y << " " << val << "\n";
        }
        ofs << "\n";
    }

}


