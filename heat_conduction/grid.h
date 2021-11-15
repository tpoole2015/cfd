#pragma once
#include <vector>
#include <array>

struct Point2d
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

struct Grid
{
    struct Point 
    {
        Point(int x, int y, const Grid *ptr) : X(x), Y(y), gridPtr_(ptr) {};

        const Grid &GetGrid() 
        {
            return *gridPtr_;
        }

        int64_t ToArrayIndex() const
        {
            return Y*gridPtr_->XDim + X;
        }

        bool InGrid() const
        {
            return X >= 0 && X<gridPtr_->XDim && Y>=0 && Y<gridPtr_->YDim;
        }

        bool IsBoundary() const
        {
            return X==0||X==(gridPtr_->XDim-1)||Y==0||Y==(gridPtr_->YDim-1);
        }

        std::pair<double, double> GetCoordinates() const
        {
            const auto arrIdx = this->ToArrayIndex();
            return {gridPtr_->xValues_[arrIdx], gridPtr_->yValues_[arrIdx]};
        }

        int X;
        int Y;

    private:
        const Grid *gridPtr_;
    };

    Grid(const int xDim, const int yDim);

    Point GetTopLeft() const;
    Point GetBottomLeft() const;
    Point GetTopRight() const;
    Point GetBottomRight() const;

    const int XDim;
    const int YDim;

    std::vector<double> xValues_;
    std::vector<double> yValues_;
};

