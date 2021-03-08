#pragma once
#include <vector>
#include <array>

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

struct Grid
{
    struct Point 
    {
        int I;
        int J;

        const Grid *GridPtr;

        int64_t ToArrayIndex() const
        {
            return J*GridPtr->XDim + I;
        }

        bool InGrid() const
        {
            return I >= 0 && I < GridPtr->XDim && J >= 0 && J < GridPtr->YDim;
        }

        void MoveLeft()
        {
            --I;
        }

        void MoveRight()
        {
            ++I;
        }

        void MoveDown()
        {
            --J;
        }

        void MoveUp()
        {
            ++J;
        }

        std::pair<double, double> GetCoordinates() const
        {
            const auto arrIdx = this->ToArrayIndex();
            return {GridPtr->xValues_[arrIdx], GridPtr->yValues_[arrIdx]};
        }
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




