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
        Point(int row, int col, const Grid *ptr) : Row(row), Col(col), gridPtr_(ptr) {};

        const Grid &GetGrid() 
        {
            return *gridPtr_;
        }

        int64_t ToArrayIndex() const
        {
            return Row*gridPtr_->XDim + Col;
        }

        bool InGrid() const
        {
            return Col >= 0 && Col < gridPtr_->XDim && Row >= 0 && Row < gridPtr_->YDim;
        }

        void MoveLeft()
        {
            --Col;
        }

        void MoveRight()
        {
            ++Col;
        }

        void MoveDown()
        {
            --Row;
        }

        void MoveUp()
        {
            ++Row;
        }

        std::pair<double, double> GetCoordinates() const
        {
            const auto arrIdx = this->ToArrayIndex();
            return {gridPtr_->xValues_[arrIdx], gridPtr_->yValues_[arrIdx]};
        }

        int Row;
        int Col;

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

