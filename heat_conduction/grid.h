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
        Point(int x, int y, const Grid *ptr) : x_(x), y_(y), gridPtr_(ptr) {};

        const Grid &GetGrid() 
        {
            return *gridPtr_;
        }

        int64_t ToArrayIndex() const
        {
            return y_*gridPtr_->XDim + x_;
        }

        bool InGrid() const
        {
            return x_ >= 0 && x_ < gridPtr_->XDim && y_ >= 0 && y_ < gridPtr_->YDim;
        }

        bool IsBoundary() const
        {
            return x_==0||x_==(gridPtr_->XDim-1)||y_==0||y_==(gridPtr_->YDim-1);
        }

        void MoveLeft()
        {
            --x_;
        }

        void MoveRight()
        {
            ++x_;
        }

        void MoveDown()
        {
            --y_;
        }

        void MoveUp()
        {
            ++y_;
        }

        std::pair<double, double> GetCoordinates() const
        {
            const auto arrIdx = this->ToArrayIndex();
            return {gridPtr_->xValues_[arrIdx], gridPtr_->yValues_[arrIdx]};
        }

        int X() const { return x_; }
        int Y() const { return y_; }

    private:
        int x_;
        int y_;
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

