#include "grid.h"

Grid::Grid(int xDim, int yDim)
    : XDim(xDim), 
      YDim(yDim),
      xValues_(XDim*YDim, 0),
      yValues_(XDim*YDim, 0)
{}

Grid::Point Grid::GetTopLeft() const
{
    return {0, YDim-1, this};
}

Grid::Point Grid::GetBottomLeft() const
{
    return {0, 0, this};
}

Grid::Point Grid::GetTopRight() const
{
    return {XDim-1, YDim-1, this};
}

Grid::Point Grid::GetBottomRight() const
{
    return {XDim-1, 0, this};
}

