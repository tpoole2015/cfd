#pragma once
#include <vector>
#include <array>

struct Point
{
    double X;
    double Y;
};

struct InputVariables
{
    struct ControlVolume
    {
        enum BoundaryValueType { Scalar, Flux };

        Point TopLeft;
        Point BottomRight;

        double InitialValue;
        double SourceValue;
        std::array<double, 4> K; // 0 = Left face, 1 = Right face, 2 = Up face, 3 = Down face

        BoundaryValueType BVType; // Only used for boundary volumes

        // If type is Scalar: BoundaryValues = { value at grid point }
        // If type is Flux:
        //  If position is TopLeft: BoundaryValues = { flux top, flux left }
        //  If position is TopRight: BoundaryValues = { flux top, flux right }
        //  If position is BottomLeft: BoundaryValues = { flux bottom, flux left}
        //  If position is BottomRight: BoundaryValues = { flux bottom, flux right}
        //  If position is Left: BoundaryValues = { flux left }
        //  If position is Right: BoundaryValues = { flux right }
        //  If position is Bottom: BoundaryValues = { flux bottom }
        //  If position is Top: BoundaryValues = { flux top }
        // Here position is determined by location of our ControlVolume in the sorted array Volumes
        std::vector<double> BoundaryValues;
    };
    int NumXVols; 
    int NumYVols;
    std::vector<ControlVolume> Volumes; // [Y*NumXVols + X]
};

struct Solution;
struct Grid
{
    struct Index
    {
        int XIndex;
        int YIndex;
        const Grid *GridPtr;

        int ToArrayIndex() const
        {
            return YIndex*GridPtr->NumXValues + XIndex;
        }

        bool InGrid() const
        {
            return XIndex >= 0 && XIndex < GridPtr->NumXValues && YIndex >= 0 && YIndex < GridPtr->NumYValues;
        }

        void MoveLeft()
        {
            --XIndex;
        }

        void MoveRight()
        {
            ++XIndex;
        }

        void MoveDown()
        {
            --YIndex;
        }

        void MoveUp()
        {
            ++YIndex;
        }
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

    Grid(InputVariables input);

    Point IndexToPoint(const Index &i) const;

    Coefficients GetDiscretizationCoeffs(double alpha, double dt, const Index &idx, const Solution &soln) const;

    Index GetTopLeft() const;
    Index GetBottomLeft() const;
    Index GetTopRight() const;
    Index GetBottomRight() const;

    const std::vector<double> &GetInitialValues() const;

    const int NumXValues;
    const int NumYValues;

private:
    // we can cache at construction time some of the parameters used in the discretization coefficients
    std::vector<double> xValues_;
    std::vector<double> yValues_;
    std::vector<double> leftCoeffs_;
    std::vector<double> rightCoeffs_;
    std::vector<double> upCoeffs_;
    std::vector<double> downCoeffs_;
    std::vector<double> constantCoeffs_;
    std::vector<double> deltaXdeltaY_;
    std::vector<double> sources_;
    std::vector<double> initialValues_;
};

