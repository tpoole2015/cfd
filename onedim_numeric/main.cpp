#include <vector>
#include "matrix.h"

using namespace std;
constexpr double gamma=1.4; // normally globals are bad but I think we can make an exception here

/////// Helper functions /////// 
double p(int i, const vector<double>& Q);
double pho(int i, const vector<double>& Q);
double u(int i, const vector<double>& Q);
double c(int i, const vector<double>& Q);
////////////////////////////////

Matrix fluxJacobian(int i, const vector<double>& Q)
{
    const double Q1 = Q[3*i];
    const double Q2 = Q[3*i+1];
    const double Q3 = Q[3*i+2];

    Matrix A(3);
    A(0,0) = 0;
    A(0,1) = 1;
    A(0,2) = 0;

    A(1,0) = (gamma-3)*0.5*(Q2*Q2)/(Q1*Q1);
    A(1,1) = (3-gamma)*Q2/Q1;
    A(1,2) = gamma-1;

    A(2,0) = ((gamma-1)*(Q2*Q2*Q2)/(Q1*Q1*Q1)) - (gamma*(Q3*Q2)/(Q1*Q1));
    A(2,1) = (gamma*Q3/Q1) - (3*(gamma-1)*0.5*(Q2*Q2)/(Q1*Q1));
    A(2,2) = gamma*Q2/Q1;

    return A;
}

Matrix phi_1(int i,
             double h,
             int I,
             double artificalDissipationConst,
             const vector<double>& Q)
{
    assert (i >= 0 && i <= I);
    if (i >= 1 && i <= I-1) {
        const double X = artificalDissipationConst*0.5*(u(i,Q) + c(i,Q) + u(i+1,Q) + c(i+1,Q));
        return (1 + 3*X*h)*Matrix::identity(3) + h*dSdX[i]*fluxJacobian(i,Q);
    } else if (i == I) {
    }
    return Matrix(3); // not needed for LHS grid point
}

Matrix phi_2(int i,
             double h,
             double dX,
             int I,
             vector<double> dSdX,
             double artificalDissipationConst,
             const Matrix &boundaryMatrix,
             const vector<double>& Q)
{
    assert (i >= 0 && i <= I);
    if (i >= 1 && i <= I-1) {
        const double X = artificalDissipationConst*0.5*(u(i,Q) + c(i,Q) + u(i+1,Q) + c(i+1,Q));
        return (-h/dX)*fluxJacobian(i-1,Q) - h*X*Matrix::identity(3);
    } else if (i == 0) {
        const double X = artificalDissipationConst*0.5*(u(0,Q) + c(0,Q) + u(1,Q) + c(1,Q));
        return Matrix::identity(3) + ((-h/dX) + h*SdX[0])*boundaryMatrix*fluxJacobian(0,Q) + (3*X*h)*Matrix::identity(3);
    }
    // TODO: i=I
}

Matrix phi_3(int i,
             double h,
             double dX,
             int I,
             double artificalDissipationConst,
             const Matrix &boundaryMatrix,
             const vector<double>& Q)
{
    assert (i >= 0 && i <= I);
    if (i <= I-1) {
        const double X = artificalDissipationConst*0.5*(u(i,Q) + c(i,Q) + u(i+1,Q) + c(i+1,Q));
        return (h/dX)*fluxJacobian(i+1,Q) - 3*h*X*Matrix::identity(3);
    }
    // TODO: i=I
}

Matrix phi_4(int i,
             double h,
             double dX,
             double artificalDisapationConst,
             const vector<double>& Q)
{
    assert (i >= 0 && i <= I);
    if (i <= I-1) {
        const double X = artificalDisapationConst*0.5*(u(i,Q) + c(i,Q) + u(i+1,Q) + c(i+1,Q));
        return h*X*Matrix::identity(3);
    }
    // TODO: i=I
}

