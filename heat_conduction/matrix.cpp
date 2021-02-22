#include <cassert>
#include "matrix.h"
#include "lapack_interface.h"

using namespace std;

int Matrix::NumScalarMultiplications = 0;
int Matrix::NumMatrixMultiplications = 0;
int Matrix::NumMatrixAdditions = 0;

Matrix::Matrix(int numRows, int numColumns) : 
    NumRows(numRows), 
    NumColumns(numColumns), 
    data_(numRows*numColumns, 0)
{
}
 
double& Matrix::operator()(int i, int j) 
{
    return data_.data()[j*NumRows + i];  // fortran style (column-major order)
}

double Matrix::operator()(int i, int j) const 
{
    return data_.data()[j*NumRows + i];  // fortran style (column-major order)
}

double* Matrix::Data() 
{
    return data_.data();
}

const double* Matrix::Data() const
{
    return data_.data();
}

void Matrix::SolveLinear(std::vector<double> *b) const
{
    (void)(b);
    assert(false);
}

Matrix& Matrix::operator+=(const Matrix &rhs)
{
    assert (NumRows == rhs.NumRows);
    assert (NumColumns == rhs.NumColumns);

    const int n = NumRows*NumColumns;
    const double alpha = 1;
    const int incx = 1;
    const int incy = 1;
    daxpy_(&n, &alpha, rhs.Data(), &incx, this->Data(), &incy);

    Matrix::NumMatrixAdditions += 1;
    return *this;
}

Matrix& Matrix::operator*=(double a) 
{
    // unorthodox but forced on us by the structure of the BLAS routines
    const Matrix result = a * (*this);
    this->data_ = result.data_;
    return *this;
}

//////////////// GeneralMatrix //////////////// 
GeneralMatrix::GeneralMatrix(int numRows, int numColumns) :
    Matrix(numRows, numColumns)
{
}

GeneralMatrix GeneralMatrix::Identity(int M)
{
    GeneralMatrix m(M, M);
    for (int i = 0; i < M; ++i)
    {
        m(i,i) = 1;
    }
    return m;
}

void GeneralMatrix::UpdateLUFactorization() 
{
    int info;
    pivotIndicies_.reserve(min(NumRows, NumColumns));
    plu_ = data_; // dgetrf forces us to make a copy
    dgetrf_(&NumRows, &NumColumns, plu_.data(), &NumRows, pivotIndicies_.data(), &info);
    assert(info >= 0);
}

void GeneralMatrix::SolveLinear(vector<double> *b) const
{
    assert(NumRows == NumColumns);
    assert(static_cast<int>(b->size()) == NumRows);

    const int nrhs = 1; 
    const char trans = 'n';
    int info;
    dgetrs_(&trans, &NumRows, &nrhs, plu_.data(), &NumRows, pivotIndicies_.data(), b->data(), &NumRows, &info);
    assert(info == 0);
}

//////////////// BandedMatrix //////////////// 
TridiagonalMatrix::TridiagonalMatrix(int order) :
    Matrix(order, 3), // col 0 = diagonal, col 1 = sub diagonal, col 2 = super diagonal
    Order(order)
{
}

const double* TridiagonalMatrix::GetDiagonal() const
{
    return this->Data();
}

double* TridiagonalMatrix::GetDiagonal()
{
    return this->Data();
}

const double* TridiagonalMatrix::GetSubDiagonal() const
{
    return this->Data() + Order;
}

double* TridiagonalMatrix::GetSubDiagonal()
{
    return this->Data() + Order;
}

const double* TridiagonalMatrix::GetSuperDiagonal() const
{
    return this->Data() + 2*Order;
}

double* TridiagonalMatrix::GetSuperDiagonal()
{
    return this->Data() + 2*Order;
}

void TridiagonalMatrix::SolveLinear(vector<double> *b) const
{
    assert(static_cast<int>(b->size()) == Order);

    TridiagonalMatrix copy(*this);
    const int nrhs = 1; 
    int info;
    dgtsv_(&Order, 
           &nrhs, 
           copy.GetSubDiagonal(), 
           copy.GetDiagonal(), 
           copy.GetSuperDiagonal(), 
           b->data(), 
           &Order, 
           &info); 
    assert(info == 0);
}


//////////////// friends //////////////// 
Matrix operator+(Matrix lhs, const Matrix &rhs)
{
    return lhs += rhs;
}

Matrix operator*(double a, const Matrix &rhs)
{
    Matrix out(rhs.NumRows, rhs.NumColumns);
    const int n = rhs.NumRows * rhs.NumColumns;
    const int incx = 1;
    const int incy = 1;
    daxpy_(&n, &a, rhs.Data(), &incx, out.Data(), &incy);

    Matrix::NumScalarMultiplications += 1;
    return out;
}

