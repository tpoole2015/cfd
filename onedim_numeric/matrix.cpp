#include <cassert>
#include "matrix.h"
#include "lapack_interface.h"

using namespace std;
namespace {
vector<pair<int, double>> finiteDifferenceCoeffs(int i, int M, const vector<double> &c)
{
    assert(c.size()%2);  // need c to have odd size
    vector<pair<int, double>> mask;
    const int k=c.size()/2;
    auto cIt = c.begin();
    for (int I=i-k; I<=i+k; ++I) {
        const int gridIndex=(I < 0)?(M + I%M)%M:I%M;
        mask.push_back({gridIndex, *cIt++});
    }
    return mask;
}
}

Matrix::Matrix(int M) : 
    M_(M), 
    data_(M_*M_, 0),
    plu_(M_*M_, 0),
    pivots_(M_, 0)
{
}
 
Matrix::Matrix(int M, const vector<double> &finiteDiffMask) : Matrix(M) 
{
    for (int i = 0; i < M_; ++i) {
        for (const auto &m : finiteDifferenceCoeffs(i, M_, finiteDiffMask))
            (*this)(i,m.first) = m.second; 
    }
}

double& Matrix::operator()(int i, int j) 
{
    return data_.data()[j*M_ + i];  // fortran style (column-major order)
}

double Matrix::operator()(int i, int j) const 
{
    return data_.data()[j*M_ + i];  // fortran style (column-major order)
}

double* Matrix::data() 
{
    return data_.data();
}
 
int Matrix::dim() const 
{
    return M_;
}

void Matrix::updateLUFactorization() 
{
    int info;
    plu_ = data_; // dgetrf forces us to make a copy
    dgetrf_(&M_, &M_, plu_.data(), &M_, pivots_.data(), &info);
    assert(info >= 0);
}

void Matrix::solveLinear(vector<double> *x) const
{
    assert(static_cast<int>(x->size()) == M_);
    const int nrhs = 1; 
    const char trans = 'n';
    int info;
    dgetrs_(&trans, &M_, &nrhs, plu_.data(), &M_, pivots_.data(), x->data(), &M_, &info);
    assert(info == 0);
}

Matrix Matrix::identity(int M)
{
    Matrix m(M, {1}); // should work, double check
    return m;
}
 
Matrix& Matrix::operator+=(const Matrix &rhs)
{
    assert (M_ == rhs.dim());
    const int n = M_*M_;
    const double alpha = 1;
    const int incx = 1;
    const int incy = 1;
    daxpy_(&n, &alpha, rhs.data(), &incx, this->data(), &incy);
    return *this;
}

Matrix& Matrix::operator*=(double a) 
{
    // unorthodox but forced on us by the structure of the BLAS routines
    const Matrix result = a * (*this);
    this->data_ = result.data_;
    return *this;
}

// A*=B => A<-A*B
Matrix& Matrix::operator*=(const Matrix &rhs);
{
    // unorthodox but forced on us by the structure of the BLAS routines
    const Matrix result = (*this) * rhs;
    this->data_ = result.data_;
    return *this;
}

Matrix operator+(Matrix lhs, const Matrix &rhs)
{
    return lhs += rhs;
}

Matrix operator*(double a, const Matrix &rhs)
{
    Matrix out(M_);
    const int n = M_*M_;
    const int incx = 1;
    const int incy = 1;
    daxpy_(&n, &a, rhs.data(), &incx, out.data(), &incy);
    return out;
}

Matrix operator*(const Matrix &lhs, const Matrix &rhs)
{
    Matrix out(M_);
    const char transa = 'n';
    const char transb = 'n';
    const int alpha = 1;
    const int beta = 0;
    dgemm_(&transa, &transb, &M_, &M_, &M_, &alpha, lhs.data(), &M_, rhs.data(), &M_, &beta, out.data(), &M_);
    return out;
}

