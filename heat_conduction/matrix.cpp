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

int Matrix::NumScalarMultiplications = 0;
int Matrix::NumMatrixMultiplications = 0;
int Matrix::NumMatrixAdditions = 0;

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

const double* Matrix::data() const
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
    Matrix m(M, {1});
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

// A*=B => A<-A*B
Matrix& Matrix::operator*=(const Matrix &rhs)
{
    // unorthodox but forced on us by the structure of the BLAS routines
    const Matrix result = (*this) * rhs;
    this->data_ = result.data_;
    return *this;
}

Matrix Matrix::inverse() const 
{
    Matrix inv(*this);

    // first compute LU decomposition
    std::vector<int> ipiv(N);
    int info;
    dgetrf_(&M_, &M_, inv.data(), &M_, ipiv.data(), &info);
    assert(info == 0); 

    // now computer the inverse
    std::vector<double> work(M_);
    dgetri_(&M_, inv.data(), &M_, ipiv.data(), work.data(), &M_, &info);
    assert(info == 0);

    return inv;
}

Matrix operator+(Matrix lhs, const Matrix &rhs)
{
    return lhs += rhs;
}

Matrix operator*(double a, const Matrix &rhs)
{
    const int M = rhs.dim();
    Matrix out(M);
    const int n = M*M;
    const int incx = 1;
    const int incy = 1;
    daxpy_(&n, &a, rhs.data(), &incx, out.data(), &incy);

    Matrix::NumScalarMultiplications += 1;
    return out;
}

Matrix operator*(const Matrix &lhs, const Matrix &rhs)
{
    assert(lhs.dim()==rhs.dim());
    const int M = lhs.dim();
    Matrix out(M);
    const char transa = 'n';
    const char transb = 'n';
    const double alpha = 1;
    const double beta = 0;
    dgemm_(&transa, &transb, &M, &M, &M, &alpha, lhs.data(), &M, rhs.data(), &M, &beta, out.data(), &M);

    Matrix::NumMatrixMultiplications += 1;
    return out;
}

