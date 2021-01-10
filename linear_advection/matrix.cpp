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
    int INFO;
    plu_ = data_; // dgetrf forces us to make a copy
    dgetrf_(&M_, &M_, plu_.data(), &M_, pivots_.data(), &INFO);
    assert(INFO >= 0);
}

const double* Matrix::getLUFactorization() const 
{
    return plu_.data();
}

const int* Matrix::getPivotIndicies() const 
{
    return pivots_.data();
}


