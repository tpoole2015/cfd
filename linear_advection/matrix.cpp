#include "matrix.h"

using namespace std;

template <template<typename> class Allocator>
Matrix<Allocator>::Matrix(int M) : 
    M_(M), 
    data_(M_*M_, 0, Allocator<double>()),
    plu_(M_*M_, 0, Allocator<double>()),
    pivots_(M_, 0, Allocator<int>()) 
{
}
 
template <template<typename> class Allocator>
Matrix<Allocator>::Matrix(int M, const vector<double> &c) : 
    Matrix(M) 
{
    for (int i = 0; i < M_; ++i) {
        for (const auto &m : finiteDifferenceCoeffs(i, M_, c))
            (*this)(i,m.first) = m.second; 
    }
}

template <template<typename> class Allocator>
double& Matrix<Allocator>::operator()(int i, int j) 
{
    return this->data()[j*M_ + i];  // fortran style (column-major order)
}

template <template<typename> class Allocator>
double Matrix<Allocator>::operator()(int i, int j) const 
{
    return this->data()[j*M_ + i];  // fortran style (column-major order)
}

template <template<typename> class Allocator>
double* Matrix<Allocator>::data() 
{
    return data_.data();
}
 
template <template<typename> class Allocator>
int Matrix<Allocator>::dim() const 
{
    return M_;
}

template <template<typename> class Allocator>
void Matrix<Allocator>::updateLUFactorization() 
{
    int INFO;
    auto plu_ = data_; // dgetrf forces us to make a copy
    dgetrf_(&amp;M_, &amp;M_, plu_.data(), &amp;M_, pivots_.data(), &amp;INFO);
    assert(INFO >= 0);
}

template <template<typename> class Allocator>
const double* Matrix<Allocator>::getLUFactorization() const 
{
    return plu_.data();
}

template <template<typename> class Allocator>
const int* Matrix<Allocator>::getPivotIndicies() const 
{
    return pivots_.data();
}

template <template<typename> class Allocator>
vector<pair<int, double>> Matrix<Allocator>::finiteDifferenceCoeffs(int i, int M, const vector<double> &c)
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


