#pragma once

#include <vector>

template <template<typename> class Allocator>
struct Matrix 
{
  Matrix(int M); 

  // used to initialize the matrix according to a finite difference scheme
  Matrix(int M, const vector<double> &c);

  double& operator()(int i, int j);
  double operator()(int i, int j) const;
  double *data(); 
  int dim() const; 
  void updateLUFactorization();
  const double *getLUFactorization() const;
  const int *getPivotIndicies() const;

private:
  static std::vector<std::pair<int, double>> finiteDifferenceCoeffs(int i, int M, const std::vector<double> &c);

  const int M_;
  std::vector<double, Allocator<double>> data_;

  // LU decomposition storage
  std::vector<double, Allocator<double>> plu_;
  std::vector<int, Allocator<int>> pivots_;
};


