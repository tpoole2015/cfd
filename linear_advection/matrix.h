#pragma once

#include <vector>

struct Matrix 
{
  Matrix(int M); 

  // used to initialize the matrix according to a finite difference scheme
  Matrix(int M, const std::vector<double> &finiteDiffMask);

  double& operator()(int i, int j);
  double operator()(int i, int j) const;
  double *data(); 
  int dim() const; 
  void updateLUFactorization();
  const double *getLUFactorization() const;
  const int *getPivotIndicies() const;

private:
  const int M_;
  std::vector<double> data_;

  // LU decomposition storage
  std::vector<double> plu_;
  std::vector<int> pivots_;
};


