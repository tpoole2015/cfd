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

  // x <- A^{-1}x
  void solveLinear(std::vector<double> *x) const;

  Matrix &operator+=(const Matrix &rhs);
  Matrix &operator*=(double a);
  Matrix &operator*=(const Matrix &rhs);

  static Matrix identity(int M);
private:
  const int M_;
  std::vector<double> data_;

  // LU decomposition storage
  std::vector<double> plu_;
  std::vector<int> pivots_;
};

Matrix operator+(Matrix lhs, const Matrix &rhs);
Matrix operator*(double a, Matrix rhs);
Matrix operator*(Matrx lhs, const Matrix &rhs);



