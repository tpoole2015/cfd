#pragma once

#include <vector>

struct Matrix 
{
  Matrix(int numRows, int numColumns); 

  double& operator()(int i, int j);
  double operator()(int i, int j) const;
  double *Data(); 
  const double *Data() const;

  // x <- A^{-1}b
  virtual void SolveLinear(std::vector<double> *b) const;

  Matrix &operator+=(const Matrix &rhs);
  Matrix &operator*=(double a);

  void Print() const;

  // keep track of how much work we're doing
  static int NumScalarMultiplications;
  static int NumMatrixMultiplications;
  static int NumMatrixAdditions;

  const int NumRows, NumColumns;

protected:
  std::vector<double> data_;

};

struct GeneralMatrix : public Matrix
{
    GeneralMatrix(int numRows, int numColumns);

    static GeneralMatrix Identity(int M);

    void UpdateLUFactorization();
    void SolveLinear(std::vector<double> *x) const;

private:
    // LU decomposition storage
    std::vector<double> plu_;
    std::vector<int> pivotIndicies_; // row i was interchanged with row pivotIndicies[i]
};

struct TridiagonalMatrix : public Matrix
{
    TridiagonalMatrix(int order);

    const double* GetDiagonal() const;
    double* GetDiagonal();
    const double* GetSubDiagonal() const;
    double* GetSubDiagonal();
    const double* GetSuperDiagonal() const;
    double* GetSuperDiagonal();

    void SolveLinear(std::vector<double> *x) const;

    const int Order;
};

Matrix operator+(Matrix& lhs, const Matrix &rhs);
Matrix operator*(double a, const Matrix &rhs);




