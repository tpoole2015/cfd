#pragma once

extern "C" {
void dgetrf_(const int *M,
             const int *N,
             double *A,
             const int *LDA,
             int *IPIV,
             int *INFO);

void dgetrs_(const char *TRANS,
             const int *N,
             const int *NRHS,
             const double *A,
             const int *LDA,
             const int *IPIV,
             double *B,
             const int *LDB,
             int *INFO);
}


