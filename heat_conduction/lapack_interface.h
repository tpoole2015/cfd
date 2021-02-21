#pragma once

extern "C" {
void dgetrf_(const int *M,
             const int *N,
             double *A,
             const int *LDA,
             int *IPIV,
             int *INFO);

void dgetri_(const int *N,
             double *A,
             int *LDA,
             const int *IPIV,
             double *WORK,
             const int *LWORK,
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

void daxpy_(const int *N,
            const double *CA,
            const double *CX,
            const int *INCX,
            double *CY,
            const int *INCY);

void dgemm_(const char *TRANSA,
            const char *TRANSB,
            const int *M,
            const int *N,
            const int *K,
            const double *ALPHA,
            const double *A,
            const int *LDA,
            const double *B,
            const int *LDB,
            const double *BETA,
            double *C,
            const int *LDC);

void dgtsv_(const int *N,
            const int *NRHS,
            double *DL,
            double *D,
            double *DU,
            double *B,
            int *LDB,
            int *INFO);
}


