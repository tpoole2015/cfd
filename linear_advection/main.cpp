#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <chrono>
#include <cassert>
#include "matrix.h"
#include "lapack_interface.h"

using namespace std;

vector<double> generateGridPoints(int M) 
{
   const double deltaX = 1.0 / static_cast<double>(M);
   vector<double> x(M);
   for (int i = 0; i < M; ++i) 
       x[i] = deltaX*static_cast<double>(i);
   return x;
}

// this function is the Fourier series of
// f(x) = 1    0.25 <= x <= 0.75
//        0    otherwise
//  we cut the Fourier series off at the harmonic frequency kMax
double initialCondition(double x, int kMax)
{
    double y = .5;
    for (int m = 0; m < (kMax - 1)/2; ++m) {
        const int k = 2*m + 1; // only do odd k
        const double a_k = (m%2 ? 2.0 : -2.0) / static_cast<double>(M_PI*k);
        y += a_k*cos(2*M_PI*k*x);
    }
    return y;
}

void writeData(const string &fn, const vector<double> &x, const vector<double> &y)
{
    ofstream fs(fn);
    if (!fs.is_open()) {
        cout << "error opening file " << fn << " for writing\n";
        return;
    }

    assert(x.size() == y.size());
    for (int i = 0; i < static_cast<int>(x.size()); ++i)
        fs << x[i] << "\t" << y[i] << "\n";
}

int main()
{
    cout << "Enter N (number of time steps)\n";
    int N;
    cin >> N; 
    assert (N >= 0);

    cout << "Enter M (number of space steps)\n";
    int M;
    cin >> M;
    assert (M >= 2);

    cout << "Enter h (dt)\n";
    double h;
    cin >> h;
    assert (h > 0);

    cout << "Enter a\n";
    double a;
    cin >> a;

    cout << "Enter beta (beta=0 centered finite difference operator, beta=-1 forward f.d operator, beta=1 backward f.d operator)\n";
    double beta;
    cin >> beta;
    assert (beta >= -1 && beta <= 1);

    cout << "Enter kMax (maximum harmonic used in building our initial condition)\n";
    int kMax;
    cin >> kMax;
  
    const double alpha=-0.5*M*a;
    const vector<double> c{-alpha*(1+beta), 2*alpha*beta, alpha*(1-beta)};
    const Matrix A(M, c);

    Matrix K(M);   // K = I - h*A
    for (int i = 0; i < A.dim(); ++i)
        for (int j = 0; j < A.dim(); ++j) {
            if (i == j) K(i, j) = 1;
            K(i,j) -= h*A(i,j);
        }
 
    K.updateLUFactorization();

    // generate our initial condition u_0(x)
    const vector<double> gridPoints = generateGridPoints(M);
    vector<double> u;
    for (double p : gridPoints)
      u.push_back(initialCondition(p, kMax));

    const auto start = chrono::steady_clock::now();

    for (int n = 0; n < N; ++n) {
        K.solveLinear(&u);
    }

    const auto end = chrono::steady_clock::now();
    cout << "took " << 1000*chrono::duration<double>(end - start).count() << "ms\n";

    cout << "Enter output file\n";
    string outFile;
    cin >> outFile;
    writeData(outFile, gridPoints, u);
}


