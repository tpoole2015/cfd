#include <cmath>
#include <cassert>
#include <iostream>
#include <fstream>
#include <vector>

#define MAXITERS 40
#define FACTOR   1.6
#define NTRY     1.6

using namespace std;

template <typename F>
bool rootBracket(F func, double *x1, double *x2) 
{
    if (*x1==*x2) assert("bad initial range in rootBracket");
    double f1=func(*x1);
    double f2=func(*x2);

    for (int i=1;i<=NTRY;++i) {
        if (f1*f2<0) return true;
        if (fabs(f1)<fabs(f2)) {
            *x1+=FACTOR*(*x1-*x2);
            f1=func(*x1);
        } else {
            *x2+=FACTOR*(*x2-*x1);
            f2=func(*x2);
        }
    }

    return false;
}

template <typename F>
double rootBisection(F func, double x1, double x2, double acc)
{
    double f=func(x1);
    double fmid=func(x2);
    if (f*fmid>=0) assert("root must be bracketed for bisection in rootBisection");
   
    // orient search so that f>0 lies at x+dx 
    double dx, rtb;
    if (f<0) {
        dx=x2-x1;
        rtb=x1;
    } else {
        dx=x1-x2;
        rtb=x2;
    }

    for (int i=1;i<=MAXITERS;++i) {
        dx*=0.5;
        const double xmid=rtb+dx;
        fmid=func(xmid);
        if (fmid<=0) rtb=xmid;
        if (fabs(dx)<acc||fmid==0.0) return rtb;
    }
    assert("too many bisections in rootBisection");
    return 0; // never get here
}

/////////////////////////////////////////////////////////

double areaRatio(double M, double gamma)
{
    const double x=2*(1+(gamma-1)*0.5*M*M)/(gamma+1);
    return pow(x,(gamma+1)*0.5/(gamma-1))/M;
}

double tempRatio(double M, double gamma)
{
    const double x=1+(gamma-1)*0.5*M*M;
    return pow(x,-1);
}

double pressureRatio(double M, double gamma)
{
    const double x=1+0.5*(gamma-1)*M*M;
    return pow(x,-gamma/(gamma-1));
}

int main()
{
    const double gamma=1.4;

    cout << "Enter p0 \n";
    double p0;
    cin >> p0;

    cout << "Enter T0 \n";
    double T0;
    cin >> T0;

    cout << "Enter M0 \n";
    double M0;
    cin >> M0;

    auto A = [](double x) -> double {
        if (x < 0) {
            return 1;
        } else if (x > 10) {
            return 0.1;
        }
        return (1-x/10.0) + x*0.1/10.0;
    };

    const double dx=0.1;
    const int N=static_cast<int>(10.0/dx);
    vector<double> x_Values(N+1,0), 
                   M_Values(N+1,0), 
                   p_Values(N+1,0),
                   T_Values(N+1,0);
    x_Values[0]=0;
    M_Values[0]=M0;
    p_Values[0]=p0;
    T_Values[0]=T0;
    for (int i=1;i<=N;++i) {
        const double x=i*dx;
        // first find M(x) (only look for subsonic solution)
        const double k=areaRatio(M0,gamma)*A(x)/A(0);
        const double M=rootBisection(
            [k,gamma](double m) -> double { return areaRatio(m,gamma) - k; },0.001,2,10e-6
        );
        const double p=p0*pressureRatio(M,gamma)/pressureRatio(M0,gamma);
        const double T=T0*tempRatio(M,gamma)/tempRatio(M0,gamma);

        x_Values[i]=x;
        M_Values[i]=M;
        p_Values[i]=p;
        T_Values[i]=T;
    }

    cout << "Enter output file\n";
    string outFile;
    cin >> outFile;
    ofstream fs(outFile);
    if (!fs.is_open()) {
        cout << "error opening file " << outFile << " for writing\n";
        return 1;
    }

    for (int i=0;i<=N;++i) {
        fs << x_Values[i] << "\t" 
           << M_Values[i] << "\t" 
           << p_Values[i] << "\t" 
           << T_Values[i] << "\n";
    }
}

