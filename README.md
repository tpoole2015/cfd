# cfd
A program to solve the 2d heat equation using the algorithm described in Chapter 4 of "Numerical Heat Transfer and Fluid Flow" by Patankar. The main source code is in heat_equation/main.cpp. For this program I assume we're solving the heat equation (with uniform conductivity k) on the unit rectangle [0,1]x[0,1] with

    I(x,y) = constant I
    S(x,y) = S in the bottom half of the grid
             0 in the top half

Here I and S are both specified via the command prompt. Also the top, bottom, left and right boundary values are constant and also specified via the command prompt. The program dumps the final solution in a text file. To plot the solution using gnuplot use the command

    plot '<output file>' with image
    
    Usage: ./heat_equation <output file> <t> <dx> <dy> <dt> <S> <I> <bv top> <bv bottom> <bv left> <bv right>

To build
1. Check out the source code
2. make a new directory where you want to build the code and cd into that directory 
3. cmake <directory where code was cloned into>
4. make 

The core algorithm is implemented in the code heat_equation/solver.h by the routine

    <template class T> 
    int Solve(const T &coeffs, const Grid &grid, Solution &soln);

This routine solves the discretization equation

    a_C*T_p = a_N*T_N + a_S*T_S + a_W*T_W + a_E*T_E + b 

for any grid point p where E,N,S and W denote the grid points east,north,south and west of p. Here a_C,a_N,a_S,a_W and a_E denote the center,north,south,west and east cofficients respectively. In code these are represented as

    T_p is soln(p) 
    a_C,a_N,a_S,a_W,a_E,b are represented by the coeffs(p). 

In particular coeffs(p) returns a **struct Coefficients** (defined in heat_equation/grid.h). In our case the coeffs object is given to us by the **struct HeatEquation** defined in heat_equation/heat_equation.h. The purpose of this class is to take an object representing the initial conditions (**struct InputVariables**) and generate the cofficients used in the discrete form of the equations. All matrix calculations are carried out using the lapack routines (see heat_equation/matrix.h and heat_equation/matrix.cpp).

The main grid in our case is a 2-dim rectangular grid, but by pulling it out and putting it in its own separate class hopefully the generalization to other grids is not too hard.


