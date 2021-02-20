// Chapter 4 of Numerical Heat Transfer and Fluid Flow by Patankar


pair<Matrix, vector<double>> BuildTriDiagonalEquations(const Grid &g, const Grid::Point &pt, const ScalarField &sf) 
{
    const int dim = grid.YSize(pt.X);

    // We construct a banded matrix as our system is in tridiagonal form.
    // See http://www.netlib.org/lapack/lug/node124.html
    const int kl = 1; // num subdiagonals
    const int ku = 1; // num superdiagonals
    Matrix m(kl+ku+1, dim, Matrix::Banded);

    vector<double> rhs(dim);

    int i = 0;
    for ( ; grid.HasNorth(pt); pt = grid.IncNorth(pt))
    {
        const Coefficients eqnCoefs = basePoint.GetDiscretizationCoeffs(sf);

        m(1, i) = eqnCoefs.Center; // diagonal
        if (i > 0)
        {
            banded(2,i) = eqnCoefs.South; // subdiagonal
        }
        if (i < dim)
        {
            banded(0,i) = eqnCoefs.North; // superdiagonal
        }
        double c = eqnCoefs.Constant;
        if (grid.HasWest(pt))
        {
            c += eqnCoefs.West*sf(grid.IncWest(pt));
        }
        if (grid.HasEast(pt))
        {
            c += eqnCoefs.East*sf(grid.IncEast(pt));
        }
        rhs[i] = c;

        ++i;
    }

    return {m, c}; // Is m copied here? (need to check)
}

int main(int argc, char *argv[])
{
    Grid grid(xSize, ySize);

    ScalarField cur(grid), next(grid); // i.e Matrix(grid.XDim(), grid.YDim());

    PopulateBoundary(cur);
    for (int t = 0; t < numTimeSteps; ++t)
    {
        for (auto pt = grid.GetOrigin(); grid.HasEast(pt); pt = grid.IncEast(pt))
        {
            const auto triDiag = BuildTriDiagonalMatrix(grid, pt, cur);
            // TODO: solve tridiagonal matrix
        }
        cur = std::move(next);
    }    
}

