#include "functions.hpp"

// Map internal grid indices (i,j) to a single vector index k.
//
// Grid:
//   - There are M points in each direction including boundaries.
//   - Internal points have indices i,j = 1, 2, ..., M-2.
// Ordering in the vector u^n:
//   (u_{1,1}, u_{2,1}, ..., u_{M-2,1},
//    u_{1,2}, u_{2,2}, ..., u_{M-2,2},
//    ...
//    u_{1,M-2}, ..., u_{M-2,M-2})
//
// We use 0-based indexing in C++ for the vector: k = 0, 1, ..., (M-2)^2 - 1.
inline int ij_to_k(int i, int j, int M)
{
    // Safety checks: make sure (i,j) are internal points
    if (i < 1 || i > M - 2 || j < 1 || j > M - 2) {
        throw std::out_of_range("ij_to_k: (i,j) outside internal grid");
    }

    const int N_internal = M - 2;   // number of internal points per direction

    // Convert to 0-based internal indices
    const int ii = i - 1;
    const int jj = j - 1;

    // x-index (i) runs fastest, then y-index (j)
    // k = jj * N_internal + ii
    return jj * N_internal + ii;
}

//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------

// A function that prints the structure of a sparse matrix to screen.
void print_sp_matrix_structure(const arma::sp_cx_mat& A)
{
    using namespace std;
    using namespace arma;

    // Declare a C-style 2D array of strings.
    string S[A.n_rows][A.n_cols];  

    // Initialise all the strings to " ".
    for (int i =0; i < A.n_rows; i++)
    {
        for (int j = 0; j < A.n_cols; j++)
        {
            S[i][j] = " ";
        }
    }

    // Next, we want to set the string to a dot at each non-zero element.
    // To do this we use the special loop iterator from the sp_cx_mat class
    // to help us loop over only the non-zero matrix elements.
    sp_cx_mat::const_iterator it     = A.begin();
    sp_cx_mat::const_iterator it_end = A.end();

    int nnz = 0;
    for(it; it != it_end; ++it)
    {
        S[it.row()][it.col()] = "•";
        nnz++;
    }

    // Finally, print the matrix to screen.
    cout << endl;
    for (int i =0; i < A.n_rows; i++)
    {
        cout << "| ";
        for (int j = 0; j < A.n_cols; j++)
        {
            cout << S[i][j] << " ";
        }
        cout <<  "|\n";
    }

    cout << endl;
    cout << "matrix size: " << A.n_rows << "x" << A.n_cols << endl;
    cout << "non-zero elements: " << nnz << endl ;
    cout << endl;
}

//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------

// Build sparse matrices A and B for the 2D Crank–Nicolson scheme.
//
//  - M          : total number of grid points in x (and y) including boundaries
//  - a, b       : diagonal entries (length (M-2)^2), ordered as in u^n
//  - r          : complex Crank–Nicolson parameter  r = i*Δt/(2*h^2)
//  - A, B       : output sparse matrices of size (M-2)^2 × (M-2)^2
//
// The pattern is:
//   * diagonal = a_k (for A) and b_k (for B)
//   * ±r connects each grid point to its four nearest neighbours
//     (left/right/up/down) in the 2D lattice.
void construct_AB_matrices(const arma::cx_vec& a, const arma::cx_vec& b, std::complex<double> r, int M, arma::sp_cx_mat& A, arma::sp_cx_mat& B)
{
    const int N_internal = M - 2;              // # of internal points per direction
    const int N = N_internal * N_internal;     // total # of internal points

    if ((int)a.n_elem != N || (int)b.n_elem != N) {
        throw std::runtime_error("construct_AB_matrices: wrong length of a or b.");
    }

    // Allocate (and clear) sparse matrices
    A = arma::sp_cx_mat(N, N);
    B = arma::sp_cx_mat(N, N);

    // Loop over internal grid points (i,j)
    for (int j = 1; j <= N_internal; ++j) {
        for (int i = 1; i <= N_internal; ++i) {

            // Convert (i,j) -> single index k (0-based)
            const int k = ij_to_k(i, j, M);

            // Set diagonal elements
            A(k, k) = a(k);
            B(k, k) = b(k);

            // --- Neighbour to the right: (i+1, j) ---
            if (i < N_internal) {
                const int kR = ij_to_k(i + 1, j, M);

                // Symmetric coupling with weight ±r
                A(k,  kR) = -r;
                A(kR, k ) = -r;

                B(k,  kR) =  r;
                B(kR, k ) =  r;
            }

            // --- Neighbour above: (i, j+1) ---
            if (j < N_internal) {
                const int kU = ij_to_k(i, j + 1, M);

                A(k,  kU) = -r;
                A(kU, k ) = -r;

                B(k,  kU) =  r;
                B(kU, k ) =  r;
            }
        }
    }
}

//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------

// Build the diagonal vectors a(N) and b(N) from M, h, Δt and V.
//
// M      : total number of grid points (including boundaries)
// h      : spatial step
// dt     : time step Δt
// V      : potential matrix (size M×M, real-valued)
// a, b   : output complex vectors of length N = (M-2)^2
void build_ab_vectors(int M, double h, double dt, const arma::mat& V, arma::cx_vec& a, arma::cx_vec& b)
{
    const int N_internal = M - 2;
    const int N          = N_internal * N_internal;

    a.set_size(N);
    b.set_size(N);

    const std::complex<double> I(0.0, 1.0);
    const std::complex<double> r = I * dt / (2.0 * h * h);

    // Loop over internal points (i,j) = 1..M-2
    for (int j = 1; j <= N_internal; ++j) {
        for (int i = 1; i <= N_internal; ++i) {

            const int k = ij_to_k(i, j, M);   // 0-based index in [0, N-1]

            const double vij = V(i, j);      // internal potential value

            const std::complex<double> vterm = I * (dt * 0.5) * vij;

            a(k) = 1.0 + 4.0 * r + vterm;
            b(k) = 1.0 - 4.0 * r - vterm;
        }
    }
}
