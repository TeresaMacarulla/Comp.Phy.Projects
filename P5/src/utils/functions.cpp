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
void build_ab_vectors(int M, double dt, const arma::mat& V, arma::cx_vec& a, arma::cx_vec& b)
{
    const int N_internal = M - 2;
    const int N          = N_internal * N_internal;
    const double h = 1.0 / (M - 1);     // grid spacing in [0,1]

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

//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------

// Perform the matrix multiplication B·u_n=b
void matrix_mult(const arma::sp_cx_mat& B, const arma::cx_vec& u_n, arma::cx_vec& b)
{
    // Optional safety check
    if (B.n_cols != u_n.n_rows) {
        throw std::runtime_error("matrix_mult: size mismatch between B and u_n");
    }

    b = B * u_n;
}

//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------

// Perform one Crank–Nicolson time step:
// Given current u, compute u_next from A*u_next = B*u.
//
// A, B   : sparse complex matrices
// u      : current solution vector (at time step n)
// u_next : output vector (at time step n+1)
//
// Returns true if spsolve succeeded, false otherwise.
bool cn_step(const arma::sp_cx_mat& A, const arma::sp_cx_mat& B, const arma::cx_vec& u, arma::cx_vec& u_next)
{
    // Basic size checks
    if (B.n_cols != u.n_rows) {
        throw std::runtime_error("cn_step: size mismatch between B and u");
    }
    if (A.n_rows != A.n_cols || A.n_rows != B.n_rows) {
        throw std::runtime_error("cn_step: incompatible sizes for A and B");
    }

    // 1) Compute b = B * u (right-hand side)
    arma::cx_vec b;
    matrix_mult(B, u, b);

    // 2) Solve A * u_next = b using sparse solver (SuperLU)
    bool ok = arma::spsolve(u_next, A, b, "superlu");

    return ok;
}

//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------

// Set up initial state u(x,y,t=0) as a Gaussian normalised wave packet
//
// Domain: x,y ∈ [0,1], uniform grid with M points in each direction.
// Dirichlet boundaries: u = 0 at i=0, i=M-1, j=0, j=M-1.
//
// Parameters:
//   M         : number of grid points in x and y (including boundaries)
//   xc, yc    : centre of the packet
//   sigma_x   : width in x
//   sigma_y   : width in y
//   p_x, p_y  : momenta in x and y
//
// Output:
//   U0        : complex M×M matrix with U0(i,j) = u_ij^0
//
void init_gaussian_packet(arma::cx_mat& U0, int M, double xc, double yc, double sigma_x, double sigma_y, double p_x, double p_y)
{
    U0.set_size(M, M);
    U0.zeros();

    const std::complex<double> I(0.0, 1.0);
    const double h = 1.0 / (M - 1);     // grid spacing in [0,1]

    double norm = 0.0;                  // to accumulate Σ |u_ij|^2

    // Loop over all grid points
    for (int j = 0; j < M; ++j) {
        for (int i = 0; i < M; ++i) {

            // Enforce Dirichlet boundary conditions explicitly
            if (i == 0 || i == M-1 || j == 0 || j == M-1) {
                U0(i,j) = std::complex<double>(0.0, 0.0);
                continue;
            }

            double x = i * h;
            double y = j * h;

            // Unnormalised Gaussian envelope
            double exponent =
                - (x - xc)*(x - xc) / (2.0 * sigma_x * sigma_x)
                - (y - yc)*(y - yc) / (2.0 * sigma_y * sigma_y);

            double envelope = std::exp(exponent);

            // Plane-wave phase factor e^{i (p_x x + p_y y)}
            std::complex<double> phase = std::exp(I * (p_x * x + p_y * y));

            std::complex<double> u = envelope * phase;

            U0(i,j) = u;
            norm += std::norm(u);   // |u|^2
        }
    }

    // Normalise so that Σ_{i,j} |u_ij|^2 = 1
    double inv_sqrt_norm = 1.0 / std::sqrt(norm);
    #pragma omp parallel for collapse(2)
    for (int j = 0; j < M; ++j) {
        for (int i = 0; i < M; ++i) {
            U0(i,j) *= inv_sqrt_norm;   // boundaries are zero, remain zero
        }
    }
}

//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------

// Initialise potential V(x,y) for a vertical wall with slits.
//
// Domain: [0,1] × [0,1], uniform grid with M points in each direction.
// V is real-valued here (arma::mat)
//
// Inputs:
//   V                : output matrix, size M×M
//   M                : number of grid points in x and y (including boundaries)
//   v0               : barrier height inside the walls
//   wall_thickness_x : thickness of the wall in x (physical units in [0,1])
//   wall_x_pos       : x-position of the wall centre (typically 0.5)
//   wall_sep_length  : length of the wall segment between neighbouring slits
//   slit_aperture    : slit opening in y (height of each slit, in [0,1])
//   n_slits          : number of slits (1, 2, 3, …)
//
// The slit pattern is constructed symmetric around y = 0.5.
void init_potential(arma::mat& V, int    M, double v0, double wall_thickness_x, double wall_x_pos, double wall_sep_length, double slit_aperture, int n_slits)
{
    V.set_size(M, M);
    V.zeros();                    // zero potential everywhere as baseline

    const double h = 1.0 / (M - 1);   // grid spacing
    const double y_center = 0.5;

    // --- 1. Identify x–indices belonging to the vertical wall ----------------

    std::vector<int> wall_i;
    wall_i.reserve(M);

    for (int i = 0; i < M; ++i) {
        double x = i * h;
        if (std::abs(x - wall_x_pos) <= 0.5 * wall_thickness_x) {
            wall_i.push_back(i);
        }
    }

    // If no index falls inside the requested thickness, nothing to do
    if (wall_i.empty()) {
        std::cout << "\n break \n";
        return;
    }

    // --- 2. Compute slit positions, symmetric around y = 0.5 -----------------

    // Total vertical span occupied by all slits and the separating wall pieces
    double total_span =
        n_slits * slit_aperture + (n_slits - 1) * wall_sep_length;

    double y_start = y_center - 0.5 * total_span;   // lower edge of the first slit

    // Precompute slit intervals [y_low, y_high] for each slit
    std::vector<std::pair<double,double>> slit_intervals;
    slit_intervals.reserve(n_slits);

    double current_y = y_start;
    for (int s = 0; s < n_slits; ++s) {
        double y_low  = current_y;
        double y_high = current_y + slit_aperture;
        slit_intervals.emplace_back(y_low, y_high);
        current_y += slit_aperture + wall_sep_length;
    }

    // --- 3. Fill wall region with v0, then carve out slits -------------------

    for (int idx_i : wall_i) {
        for (int j = 0; j < M; ++j) {
            double y = j * h;

            // Check if (x_i, y_j) lies inside any slit interval
            bool in_slit = false;
            for (const auto& interval : slit_intervals) {
                if (y >= interval.first && y <= interval.second) {
                    in_slit = true;
                    break;
                }
            }

            if (!in_slit) {
                // Inside the wall but outside the slits: set high potential
                V(idx_i, j) = v0;
            }
            // else: point lies in a slit -> keep V = 0
        }
    }
}

//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------

// Print an representation of the potential V.
// Cells with V > threshold are shown as '#', others as '.'.
// To avoid huge outputs, we downsample if M is large.
void print_potential_structure(const arma::mat& V, double threshold)
{
    int Mx = static_cast<int>(V.n_rows); // x index (i)
    int My = static_cast<int>(V.n_cols); // y index (j)

    // Choose a step to keep width reasonable in terminal
    int step_x = std::max(1, Mx / 80);  // max ~80 chars wide
    int step_y = std::max(1, My / 40);  // max ~40 lines high

    std::cout << "Potential structure ( '#' = barrier, '.' = free )\n";

    // Print y from top to bottom so it looks like a usual plot
    for (int j = My - 1; j >= 0; j -= step_y) {
        std::cout << "|";
        for (int i = 0; i < Mx; i += step_x) {
            if (V(i, j) > threshold) {
                std::cout << "#";
            } else {
                std::cout << ".";
            }
        }
        std::cout << "|\n";
    }
    std::cout << std::endl;
}

//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------

// Save a grayscale image of the potential V to a PGM file.
//
// High potential -> dark, low potential -> light.
// The image is oriented so that y=0 is at the bottom.
void save_potential_image(const arma::mat& V, const std::string& filename)
{
    int Mx = static_cast<int>(V.n_rows); // x index (i)
    int My = static_cast<int>(V.n_cols); // y index (j)

    // Find range of V for scaling
    double Vmin = V.min();
    double Vmax = V.max();
    if (Vmax == Vmin) {
        Vmax = Vmin + 1.0; // avoid division by zero
    }

    // Image: rows = y (vertical), cols = x (horizontal)
    // We'll store img(row=j_img, col=i) with j_img=0 at the TOP,
    // so we flip the y index to get a standard Cartesian view.
    arma::Mat<unsigned char> img(My, Mx);

    for (int j = 0; j < My; ++j) {
        for (int i = 0; i < Mx; ++i) {

            double v = V(i, j);
            double t = (v - Vmin) / (Vmax - Vmin); // normalise to [0,1]
            t = std::clamp(t, 0.0, 1.0);

            // High V -> dark (0), low V -> light (255)
            unsigned char gray = static_cast<unsigned char>(255.0 * (1.0 - t));

            int j_img = My - 1 - j;  // flip vertically: y=0 at bottom
            img(j_img, i) = gray;
        }
    }

    // Save as binary PGM image
    img.save(filename, arma::pgm_binary);
}

//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------

// Extract internal values of U (M×M) into a vector u of length (M-2)^2,
arma::cx_vec pack_internal_to_vec(const arma::cx_mat& U, int M)
{
    const int N_internal = M - 2;
    const int N          = N_internal * N_internal;

    arma::cx_vec u(N);

    for (int j = 1; j <= N_internal; ++j) {      // internal j: 1..M-2
        for (int i = 1; i <= N_internal; ++i) {  // internal i: 1..M-2
            int k = ij_to_k(i, j, M);           // 0-based index in [0, N-1]
            u(k) = U(i, j);                     // take value from U0
        }
    }

    return u;
}

//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------

// Fill an M×M matrix U from the internal vector u (length (M-2)^2).
// - U will include the boundary points and is initialised to zero.
// - Internal points (i,j = 1..M-2) are taken from u(k) using ij_to_k.
void unpack_vec_to_internal(arma::cx_mat& U, const arma::cx_vec& u, int M)
{
    const int N_internal = M - 2;
    const int N          = N_internal * N_internal;

    if ((int)u.n_elem != N) {
        throw std::runtime_error("unpack_vec_to_internal: wrong length of u");
    }

    // Allocate M×M and impose Dirichlet boundaries (all zeros)
    U.set_size(M, M);
    U.zeros();

    // Copy internal points from u back into U(i,j)
    for (int j = 1; j <= N_internal; ++j) {      // internal j: 1..M-2
        for (int i = 1; i <= N_internal; ++i) {  // internal i: 1..M-2
            int k = ij_to_k(i, j, M);           // same mapping as in pack_internal_to_vec
            U(i, j) = u(k);
        }
    }
}
