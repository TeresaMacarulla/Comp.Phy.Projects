#include "functions.hpp"

using namespace std;

// COMPILING AND LINKING FROM P5: g++ -std=c++17 -O2 -I include src/main.cpp src/utils/*.cpp -larmadillo -o src/main.exe -fopenmp

int main() {

    int option;

    // Show message in the terminal
    cout << "\n=== Double-slit setup ===\n"
          << " Choose an option\n"  
          << " 1. Create A and B matrices and check. \n"
          << " 2. Initialise the potential. \n"
          << " 3. Run a simulation. \n"
          << " 4. \n"
          << " 5. \n";

    cin >> option;

    if (option == 1){

        int M;
        cout << "Choose M for the A and B matrices of size (M-2)^2 × (M-2)^2\n";
        cin >> M;
        const double dt = 0.001;
        const double h = 1.0 / (M - 1);     // grid spacing in [0,1]

        arma::mat V = arma::eye(M, M);

        arma::cx_vec a, b;
        build_ab_vectors(M, dt, V, a, b);

        // r for the Crank–Nicolson stencil 
        std::complex<double> I(0.0, 1.0);
        std::complex<double> r = I * dt / (2.0 * h * h);

        arma::sp_cx_mat A, B;
        construct_AB_matrices(a, b, r, M, A, B);

        // e.g. to inspect the sparsity pattern:
        cout << "A matrix structure \n";     
        print_sp_matrix_structure(A);
        cout << "B matrix structure\n";
        print_sp_matrix_structure(B);     
    }

    if (option == 2){

        cout << "Enter the characteristics of the potential wall: \n";
        int M; double v0; double wall_thickness_x; double wall_x_pos; double wall_sep_length; double slit_aperture; int n_slits;
        arma::mat V = arma::zeros(M, M);
        cout << "grid dimension: ";
        cin >> M; 
        cout << "potential value in the wall: ";
        cin >> v0; 
        cout << "wall thickness: ";
        cin >> wall_thickness_x; 
        cout << "wall position in x: ";
        cin >> wall_x_pos; 
        cout << "space between slites: ";
        cin >> wall_sep_length;
        cout << "slit aperture: ";
        cin >> slit_aperture; 
        cout << "number of slits: ";
        cin >> n_slits; 
          
        init_potential(V, M, v0, wall_thickness_x, wall_x_pos, wall_sep_length, slit_aperture, n_slits);
        // --- Visualise potential in terminal (optional) ---
        print_potential_structure(V, 0.5 * v0);
        save_potential_image(V, "plots/potential.pgm");
        cout << "A finer representation of your potential can be found in plots/potential.pgm \n"; 
    }
    
    if (option == 3){

        int M;
        cout << "Choose M for the grid dimension\n";
        cin >> M;
        double dt;
        cout << "Choose dt for the timestep\n";
        cin >> dt;
        const double h = 1.0 / (M - 1);     // grid spacing in [0,1]

        arma::mat V = arma::zeros(M, M);
        init_potential(V, M);

        arma::cx_vec a, b;
        build_ab_vectors(M, dt, V, a, b);

        // r for the Crank–Nicolson stencil 
        std::complex<double> I(0.0, 1.0);
        std::complex<double> r = I * dt / (2.0 * h * h);

        arma::sp_cx_mat A, B;
        construct_AB_matrices(a, b, r, M, A, B);

        cout << "Enter the characteristics of the Gaussian wave packet (u0): \n";
        double xc; double yc; double sigma_x; double sigma_y; double p_x; double p_y;
        cout << "centre of the initial wave packet x: ";
        cin >> xc; 
        cout << "centre of the initial wave packet y: ";
        cin >> yc; 
        cout << "initial widths of the wave packet x: ";
        cin >> sigma_x; 
        cout << "initial widths of the wave packet y: ";
        cin >> sigma_y; 
        cout << "wave packet momenta x: ";
        cin >> p_x;
        cout << "wave packet momenta y: ";
        cin >> p_y; 

        arma::cx_mat U0;
        init_gaussian_packet(U0, M, xc, yc, sigma_x, sigma_y, p_x, p_y);

        // pack U0 (matrix) into u (vector of internal points)
        arma::cx_vec u = pack_internal_to_vec(U0, M);

        // now you can do one CN step
        arma::cx_vec u_next;
        bool ok = cn_step(A, B, u, u_next);
    }   
    //if (option == 4){}
    //if (option == 5){} 
}