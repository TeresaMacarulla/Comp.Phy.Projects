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
          << " 3. \n"
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
        cout << "grid dimension, potential value in the wall, wall thickness, wall position in x, space between slites, slit aperture, number of slits\n";
        cin >> M >> v0 >> wall_thickness_x >> wall_x_pos >> wall_sep_length >> slit_aperture >> n_slits;

        cout << "\n The characteristics of your wall: \n"
          << " grid dimension " << M << "\n"
          << " v0 " << v0 << "\n"
          << " wall thickness " << wall_thickness_x << "\n"
          << " wall position in x " << wall_x_pos << "\n"
          << " space between slites " << wall_sep_length << "\n"
          << " slit aperture " << slit_aperture << "\n"
          << " number of slits " << n_slits << "\n";
          
        init_potential(V, M, v0, wall_thickness_x, wall_x_pos, wall_sep_length, slit_aperture, n_slits);
        //V.print("Matrix V:");
    }
    //if (option == 3){}   
    //if (option == 4){}
    //if (option == 5){} 
}