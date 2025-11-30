#include "functions.hpp"

using namespace std;

// COMPILING AND LINKING FROM P5: g++ -std=c++17 -O2 -I include src/main.cpp src/utils/*.cpp -larmadillo -o src/main.exe -fopenmp

int main() {

    int option;

    // Show message in the terminal
    cout << "\n=== ISING MODEL SIMULATION MENU ===\n"
          << " Choose an option\n"  
          << " 1. Create A and B matrices and check. \n"
          << " 2. \n"
          << " 3. \n"
          << " 4. \n"
          << " 5. \n";

    cin >> option;

    if (option == 1){

        int M;
        cout << "Choose M for the A and B matrices of size (M-2)^2 × (M-2)^2\n";
        cin >> M;
        double h = 0.001;
        double dt = 0.001;

        arma::mat V = arma::eye(M, M);
        // fill V(i,j) ...

        arma::cx_vec a, b;
        build_ab_vectors(M, h, dt, V, a, b);

        // r for the Crank–Nicolson stencil (same expression as in build_ab_vectors)
        std::complex<double> I(0.0, 1.0);
        std::complex<double> r = I * dt / (2.0 * h * h);

        arma::sp_cx_mat A, B;
        construct_AB_matrices(a, b, r, M, A, B);

        // e.g. to inspect the sparsity pattern:
        cout << "A matrix \n";     
        print_sp_matrix_structure(A);
        cout << "B matrix \n";
        print_sp_matrix_structure(B);
        
    }

    //if (option == 2){}
    //if (option == 3){}   
    //if (option == 4){}
    //if (option == 5){} 
}