#include "functions.hpp"

using namespace std;

// COMPILING AND LINKING FROM P5: g++ -std=c++17 -O2 -I include src/main.cpp src/utils/*.cpp -larmadillo -o src/main.exe -fopenmp

int main() {

    int option;

    // Show message in the terminal
    cout << "\n=== Double-slit setup ===\n"
          << " Choose an option (write a number)\n"  
          << " 1. Create A and B matrices and check. \n"
          << " 2. Initialise the potential. \n"
          << " 3. Run a simulation. \n";

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

        cout << "\nEnter the characteristics of the potential wall: \n";
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
        cout << "\nChoose M for the grid dimension: ";
        cin >> M;
        double dt;
        cout << "Choose dt for the timestep: ";
        cin >> dt;
        double T;
        std::cout << "Choose total simulation time T: ";
        std::cin  >> T;
        int Nt = static_cast<int>(T / dt);
        const double h = 1.0 / (M - 1);     // grid spacing in [0,1]

        int potential;
        cout << "\n Choose your potential (write a number)\n"
          << " 1. Predetermined\n"  
          << " 2. V = 0 for the entire grid\n"
          << " 3. Define your potential\n";

        cin >> potential;

        arma::mat V = arma::zeros(M, M);
        if (potential == 1){
            init_potential(V, M);
            cout << "\nThe characteristics of the potential are:\n";
            cout << "potential value in the wall: 70\n"; 
            cout << "wall thickness: 0.02\n";
            cout << "wall position in x: 0.5\n";
            cout << "space between slites: 0.05\n";
            cout << "slit aperture: 0.05\n";
            cout << "number of slits: 2\n";
        }
        if (potential == 2){
            cout << "\nYour don't have a wall.\n";
        }
        if (potential == 3){
            cout << "\n Enter the characteristics of the potential wall: \n";
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
        }

        arma::cx_vec a, b;
        build_ab_vectors(M, dt, V, a, b);

        // r for the Crank–Nicolson stencil 
        std::complex<double> I(0.0, 1.0);
        std::complex<double> r = I * dt / (2.0 * h * h);

        arma::sp_cx_mat A, B;
        construct_AB_matrices(a, b, r, M, A, B);

        cout << "\n Enter the characteristics of the Gaussian wave packet (u0): \n";
        double xc; double yc; double sigma_x; double sigma_y; double p_x; double p_y;
        cout << "centre of the initial wave packet x: ";
        cin >> xc; 
        cout << "centre of the initial wave packet y: ";
        cin >> yc; 
        cout << "initial width of the wave packet x: ";
        cin >> sigma_x; 
        cout << "initial width of the wave packet y: ";
        cin >> sigma_y; 
        cout << "wave packet momenta x: ";
        cin >> p_x;
        cout << "wave packet momenta y: ";
        cin >> p_y; 

        arma::cx_mat U0;
        init_gaussian_packet(U0, M, xc, yc, sigma_x, sigma_y, p_x, p_y);

        // pack U0 (matrix) into u (vector of internal points)
        arma::cx_vec u = pack_internal_to_vec(U0, M);

        // Allocate cube to store all time states
        arma::cx_cube U(M, M, Nt+1, arma::fill::zeros);
        U.slice(0) = U0;   // store initial wavefunction

        cout << "\n Now what? Choose an option\n"
          << " 1. Study the deviation of the total probability from 1.0.\n"  
          << " 2. Study time evolution of the 2D probability function.\n"
          << " 3. Analyze the real and imaginary parts of the wave function.\n"
          << " 4. Analyze the detections on a screen. \n"
          << " 5. Create an animations of your simulation.\n";

        cin >> option;

        if (option == 1){
        
            // Time-stepping loop
            arma::cx_vec u_next;

            // A file to store the probability deviation
            ofstream prob_file("data/prob_deviation.txt");

            for (int n = 0; n < Nt; ++n) {

                bool ok = cn_step(A, B, u, u_next);
                if (!ok) {
                    std::cerr << "Crank–Nicolson step failed at n = " << n << "\n";
                    break;
                }

                // ---- total probability and deviation from 1 ----

                double prob = 0.0;
                for (arma::uword k = 0; k < u_next.n_elem; ++k) {
                    prob += std::norm(u_next(k));   // |u_k|^2
                }
                
                //double prob = arma::accu( arma::square( arma::abs(u_next) ) );

                double deviation = std::abs(prob - 1.0);

                if (prob_file) {
                    prob_file << n << "  " << deviation << "\n";
                }

                // Convert u_next (internal vector) back to matrix form Unext_mat
                arma::cx_mat Unext_mat(M, M, arma::fill::zeros);
                unpack_vec_to_internal(Unext_mat, u_next, M);

                // Store as slice n+1
                U.slice(n + 1) = Unext_mat;

                // Prepare for next iteration: u <- u_next
                u = u_next;
            }
        
            prob_file.close();
            U.save("data/wavefunction_cube.bin");  
            cout << "\n Run the python script from scripts/prob_deviation.py to see your probability deviations results \n";

        }
    }   

}