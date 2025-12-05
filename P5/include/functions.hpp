#ifndef __FUNCTIONS_HPP__
#define __FUNCTIONS_HPP__

#include <stdexcept>
#include <armadillo>
#include <vector>
#include <string>
#include <complex>
#include <stdexcept>
#include <iostream>
#include <algorithm>

// Function 1: Map internal grid indices (i,j) to a single vector index k.
int ij_to_k(int i, int j, int M);

// Function 2: Print the structure of a sparse matrix to screen.
void print_sp_matrix_structure(const arma::sp_cx_mat& A);

// Function 3: Build sparse matrices A and B for the 2D Crank–Nicolson scheme.
void construct_AB_matrices(const arma::cx_vec& a, const arma::cx_vec& b, std::complex<double> r, int M, arma::sp_cx_mat& A, arma::sp_cx_mat& B);

// Function 4: Build the diagonal vectors a(N) and b(N) from M, h, Δt and V.
void build_ab_vectors(int M, double dt, const arma::mat& V, arma::cx_vec& a, arma::cx_vec& b);

// Funtion 5: Perform the matrix multiplication B·u_n=b
void matrix_mult(const arma::sp_cx_mat& B, const arma::cx_vec& u_n, arma::cx_vec& b);

//Function 6: Perform one Crank–Nicolson time step
bool cn_step(const arma::sp_cx_mat& A, const arma::sp_cx_mat& B, const arma::cx_vec& u, arma::cx_vec& u_next);

//Function 7: Set up initial state u(x,y,t=0) as a Gaussian normalised wave packet
void init_gaussian_packet(arma::cx_mat& U0, int M, double xc, double yc, double sigma_x, double sigma_y, double p_x, double p_y);

//Function 8: Initialise potential V(x,y) for a vertical wall with slits.
void init_potential(arma::mat& V, int M, double v0 = 1.0e10, double wall_thickness_x = 0.02, double wall_x_pos = 0.5, double wall_sep_length = 0.05, double slit_aperture = 0.05, int n_slits = 2);

//Function 9: Print a representation of the potential V in the terminal.
void print_potential_structure(const arma::mat& V, double threshold);

//Function 10: Save a grayscale image of the potential V to a PGM file.
void save_potential_image(const arma::mat& V, const std::string& filename);

//Function 11: Extract internal values of U (M×M) into a vector u of length (M-2)^2
arma::cx_vec pack_internal_to_vec(const arma::cx_mat& U, int M);

//Function 12: Fill an M×M matrix U from the internal vector u (length (M-2)^2)
void unpack_vec_to_internal(arma::cx_mat& U, const arma::cx_vec& u, int M);

#endif