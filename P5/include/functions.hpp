#ifndef __FUNCTIONS_HPP__
#define __FUNCTIONS_HPP__

#include <stdexcept>
#include <armadillo>
#include <vector>
#include <string>
#include <complex>

// Function 1: Map internal grid indices (i,j) to a single vector index k.
int ij_to_k(int i, int j, int M);

// Function 2: Print the structure of a sparse matrix to screen.
void print_sp_matrix_structure(const arma::sp_cx_mat& A);

// Function 3: Build sparse matrices A and B for the 2D Crank–Nicolson scheme.
void construct_AB_matrices(const arma::cx_vec& a, const arma::cx_vec& b, std::complex<double> r, int M, arma::sp_cx_mat& A, arma::sp_cx_mat& B);

// Function 4: Build the diagonal vectors a(N) and b(N) from M, h, Δt and V.
void build_ab_vectors(int M, double h, double dt, const arma::mat& V, arma::cx_vec& a, arma::cx_vec& b);
#endif