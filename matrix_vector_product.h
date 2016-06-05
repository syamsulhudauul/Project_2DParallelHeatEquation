#ifndef matrix_vector_product_h
#define matrix_vector_product_h

#include <iostream>
#include <math.h>
#include <vector>
#include "matrix_template.h"


std::vector<double> matrix_vector_product(int M, std::vector<unsigned int> const &I_csr, std::vector<unsigned int> const &J_csr, std::vector<double> const &x, std::vector<double> const &val_csr, int symmetric);

double vector_magnitude(std::vector<double> const v);

double dot_product(std::vector<double> u, std::vector<double> v);

std::vector<double> scalar_product(double h, std::vector<double> v);

std::vector<double> scalar_division(double h, std::vector<double> v);

void check_orthogonality(Matrix<double> V, FILE *file_name);

std::vector<double> vector_subtraction(std::vector<double> u, std::vector<double> v);

std::vector<double> vector_sum(
	std::vector<double> u,
	std::vector<double> v);

std::vector<double> backward_substitution(
	std::vector<double> x,
	Matrix<double> U);

std::vector<double> dense_matrix_vector_product(
	std::vector<double> x,
	Matrix<double> A);

double infinity_norm(Matrix<double> u);

double infinity_norm_master_worker(
	Matrix<double> u,
	int start,
	int end
);

Matrix<double> matrix_subtraction(
	Matrix<double> const &u,
	Matrix<double> const &v);

Matrix<double> matrix_addition(
	Matrix<double> const &u,
	Matrix<double> const &v);

Matrix<double> matrix_subtraction_absolute(
	Matrix<double> &u,
	Matrix<double> &v);

Matrix<double> opposite_matrix(Matrix<double> u);

#endif /* matrix_vector_product_h */