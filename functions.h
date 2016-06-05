#ifndef functions_h
#define functions_h

#include <math.h>
#include <vector>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include "definitions.h"
#include "mpi.h"
#include "matrix_template.h"
#include "matrix_vector_product.h"


void initialize_data(Matrix<double> &solution);


void initialize_source_data(Matrix<double> &source_term);

void initialize_thermal_conductivity(
	Matrix<double> &thermal_conductivity_matrix,
	double spatial_resolution
);

void initialize_intermediate_thermal_conductivity(
	Matrix<double> &intermediate_thermal_conductivity_matrix,
	double spatial_resolution
);

void populate_coefficients(
	Matrix<double> &alpha,
	Matrix<double> &thermal_conductivity_matrix,
	double delta_t,
	double spatial_resolution
);


void jacobi_iterations(
	Matrix<double> &thermal_conductivity_matrix,
	Matrix<double> &solution_old,
	Matrix<double> &solution,
	Matrix<double> &solution_new,
	Matrix<double> &source_term,
	Matrix<double> &alpha,
	std::vector<double> right_hand_side,
	int start,
	int end,
	int offset,
	int left_neighbor,
	int right_neighbor,
	int rows_per_worker,
	double spatial_resolution,
	double delta_t,
	int myrank,
	MPI_Comm new_communicator,
	MPI_Status status
);

#endif /* functions_h */