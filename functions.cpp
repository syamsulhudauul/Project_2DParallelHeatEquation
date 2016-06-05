
#include <math.h>
#include <vector>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include "definitions.h"
#include "mpi.h"
#include "functions.h"
#include "matrix.h"
#include "matrix_vector_product.h"

// Function to initialize the temperature distribution in the 2D plate
void initialize_data(Matrix<double> &solution)
{
	for (int i = 0; i < grid_size; i++) {
		for (int j = 0; j < grid_size; j++) {
			if (i == 0 || i == grid_size-1 || j == 0 || j == grid_size-1) {
				solution(i, j) = 300;
			}
			
			else {
				solution(i, j) = 273;
			}
		}
	}
}

// Function to initialize the heat generation grid
void initialize_source_data(Matrix<double> &source_term)
{

	int x = source_term_x_coordinate;
	int y = source_term_y_coordinate;

	for (int i = 0; i < grid_size; i++) {
		for (int j = 0; j < grid_size; j++) {
			if (i >=x-5 && i <=x+5 && j >= y-5 && j <= y+5) {
				source_term(i, j) = (double) source_term_value;
			}
			
			else {
				source_term(i, j) = 0;
			}

		}
	}
}

// Function to initialize the thermal conductivity matrix
void initialize_thermal_conductivity(
	Matrix<double> &thermal_conductivity_matrix,
	double spatial_resolution
)
{
	for (int i = 0; i < grid_size; i++) {
		for (int j = 0; j < grid_size; j++) {
			thermal_conductivity_matrix(i, j) = thermal_conductivity;
		}
	}
	
}

// Function to initialize the intermediate thermal conductivity matrix
void initialize_intermediate_thermal_conductivity(
	Matrix<double> &intermediate_thermal_conductivity_matrix,
	double spatial_resolution
)
{
	
	std::srand(std::time(0));
	for (int i = 0; i < grid_size; i++) {
		for (int j = 0; j < grid_size; j++) {
			intermediate_thermal_conductivity_matrix(i, j) = thermal_conductivity + ((std::rand()%10))/10;
		}
	}
}

// Function to initialize the coefficient matrix
void populate_coefficients(
	Matrix<double> &alpha,
	Matrix<double> &thermal_conductivity_matrix,
	double delta_t,
	double spatial_resolution
)
{
	for (int i = 0; i < grid_size; i++) {
		for (int j = 0; j < grid_size; j++) {
			alpha(i, j) = (delta_t * thermal_conductivity_matrix(i, j)) / (spatial_resolution*spatial_resolution);
		}
	}
}

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
)
{
	populate_coefficients(alpha, thermal_conductivity_matrix, delta_t, spatial_resolution);
	
	// Begin doing Jacobi iterations
	for (int t = 0; t < time_steps; t++) {
		// Initialize solution to 0
		// Ignore boundaries
		for (int i = start; i <= end; i++) {
			for (int j = 1; j < solution_new.col-1; j++) {
				solution_new(i, j) = 0.0;
			}
		}
		
		// Inner iterations
		for (int m = 0; m < 1000; m++) {
			// Send border rows to ghost rows for neighboring chunks
			if (left_neighbor != 0) {
				MPI_Send(&solution(offset, 0), grid_size, MPI_DOUBLE, left_neighbor, RIGHT_TAG, MPI_COMM_WORLD);
				MPI_Recv(&solution(offset-1, 0), grid_size, MPI_DOUBLE, left_neighbor, LEFT_TAG, MPI_COMM_WORLD, &status);
			}
			
			if (right_neighbor != 0) {
				MPI_Send(&solution(offset+rows_per_worker-1, 0), grid_size, MPI_DOUBLE, right_neighbor, LEFT_TAG, MPI_COMM_WORLD);
				MPI_Recv(&solution(offset+rows_per_worker, 0), grid_size, MPI_DOUBLE, right_neighbor, RIGHT_TAG, MPI_COMM_WORLD, &status);
			}
			
			// Update right hand side and solution_new
			// Ignore boundaries
			for (int i = start; i <= end; i++) {
				for (int j = 1; j < solution_new.col-1; j++) {
					right_hand_side[j] = solution_old(i, j) / (1.0 + 4 * alpha(i, j)) + (source_term(i, j) * delta_t / (density * specific_heat));
					
					solution_new(i, j) = right_hand_side[j] + (alpha(i, j) / (1.0 + 4 * alpha(i, j))) * (solution(i-1, j) + solution(i+1, j) + solution(i, j+1) + solution(i, j-1));
				}
			}
			
			// Update solution
			// Ignore boundaries
			for (int i = start; i <= end; i++) {
				for (int j = 1; j < solution.col-1; j++) {
					solution(i, j) = solution_new(i, j);
				}
			}
		}
		
		// Update solution_old
		solution_old = solution;
	}
		
}
