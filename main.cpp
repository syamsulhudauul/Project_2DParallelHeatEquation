#include <math.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <ctime>
#include <algorithm>
//#include <omp.h>
#include "definitions.h"
#include "functions.h"
#include <mpi.h>
#include "matrix_template.h"
#include "matrix_vector_product.h"



int main(int argc, char *argv[])
{
	// <double> Solution terms
	auto solution_old = Matrix<double>("solution_old", grid_size, grid_size);
	auto solution = Matrix<double>("solution", grid_size, grid_size);
	auto solution_new = Matrix<double>("solution_new", grid_size, grid_size);
	
	
	// Calculation terms
	std::vector<double> right_hand_side(grid_size);
	double delta_t = (double)TIME / time_steps;
	double spatial_resolution = 1.0 / grid_size;
	
	// Create and populate the thermal conductivity matrix
	auto thermal_conductivity_matrix = Matrix<double>("thermal_conductivity_matrix", grid_size, grid_size);
	initialize_thermal_conductivity(thermal_conductivity_matrix, spatial_resolution);
		
	// Declare and popualate right hand side coefficient matrix
	auto alpha = Matrix<double>("alpha", grid_size, grid_size);
	
	// Declare and populate the heat generation matrix
	auto source_term = Matrix<double>("source_term", grid_size, grid_size);
	initialize_source_data(source_term);
	
	// Declare MPI variables
	int myrank;								// This process' ID
	int number_of_processes;				// Total number of processes
	int number_of_workers;					// Number of worker processes
	int average_rows_per_worker;			// grid_size / number_of_workers
	int extra_rows;							// Extra rows if grid_size % number_of_workers != 0
	int offset;								// Offset for sending data to workers
	int rows_per_worker;					// Number of rows sent to each worker
	int left_neighbor;						// This worker's left neighbor's ID
	int right_neighbor;						// This worker's right neighbor's ID
	int destination;						// Tag for MPI_Send destination
	int source;								// Tag for MPI_Receive source
	int tag;								// Variable for holding tag
	int start;								// Indexing variable for each process' chunk of data
	int end;								// Indexing variable for each process' chunk of data
	MPI_Group original_group, new_group;	// Group handle for Allreduce excluding master
	MPI_Comm new_communicator;				// Communicator handle for new group
	MPI_Status status;						// MPI Status
	
	//double start_time = omp_get_wtime();
	
	// Start MPI, find out number of processes and find out process rank.
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &number_of_processes);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	number_of_workers = number_of_processes - 1;
	int ranks[number_of_workers];
	
	// Put worker ranks into array
		for (int i = 1; i <= number_of_workers; i++) {
			ranks[i-1] = i;
		}
		
	// Set up subgroup for performing allreduce on only the workers for the convergence
	MPI_Comm_group(MPI_COMM_WORLD, &original_group);
	MPI_Group_incl(original_group, number_of_workers, ranks, &new_group);
	MPI_Comm_create(MPI_COMM_WORLD, new_group, &new_communicator);
	
	/****************************** MASTER CODE ******************************/
	if (myrank == 0) {
		// Initialize solution
		initialize_data(solution_old);
		//initialize_data_dco(solution_temporary_dco);
		
		// Find out how many rows to send to each worker, accounting for
		// an uneven split
		average_rows_per_worker = grid_size / number_of_workers;
		extra_rows = grid_size % number_of_workers;
		offset = 0;
		
		// For each worker
		for (int i = 1; i <= number_of_workers; i++) {
			// Find out how many rows it gets
			if (i <= extra_rows) {
				rows_per_worker = average_rows_per_worker + 1;
			}
			else {
				rows_per_worker = average_rows_per_worker;
			}
			
			// Find its neighbors
			if (i == 1) {
				left_neighbor = 0;
			}
			else {
				left_neighbor = i - 1;
			}
			
			if (i == number_of_workers) {
				right_neighbor = 0;
			}
			else {
				right_neighbor = i + 1;
			}
			
			// Send initial data
			destination = i;
			tag = INITIALIZE;
			
			MPI_Send(&offset, 1, MPI_INT, destination, tag, MPI_COMM_WORLD);
			MPI_Send(&rows_per_worker, 1, MPI_INT, destination, tag, MPI_COMM_WORLD);
			MPI_Send(&left_neighbor, 1, MPI_INT, destination, tag, MPI_COMM_WORLD);
			MPI_Send(&right_neighbor, 1, MPI_INT, destination, tag, MPI_COMM_WORLD);
			MPI_Send(&solution_old(offset, 0), rows_per_worker*grid_size, MPI_DOUBLE, destination, tag, MPI_COMM_WORLD);
			//MPI_Send(&solution_temporary_dco(offset, 0), rows_per_worker*grid_size, MPI_DOUBLE, destination, tag, MPI_COMM_WORLD);
			
			// Increment offset for next worker
			offset = offset + rows_per_worker;
		}
		
		// Get results from workers
		for (int i = 1; i <= number_of_workers; i++) {
			source = i;
			tag = FINALIZE;
			
			MPI_Recv(&offset, 1, MPI_INT, source, tag, MPI_COMM_WORLD, &status);
			MPI_Recv(&rows_per_worker, 1, MPI_INT, source, tag, MPI_COMM_WORLD, &status);
			MPI_Recv(&solution_old(offset, 0), rows_per_worker*grid_size, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
		}
		
		//double runtime = omp_get_wtime() - start_time;
		//std::cout << "The simulation ran for: " << runtime << " seconds." << std::endl;
		
		// Print out the solution (for debugging; formatting is not properly handled)
		std::ofstream T_end ("final_temperatures.csv");
		for (int i = 0; i < solution_old.row; i++) {
			for (int j = 0; j < solution_old.col; j++) {
				T_end << solution_old(i, j) << ",";
				//std::cout << solution_old(i, j) << "\t";
			}
			T_end << std::endl;
			//std::cout << std::endl;
		}
		T_end.close();
	}
	
	/****************************** WORKER CODE ******************************/
	if (myrank != 0) {
		// Receive initial data from master
		source = 0;
		tag = INITIALIZE;
		
		MPI_Recv(&offset, 1, MPI_INT, source, tag, MPI_COMM_WORLD, &status);
		MPI_Recv(&rows_per_worker, 1, MPI_INT, source, tag, MPI_COMM_WORLD, &status);
		MPI_Recv(&left_neighbor, 1, MPI_INT, source, tag, MPI_COMM_WORLD, &status);
		MPI_Recv(&right_neighbor, 1, MPI_INT, source, tag, MPI_COMM_WORLD, &status);
		MPI_Recv(&solution_old(offset, 0), rows_per_worker*grid_size, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
		//MPI_Recv(&solution_temporary_dco(offset, 0), rows_per_worker*grid_size, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
		
		// Make sure solution gets boundary data
		for (int i = offset; i < offset+rows_per_worker; i++) {
			for (int j = 0; j < solution_old.col; j++) {
				solution(i, j) = solution_old(i, j);
			}
		}
		
		// Determine border elements
		start = offset;
		end = offset + rows_per_worker - 1;
		
		if (offset == 0) {
			start = 1;
		}
		
		if ((offset + rows_per_worker) == grid_size) {
			end = offset + rows_per_worker - 2;
		}
		
		// Perform Jacobi iterations
		jacobi_iterations(thermal_conductivity_matrix, solution_old, solution, solution_new, source_term, alpha, right_hand_side, start, end, offset, left_neighbor, right_neighbor, rows_per_worker, spatial_resolution, delta_t, myrank, new_communicator, status);
		
		// Send data back to master processor
		destination = 0;
		tag = FINALIZE;
		
		MPI_Send(&offset, 1, MPI_INT, destination, tag, MPI_COMM_WORLD);
		MPI_Send(&rows_per_worker, 1, MPI_INT, destination, tag, MPI_COMM_WORLD);
		MPI_Send( &solution_old(offset, 0), rows_per_worker*grid_size, MPI_DOUBLE, destination, tag, MPI_COMM_WORLD);
	}
	
	// **Gracefully** exit MPI
	MPI_Finalize();
}