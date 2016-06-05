#ifndef DEFINITIONS_H
#define DEFINITIONS_H

/*******************************
 *  CHANGE THESE VARIABLES TO  *
 * MODIFY THE SIMULATION INPUT *
 *******************************/

// Solver Definitions
#define grid_size					101		// Number of nodes per side in the square grid
#define TIME						100		// How much time should be simulated (seconds)
#define CONVERGENCE_CRITERIA		10e-8	// Convergence criteria for Jacobi iterations
#define time_steps					100	    // How many time steps to reach TIME
#define source_term_x_coordinate	50		// The x-coordinate of the heat generation node
#define source_term_y_coordinate	50		// The y-coordinate of the heat generation node
#define source_term_value			1000000	// The value of heat generation (Joules)


/***********************
 * Material Properties *
 ***********************/

// Currently using SAE 1010 steel
#define thermal_conductivity		59		// Thermal conductivity of the material (W/m*K)
#define specific_heat				434		// Specific heat of the material (J/kg*K)
#define density						7832	// Density of the material of the plate (kg/m3)

/*******************
 * MPI DEFINITIONS *
 *  DO NOT MODIFY  *
 *******************/
#define INITIALIZE	0		// Message tag for initialization data
#define LEFT_TAG	1		// Message tag for sending to left neighbor
#define RIGHT_TAG	2		// Message tag for sending to right neighbor
#define FINALIZE	3		// Message tag for sending back to master

#endif // DEFINITIONS_H