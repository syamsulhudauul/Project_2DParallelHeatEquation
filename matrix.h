#ifndef MATRIX_H_
#define MATRIX_H_

#include <string>
#include <vector>
#include <iostream>

class Matrix {
public:
	Matrix(const std::string& name);
	Matrix(const std::string& name, unsigned int row, unsigned int col);
	
	double& operator()(unsigned int row, unsigned int col);
	double operator()(unsigned int row, unsigned int col) const;
	
	std::vector<double> matrix;
	unsigned int row;
	unsigned int col;
	std::string name;
	
	std::vector<double> get_column(int col);
	
	std::vector<double> get_row(int row);
	
	void assign(
	unsigned int row,
	unsigned int col);
	
	void set_column(
	int col,
	std::vector<double> vector);
	
	void print_matrix();
};

#endif /* MATRIX_H_ */