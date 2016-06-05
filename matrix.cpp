#include "matrix.h"
#include <vector>
#include <string>
#include <iostream>

#define self (*this)

Matrix::Matrix(const std::string& name)
{
	self.name = name;
	self.row = 0;
	self.col = 0;
}

Matrix::Matrix(
	const std::string& name,
	unsigned int row,
	unsigned int col) :
 matrix(row*col)
{
	self.row = row;
	self.col = col;
	self.name = name;
}

void Matrix::assign(
	unsigned int row,
	unsigned int col)
{
	self.row = row;
	self.col = col;
	self.matrix.assign(row*col, 0);
}

double& Matrix::operator()(
	unsigned int row,
	unsigned int col)
{
	//return matrix[row + col * self.row];
	return matrix[row*self.col + col];
	
}

double Matrix::operator()(
	unsigned int row,
	unsigned int col) const
{
	//return matrix[row + col * self.row];
	return matrix[row*self.col + col];
}

std::vector<double> Matrix::get_column(int col)
{
	std::vector<double> column(self.row, 0);
	
	for (int i = 0; i < self.row; i++) {
		column[i] = self(i, col);
	}
	
	return column;
}

std::vector<double> Matrix::get_row(int row)
{
	std::vector<double> row_vector(0, self.col);
	
	for (int i = 0; i < self.col; i++) {
		row_vector[i] = self(row, i);
	}
	
	return row_vector;
}

void Matrix::set_column(
	int col,
	std::vector<double> vector)
{	
	for (int i = 0; i < self.row; i++) {
		self(i, col) = vector[i];
	}
}

void Matrix::print_matrix()
{
	for (int i = 0; i < self.row; i ++) {
		for (int j = 0; j < self.col; j++) {
			fprintf(stdout, "%g	", self(i,j));
		}
		fprintf(stdout, "\n");
	}
}