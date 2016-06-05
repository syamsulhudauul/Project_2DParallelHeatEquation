#include "matrix_template.h"
#include <vector>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#define self (*this)

template <class T>
Matrix<T>::Matrix(const std::string& name)
{
	self.name = name;
	self.row = 0;
	self.col = 0;
}

template <class T>
Matrix<T>::Matrix(
	const std::string& name,
	unsigned int row,
	unsigned int col) :
 matrix(row*col)
{
	self.row = row;
	self.col = col;
	self.name = name;
}

template <class T>
void Matrix<T>::assign(
	unsigned int row,
	unsigned int col)
{
	self.row = row;
	self.col = col;
	self.matrix.assign(row*col, 0);
}

template <class T>
T& Matrix<T>::operator()(
	unsigned int row,
	unsigned int col)
{
	//return matrix[row + col * self.row];
	return matrix[row*self.col + col];
	
}

template <class T>
T Matrix<T>::operator()(
	unsigned int row,
	unsigned int col) const
{
	//return matrix[row + col * self.row];
	return matrix[row*self.col + col];
}

template <class T>
std::vector<T> Matrix<T>::get_column(int col)
{
	std::vector<T> column(self.row, 0);
	
	for (int i = 0; i < self.row; i++) {
		column[i] = self(i, col);
	}
	
	return column;
}

template <class T>
std::vector<T> Matrix<T>::get_row(int row)
{
	std::vector<T> row_vector(0, self.col);
	
	for (int i = 0; i < self.col; i++) {
		row_vector[i] = self(row, i);
	}
	
	return row_vector;
}

template <class T>
void Matrix<T>::set_column(
	int col,
	std::vector<T> vector)
{	
	for (int i = 0; i < self.row; i++) {
		self(i, col) = vector[i];
	}
}

template <class T>
void Matrix<T>::print_matrix()
{
	for (int i = 0; i < self.row; i ++) {
		for (int j = 0; j < self.col; j++) {
			std::cout << self(i, j);
		}
		std::cout << std::endl;
	}
}

template class Matrix<double>;
//template class Matrix<dco_datatype>;