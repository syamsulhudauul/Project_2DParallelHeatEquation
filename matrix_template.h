#ifndef MATRIX_H_
#define MATRIX_H_

#include <string>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

template <class T>
class Matrix {
public:
	Matrix(const std::string& name);
	Matrix(const std::string& name, unsigned int row, unsigned int col);
	
	T& operator()(unsigned int row, unsigned int col);
	T operator()(unsigned int row, unsigned int col) const;
	
	std::vector<T> matrix;
	unsigned int row;
	unsigned int col;
	std::string name;
	
	std::vector<T> get_column(int col);
	
	std::vector<T> get_row(int row);
	
	void assign(
	unsigned int row,
	unsigned int col);
	
	void set_column(
	int col,
	std::vector<T> vector);
	
	void print_matrix();
};

#endif /* MATRIX_H_ */