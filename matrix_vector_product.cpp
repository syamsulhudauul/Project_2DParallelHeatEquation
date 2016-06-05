#include <iostream>
#include <math.h>
#include <vector>
#include "matrix_vector_product.h"
#include "matrix_template.h"


std::vector<double> matrix_vector_product(
	int M,
	std::vector<unsigned int> const &I_csr,
	std::vector<unsigned int> const &J_csr,
	std::vector<double> const &x,
	std::vector<double> const &val_csr,
	int symmetric)
{
	int i1, i2;
	std::vector<double> y(M, 0);
	
	for (int i = 0; i < M; i++) {
		i1 = I_csr[i];
		i2 = I_csr[i + 1];
		for (int j = i1; j < i2; j++) {
			y[i] += val_csr[j] * x[J_csr[j]];
			
			// If the matrix is symmetric
			if (i != J_csr[j] && symmetric == 1) {
				y[J_csr[j]] += val_csr[j] * x[i];
			}
		}
	}
	
	return y;
}

std::vector<double> dense_matrix_vector_product(
	std::vector<double> x,
	Matrix<double> A)
{
	std::vector<double> y(A.row, 0);
	

	for (int i = 0; i < A.row; i++) {
		for (int j = 0; j < A.col; j ++) {
			y[i] += A(i, j) * x[j];
		}
	}
	
	return y;
}

double vector_magnitude(std::vector<double> const v)
{
	double sum = 0;
	

	for (int i = 0; i < v.size(); i++) {
		sum += v[i]*v[i];
	}
	
	return sqrt(sum);
}

double dot_product(
	std::vector<double> u,
	std::vector<double> v)
{
	double sum = 0;
	

	for (int i = 0; i < u.size(); i++) {
		sum += u[i] * v[i];
	}
	
	return sum;
}

std::vector<double> scalar_product(
	double h,
	std::vector<double> v)
{
	std::vector<double> result(v.size(), 0);
	

	for (int i = 0; i < v.size(); i++) {
		result[i] = h * v[i];
	}
	
	return result;
}

std::vector<double> scalar_division(
	double h,
	std::vector<double> v)
{
	std::vector<double> result(v.size(), 0);
	

	for (int i = 0; i < v.size(); i++) {
		result[i] = v[i] / h;
	}
	
	return result;
}

void check_orthogonality(
	Matrix<double> V,
	FILE *file_name)
{
	double check;
	double tol = 1e-5;
	
	int passed = 1;
	int i = 0;
	//for (int i = 0; i < V.col; i++) {
		for (int j = 0; j < V.col; j ++) {
			
			check = fabs(dot_product(V.get_column(i), V.get_column(j)));
			fprintf(file_name, "%d	%lg\n", j, check);
			//check = dot_product(V.get_column(i), V.get_column(j));
			/*
			if (i == j) {
				if (fabs(check - 1) < tol) {
					//fprintf(stdout, "%d	%d	%g	PASSED\n", i, j, check);
				}
				else {
					passed = 0;
					fprintf(stdout, "%d	%d	%g	f\n", i, j, check);
				}
			}
			else {
				if (fabs(check) < tol) {
					//fprintf(stdout, "%d	%d	%g	PASSED\n", i, j, check);
				}
				else {
					passed = 0;
					fprintf(stdout, "%d	%d	%g	f\n", i, j, check);
				}
			}*/
		}
	//}
	fprintf(stdout, "%d\n", passed);
}

std::vector<double> vector_subtraction(
	std::vector<double> u,
	std::vector<double> v)
{
	std::vector<double> result(u.size(), 0);
	

	for (int i = 0; i < u.size(); i++) {
		result[i] = u[i] - v[i];
	}
	
	return result;
}

std::vector<double> vector_sum(
	std::vector<double> u,
	std::vector<double> v)
{
	std::vector<double> result(u.size(), 0);
	

	for (int i = 0; i < u.size(); i++) {
		result[i] = u[i] + v[i];
	}
	
	return result;
}

std::vector<double> backward_substitution(
	std::vector<double> x,
	Matrix<double> U)
{
	std::vector<double> y(U.col, 0);
	
	for (int i = y.size() - 1; i >= 0; i--) {
		y[i] = x[i];
		
		for (int j = i + 1; j < y.size(); j++) {
			y[i] = y[i] - U(i, j) * y[j];
		}
		
		y[i] = y[i] / U(i, i);
	}
	
	return y;
}

double infinity_norm(Matrix<double> u)
{
	std::vector<double> row_sum(u.row, 0);
	double ret = 0;
	
	for (int i = 0; i < u.row; i++) {
		double sum = 0;

		for (int j = 0; j < u.col; j++) {
			sum += fabs(u(i, j));
		}
		
		if (i == 0)
			ret = sum;
		else
			ret = std::max(ret, sum);
	}
	return ret;
}

double infinity_norm_master_worker(
	Matrix<double> u,
	int start,
	int end
)
{
	std::vector<double> row_sum(u.row, 0);
	double ret = 0;
	
	for (int i = start; i <= end; i++) {
		double sum = 0;

		for (int j = 1; j < u.col-1; j++) {
			sum += fabs(u(i, j));
		}
		
		if (i == 0)
			ret = sum;
		else
			ret = std::max(ret, sum);
	}
	return ret;
}

Matrix<double> matrix_subtraction(
	Matrix<double> const &u,
	Matrix<double> const &v)
{
	auto result = Matrix<double>("result");
	result.assign(u.row, u.col);
	

	for (int j = 0; j < u.row; j++) {
		for (int i = 0; i < u.col; i++) {
			result(i, j) = u(i, j) - v(i, j);
		}
	}
	
	return result;
}

Matrix<double> matrix_addition(
	Matrix<double> const &u,
	Matrix<double> const &v)
{
	auto result = Matrix<double>("result");
	result.assign(u.row, u.col);
	

	for (int j = 0; j < u.row; j++) {
		for (int i = 0; i < u.col; i++) {
			result(i, j) = u(i, j) + v(i, j);
		}
	}
	
	return result;
}

Matrix<double> matrix_subtraction_absolute(
	Matrix<double> &u,
	Matrix<double> &v)
{
	auto result = Matrix<double>("result");
	result.assign(u.row, u.col);
	

	for (int j = 0; j < u.row; j++) {
		for (int i = 0; i < u.col; i++) {
			result(i, j) = fabs(u(i, j) - v(i, j));
		}
	}
	
	return result;
}

Matrix<double> opposite_matrix(Matrix<double> u)
{
	for (int i = 0; i < u.row; i++) {
		for (int j = 0; j < u.col; j++) {
			u(i, j) = -u(i, j);
		}
	}
	
	return u;
}