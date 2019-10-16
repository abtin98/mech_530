/*
 * utilities.h
 *
 *  Created on: Sep 23, 2019
 *      Author: abtinameri
 */

#ifndef UTILITIES_H_
#define UTILITIES_H_

#include <iostream>
#include <array>
#include <vector>
#include <assert.h>
#include "properties.h"


using Matrix = std::vector<std::vector<double>>;

namespace MatrixOperations
{
void output_matrix (Matrix matrix, std::string name, std::string units);

void output_vector (std::vector<double> vector, std::string name, std::string units);

std::vector<double> vector_mult (Matrix matrix, std::vector<double> vector);

Matrix matrix_mult (Matrix matrix_left, Matrix matrix_right);

Matrix calculate_adjoint (Matrix matrix);

Matrix obtain_submatrix (Matrix matrix, int row, int column);

double calculate_determinant (Matrix matrix);

Matrix matrix_inverse (Matrix matrix);

Matrix calculate_co_factor_matrix (Matrix matrix);

std::vector<double> solve_linear_system (Matrix A, std::vector<double> b);

std::vector<double> v_add (std::vector<double> v1, std::vector<double> v2);
}

namespace StressTransformation
{
Matrix calculate_stress_transform_off_to_on(double angle);
Matrix calculate_stress_transform_on_to_off(double angle);
}

namespace StrainTransformation
{
Matrix calculate_strain_transform_off_to_on(double angle);
Matrix calculate_strain_transform_on_to_off(double angle);
}

namespace ModulusTransform
{
Matrix calculate_modulus_transform_on_to_off(double angle, Matrix Q);
}

namespace ComplianceTransform
{
Matrix calculate_compliance_transform_on_to_off(double angle, Matrix S);
}

namespace OverallModulus
{
Matrix calculate_overall_in_plane_modulus(std::vector<AllParameters> parameters_vector, Matrix Q);
Matrix calculate_overall_in_plane_compliance(std::vector<AllParameters> parameters_vector, Matrix Q);
}

namespace FlexuralModulus
{
Matrix calculate_flexural_modulus (std::vector<AllParameters> parameters_vector, Matrix Q, double z_c);
}


#endif /* UTILITIES_H_ */
