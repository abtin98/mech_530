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


using Matrix = std::vector<std::vector<double>>;

namespace MatrixOperations
{

void output_matrix (Matrix matrix, std::string name, std::string units);

void output_vector (std::vector<double> vector, std::string name, std::string units);

std::vector<double> vector_mult (Matrix matrix, std::vector<double> vector);

Matrix matrix_mult (Matrix matrix_left, Matrix matrix_right);

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


#endif /* UTILITIES_H_ */
