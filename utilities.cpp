/*
 * utilities.cpp
 *
 *  Created on: Sep 23, 2019
 *      Author: abtinameri
 */


#include "utilities.h"
#include <math.h>

#ifndef PI
#define PI 3.14159265358979323846
#endif


using Matrix = std::vector<std::vector<double>>;

namespace MatrixOperations
{
void output_matrix(Matrix matrix, std::string name, std::string units)
{
    std::cout << "The " << name << " matrix is (" << units << "):" <<  std::endl;
    for (int i = 0; i < matrix.size() ; ++i)
    {
        for (int j = 0; j < matrix[0].size(); ++j)
        {
            std::cout << matrix[i][j] << " | ";
        }
        std::cout << std::endl;
    }
    std::cout << "--------------------------------------------" << std::endl;
}

void output_vector (std::vector<double> vector, std::string name, std::string units)
{
    std::cout << "The " << name << " vector is (" << units << "):" <<  std::endl;
    for (int i = 0; i < vector.size() ; ++i)
    {
        std::cout << vector[i] << std::endl;
    }
    std::cout << "--------------------------------------------" << std::endl;
}

std::vector<double> vector_mult (Matrix matrix, std::vector<double> vector)
{
    assert (matrix[0].size() == vector.size());
    std::vector<double> result(vector.size());

    for (int i = 0; i < matrix.size(); ++i)
    {
        result[i] = 0;
        for (int j = 0 ; j < vector.size(); ++j)
        {
            result[i] += matrix[i][j] * vector[j];
        }
    }

    return result;
}

Matrix matrix_mult (Matrix matrix_left, Matrix matrix_right)
{
    Matrix result;

    assert(matrix_left[0].size() == matrix_right.size());

    result.resize(matrix_left.size(),std::vector<double>(matrix_right[0].size()));

    for (int i = 0; i < matrix_left.size(); ++i)
    {
        for (int j = 0; j < matrix_right[0].size(); ++ j)
        {
            result[i][j] = 0;
            for (int k = 0; k < matrix_left[0].size(); ++k)
            {
                result[i][j] += matrix_left[i][k] * matrix_right[k][j];
            }
        }
    }

    return result;
}

//Matrix matrix_inverse (Matrix matrix)
//{
//    const unsigned int dim = matrix[0].size();
//    assert (dim == matrix.size());
//
//    Matrix inv_matrix;
//    inv_matrix.resize(dim,std::vector<double>(dim));
//
//
//
//    return inv_matrix;
//}

//std::vector<double> solve_linear_system (Matrix A, std::vector<double> b)
//{
//    std::vector<double> x;
//    assert (A.size() == A[0].size());
//    assert (A[0].size() == b.size());
//
//    x.resize (A.size());
//
//
//
//    return x;
//}

Matrix obtain_submatrix (Matrix matrix, int row, int column)
{
    Matrix submatrix;

    submatrix.resize( matrix.size() - 1,std::vector<double>( matrix[0].size() - 1 )  );

    int submatrix_row = 0;



    for (int i = 0; i < matrix.size(); ++i)
    {
        if (i == row) continue;
        int submatrix_col = 0;
        for (int j = 0 ; j < matrix[0].size(); ++j)
        {
            if (j == column) continue;
            submatrix[submatrix_row][submatrix_col] = matrix[i][j];
            ++submatrix_col;
        }
        ++submatrix_row;
    }

    return submatrix;

}

double calculate_determinant (Matrix matrix)
{
    double determinant = 0;

    double dim = matrix.size();
    assert (matrix.size() == matrix[0].size());
    assert (dim >= 1);

    if (dim == 1)
    {
        return matrix[0][0];
    }
    else if (dim == 2)
    {
        return (matrix[0][0] * matrix[1][1]) - (matrix[0][1]*matrix[1][0]);
    }

    else
    {
        for (int i = 0; i < matrix.size(); ++i)
        {
            int sign = (i % 2 == 0) ? 1 : -1;
            if (matrix[0][i] == 0) continue;
            Matrix submatrix = obtain_submatrix(matrix, 0, i);
            determinant += sign * matrix[0][i]* calculate_determinant(submatrix);
           // std::cout << "submatrix determinant is " << calculate_determinant(submatrix) << std::endl;
        }
        return determinant;
    }

}

Matrix calculate_co_factor_matrix (Matrix matrix)
{
    Matrix co_factor_matrix;
    assert (matrix.size() == matrix[0].size());
    co_factor_matrix.resize( matrix.size(),std::vector<double>( matrix[0].size() ) );

    for (int i = 0 ; i < matrix.size(); ++i)
    {
        for (int j = 0; j < matrix.size(); ++j)
        {
            int sign = ((i + j) % 2 == 0) ? 1: -1;
            co_factor_matrix[i][j] = sign * calculate_determinant(obtain_submatrix(matrix,i,j));
        }
    }


    return co_factor_matrix;
}

Matrix calculate_adjoint (Matrix matrix)
{
    Matrix adjoint;
    assert (matrix.size() == matrix[0].size());
    adjoint.resize( matrix.size(),std::vector<double>( matrix[0].size() ) );

    for (int i = 0; i < matrix.size(); ++i)
    {
        for (int j = 0; j < matrix.size(); ++j)
        {
            adjoint[j][i] = matrix[i][j];
        }
    }
    return adjoint;
}

Matrix matrix_inverse (Matrix matrix)
{
    Matrix inverse;
    assert (matrix.size() == matrix[0].size());
    inverse.resize( matrix.size(),std::vector<double>( matrix[0].size() ) );

    double tol = 1e-7;
    double matrix_determinant = calculate_determinant(matrix);
    assert(std::abs(matrix_determinant) > tol);

    Matrix co_factor_matrix = calculate_co_factor_matrix(matrix);

    Matrix inverse_unscaled = calculate_adjoint(co_factor_matrix);

    for (int i = 0; i < matrix.size(); ++i)
    {
        for (int j = 0; j < matrix.size(); ++j)
        {
            inverse[i][j] = inverse_unscaled[i][j]/matrix_determinant;
        }
    }

    return inverse;

}

}

namespace StressTransformation
{
Matrix calculate_stress_transform_off_to_on(double angle)
{
    Matrix transform_matrix;

    transform_matrix.resize(3,std::vector<double>(3));
    //convert angle to radians:
    angle = angle * PI / 180.0;

    double m = std::cos(angle);
    double n = std::sin(angle);

    transform_matrix[0][0] = m * m;
    transform_matrix[0][1] = n * n;
    transform_matrix[0][2] = 2.0 * m * n;
    transform_matrix[1][0] = n * n;
    transform_matrix[1][1] = m * m;
    transform_matrix[1][2] = -2.0 * m * n;
    transform_matrix[2][0] = -1.0*m*n;
    transform_matrix[2][1] = m*n;
    transform_matrix[2][2] = (m*m) - (n*n);

    return transform_matrix;

}

Matrix calculate_stress_transform_on_to_off(double angle)
{
    Matrix transform_matrix;

    transform_matrix.resize(3,std::vector<double>(3));

    transform_matrix = calculate_stress_transform_off_to_on(-angle);

    return transform_matrix;
}
}

namespace StrainTransformation
{
Matrix calculate_strain_transform_off_to_on(double angle)
{
    Matrix transform_matrix;

    transform_matrix.resize(3,std::vector<double>(3));
    //convert angle to radians:
    angle = angle * PI / 180.0;

    double m = std::cos(angle);
    double n = std::sin(angle);

    transform_matrix[0][0] = m * m;
    transform_matrix[0][1] = n * n;
    transform_matrix[0][2] = m * n;
    transform_matrix[1][0] = n * n;
    transform_matrix[1][1] = m * m;
    transform_matrix[1][2] = -1.0 * m * n;
    transform_matrix[2][0] = -2.0*m*n;
    transform_matrix[2][1] = 2.0*m*n;
    transform_matrix[2][2] = (m*m) - (n*n);

    return transform_matrix;

}

Matrix calculate_strain_transform_on_to_off(double angle)
{
    Matrix transform_matrix;

    transform_matrix.resize(3,std::vector<double>(3));

    transform_matrix = calculate_strain_transform_off_to_on(-angle);

    return transform_matrix;
}
}

namespace ModulusTransform
{

Matrix calculate_modulus_transform_on_to_off(double angle, Matrix Q)
{
    double U1, U2, U3, U4, U5;
    U1 = 1./8. * (3.*Q[0][0] + 3. * Q[1][1] + 2. * Q[0][1] + 4. * Q[2][2]);
    U2 = 1./2. * (Q[0][0] - Q[1][1]);
    U3 = 1./8. * (Q[0][0] + Q[1][1] - 2.* Q[0][1] - 4. * Q[2][2]);
    U4 = 1./8. * (Q[0][0] + Q[1][1] + 6.* Q[0][1] - 4. * Q[2][2]);
    U5 = 1./8. * (Q[0][0] + Q[1][1] - 2.* Q[0][1] + 4. * Q[2][2]);

    Matrix Q_out;
    Q_out.resize(3,std::vector<double>(3));

    //convert angle to radians:
    angle = angle * PI / 180.0;


    Q_out[0][0] = U1 + U2 * std::cos(2*angle) + U3 * std::cos(4*angle);
    Q_out[0][1] = U4 - U3 * std::cos(4*angle);
    Q_out[0][2] = U2 * 0.5 * std::sin(2*angle) + U3 * std::sin(4*angle);
    Q_out[1][0] = Q_out[0][1];
    Q_out[1][1] = U1 - U2 * std::cos(2*angle) + U3 * std::cos(4*angle);
    Q_out[1][2] = U2 * 0.5 * std::sin(2*angle) - U3 * std::sin(4*angle);
    Q_out[2][0] = Q_out[0][2];
    Q_out[2][1] = Q_out[1][2];
    Q_out[2][2] = U5 - U3 * std::cos(4*angle);

    return Q_out;
}

}

namespace ComplianceTransform
{

Matrix calculate_compliance_transform_on_to_off(double angle, Matrix S)
{
    double U1, U2, U3, U4, U5;
    U1 = 1./8. * (3*S[0][0] + 3 * S[1][1] + 2 * S[0][1] + S[2][2]);
    U2 = 1./2. * (S[0][0] - S[1][1]);
    U3 = 1./8. * (S[0][0] + S[1][1] - 2* S[0][1] - S[2][2]);
    U4 = 1./8. * (S[0][0] + S[1][1] + 6* S[0][1] - S[2][2]);
    U5 = 1./2. * (S[0][0] + S[1][1] - 2* S[0][1] + S[2][2]);

    Matrix S_out;
    S_out.resize(3,std::vector<double>(3));

    //convert angle to radians:
    angle = angle * PI / 180.0;


    S_out[0][0] = U1 + U2 * std::cos(2*angle) + U3 * std::cos(4*angle);
    S_out[0][1] = U4 - U3 * std::cos(4*angle);
    S_out[0][2] = U2  * std::sin(2*angle) + 2 * U3 * std::sin(4*angle);
    S_out[1][0] = S_out[0][1];
    S_out[1][1] = U1 - U2 * std::cos(2*angle) + U3 * std::cos(4*angle);
    S_out[1][2] = U2 * std::sin(2*angle) - 2 * U3 * std::sin(4*angle);
    S_out[2][0] = S_out[0][2];
    S_out[2][1] = S_out[1][2];
    S_out[2][2] = U5 - 4 * U3 * std::cos(4*angle);

    return S_out;
}

}

namespace OverallModulus
{
Matrix calculate_overall_in_plane_modulus(std::vector<AllParameters> parameters_vector, Matrix Q)
{
    double V1, V2, V3, V4;
    double total_height = 0;
    for (int i = 0; i < parameters_vector.size(); ++i)
    {
        double angle = parameters_vector[i].angle * PI/180.0;
        double h = parameters_vector[i].thickness;
        V1 += std::cos(2* angle) * h;
        V2 += std::cos(4*angle)*h;
        V3 += std::sin(2*angle)*h;
        V4 += std::sin(4*angle)*h;
        total_height += h;
    }


    double U1, U2, U3, U4, U5;
    U1 = 1./8. * (3.*Q[0][0] + 3. * Q[1][1] + 2. * Q[0][1] + 4. * Q[2][2]);
    U2 = 1./2. * (Q[0][0] - Q[1][1]);
    U3 = 1./8. * (Q[0][0] + Q[1][1] - 2.* Q[0][1] - 4. * Q[2][2]);
    U4 = 1./8. * (Q[0][0] + Q[1][1] + 6.* Q[0][1] - 4. * Q[2][2]);
    U5 = 1./8. * (Q[0][0] + Q[1][1] - 2.* Q[0][1] + 4. * Q[2][2]);

    Matrix A;
    A.resize(3,std::vector<double>(3));
    A[0][0] = (U1 * total_height) + (U2 * V1) + (U3 * V2);
    A[1][1] = (U1 * total_height) - (U2 * V1) + (U3 * V2);
    A[0][1] = (U4 * total_height) - (U3 * V2);
    A[2][2] = (U5 * total_height) - (U3 * V2);
    A[0][2] = (0.5 * U2 * V3) + (U3 * V4);
    A[1][2] = (0.5 * U2 * V3) - (U3 * V4);
    A[2][1] = A[1][2];
    A[1][0] = A[0][1];
    A[2][0] = A[0][2];

    return A;
}

Matrix calculate_overall_in_plane_compliance(std::vector<AllParameters> parameters_vector, Matrix S)
{
    Matrix A = calculate_overall_in_plane_modulus(parameters_vector,S);
    Matrix a = MatrixOperations::matrix_inverse(A);
    return a;
}

}

