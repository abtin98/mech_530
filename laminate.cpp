/*
 * laminate.cpp
 *
 *  Created on: Sep 24, 2019
 *      Author: abtinameri
 */

#include "laminate.h"
#include <assert.h>
#include "utilities.h"

Laminate::Laminate(std::vector<AllParameters> all_parameters_vector, double z_c, std::vector<double> applied_stress_vector, std::vector<double> applied_moment_vector)
:
z_c(z_c),
all_parameters_input(all_parameters_vector),
applied_stress_vector(applied_stress_vector),
applied_moment_vector(applied_moment_vector)
{
    double height = - z_c;

    for (int i = 0; i < all_parameters_vector.size()/2; ++i)
    {
        height-=all_parameters_vector[i].thickness;
    }

    for (int i = 0; i < all_parameters_vector.size(); ++i)
    {
        Composite composite(all_parameters_input[i]);
        layer.push_back(composite);

    }

    A = OverallModulus::calculate_overall_in_plane_modulus(all_parameters_input, layer[0].on_axis_Q_matrix);
    D = FlexuralModulus::calculate_flexural_modulus(all_parameters_input,layer[0].on_axis_Q_matrix,z_c);
    a = MatrixOperations::matrix_inverse(A);
    d = MatrixOperations::matrix_inverse(D);

    ply_z_coordinate.resize(all_parameters_vector.size());

    for (int i = 0; i < all_parameters_vector.size(); ++i)
    {
        if (i != all_parameters_vector.size()/2.0 - 1)
        {
            ply_z_coordinate[i] = height;
            height += all_parameters_vector[i].thickness;
        }

        else if (i == all_parameters_vector.size()/2.0 - 1)
        {
            ply_z_coordinate[i] = height;
            height += all_parameters_vector[i].thickness + 2*z_c;
        }
    }

    off_axis_in_plane_strain = MatrixOperations::vector_mult (a, applied_stress_vector);
    curvature = MatrixOperations::vector_mult(d, applied_moment_vector);

    for (int i = 0; i < all_parameters_vector.size(); ++i)
    {
        layer[i].calculate_stresses_and_strains(off_axis_in_plane_strain, curvature, ply_z_coordinate[i]);
    }
}

void Laminate::output_all()
{
    std::cout << "OUTPUT" << std::endl;
    std::cout << "z_c = " << z_c << "mm" << std::endl;
    std::cout << "Number of layers = " << all_parameters_input.size() << std::endl;

    MatrixOperations::output_matrix(A,"[A] matrix", "kN/mm");
    MatrixOperations::output_matrix(a,"[a] matrix", "mm/kN");
    MatrixOperations::output_matrix(D,"[D] matrix", "?????????");
    MatrixOperations::output_matrix(d,"[d] matrix", "?????????");

    MatrixOperations::output_vector (applied_stress_vector, "applied stress vector", "kN/mm");
    MatrixOperations::output_vector (applied_moment_vector, "applied moment vector", "kN-mm/mm");

    MatrixOperations::output_vector(off_axis_in_plane_strain, "off-axis strain", "");

    for (int i = 0; i < all_parameters_input.size(); ++i)
    {
        layer[i].output_all(i);
    }
}


