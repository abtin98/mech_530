/*
 * laminate.cpp
 *
 *  Created on: Sep 24, 2019
 *      Author: abtinameri
 */

#include "laminate.h"
#include <assert.h>
#include "utilities.h"

Laminate::Laminate(std::vector<AllParameters> all_parameters_vector, double z_c, std::vector<double> applied_stress_vector)
:
z_c(z_c),
all_parameters_input(all_parameters_vector),
applied_stress_vector(applied_stress_vector)
{

    for (int i = 0; i < all_parameters_vector.size(); ++i)
    {
        Composite composite(all_parameters_input[i]);
        layer.push_back(composite);
    }

    A = OverallModulus::calculate_overall_in_plane_modulus(all_parameters_input, layer[0].on_axis_S_matrix);
    a = MatrixOperations::matrix_inverse(A);

    off_axis_strain = MatrixOperations::vector_mult (a, applied_stress_vector);

    for (int i = 0; i < all_parameters_vector.size(); ++i)
    {
        layer[i].calculate_stresses_and_strains(off_axis_strain);
    }

}

void Laminate::output_all()
{
    std::cout << "OUTPUT" << std::endl;
    std::cout << "z_c = " << z_c << std::endl;
    std::cout << "Number of layers = " << all_parameters_input.size() << std::endl;

    MatrixOperations::output_matrix(A,"[A] matrix", "N/m");
    MatrixOperations::output_matrix(a,"[a] matrix", "m/N");

    MatrixOperations::output_vector (applied_stress_vector, "applied stress vector", "MPa-m");

    MatrixOperations::output_vector(off_axis_strain, "off-axis strain", "");

    for (int i = 0; i < all_parameters_input.size(); ++i)
    {
        layer[i].output_all(i);
    }
}


