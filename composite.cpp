/*
 * composite.cpp
 *
 *  Created on: Sep 15, 2019
 *      Author: abtinameri
 */

#include "composite.h"
#include "utilities.h"
#include "properties.h"

Composite::Composite(AllParameters parameters_input)
:
all_parameters(parameters_input)
{
    on_axis_S_matrix.resize(3,std::vector<double>(3));
    on_axis_Q_matrix.resize(3,std::vector<double>(3));
    off_axis_S_matrix.resize(3,std::vector<double>(3));
    off_axis_Q_matrix.resize(3,std::vector<double>(3));

    on_axis_S_matrix = calculate_on_axis_S_matrix();
    on_axis_Q_matrix = calculate_on_axis_Q_matrix();

    //need to change this later
    off_axis_Q_matrix = ModulusTransform::calculate_modulus_transform_on_to_off(all_parameters.angle,on_axis_Q_matrix);
    off_axis_S_matrix = ComplianceTransform::calculate_compliance_transform_on_to_off(all_parameters.angle,on_axis_S_matrix);

}

std::vector<std::vector<double>> Composite::calculate_on_axis_S_matrix()
{
	double S_xx, S_xy, S_yy, S_ss;


	S_xx = 1./all_parameters.E_x;
	S_xy = -all_parameters.nu_x/all_parameters.E_x;
	S_yy = 1./all_parameters.E_y;
	S_ss = 1./all_parameters.E_s;


    on_axis_S_matrix[0][0] = S_xx;
    on_axis_S_matrix[0][1] = S_xy;
    on_axis_S_matrix[0][2] = 0;
    on_axis_S_matrix[1][0] = S_xy;
    on_axis_S_matrix[1][1] = S_yy;
    on_axis_S_matrix[1][2] = 0;
    on_axis_S_matrix[2][0] = 0;
    on_axis_S_matrix[2][1] = 0;
    on_axis_S_matrix[2][2] = S_ss;



//	std::cout << "The S matrix is (GPa):" << std::endl
//			  << S_xx << " | " << S_xy << " | " << 0 << std::endl
//			  << S_xy << " | " << S_yy << " | " << 0 << std::endl
//			  <<  0 << " | " << 0 << " | " << S_ss << std::endl
//			  << "--------------------------------------------" << std::endl;
//
//	std::cout << "The Q matrix is (1/GPa):" << std::endl
//				  << Q_xx << " | " << Q_xy << " | " << 0 << std::endl
//				  << Q_xy << " | " << Q_yy << " | " << 0 << std::endl
//				  <<  0 << " | " << 0 << " | " << Q_ss << std::endl
//				  << "--------------------------------------------" << std::endl;

    return on_axis_S_matrix;
}

std::vector<std::vector<double>> Composite::calculate_on_axis_Q_matrix()
{
    double Q_xx, Q_xy, Q_yy, Q_ss;

    double m = 1./(1.-(all_parameters.nu_x * all_parameters.nu_x * all_parameters.E_y/all_parameters.E_x));

    Q_xx = m * all_parameters.E_x;
    Q_xy = m * all_parameters.nu_x * all_parameters.E_y;
    Q_yy = m * all_parameters.E_y;
    Q_ss = all_parameters.E_s;

    on_axis_Q_matrix[0][0] = Q_xx;
    on_axis_Q_matrix[0][1] = Q_xy;
    on_axis_Q_matrix[0][2] = 0;
    on_axis_Q_matrix[1][0] = Q_xy;
    on_axis_Q_matrix[1][1] = Q_yy;
    on_axis_Q_matrix[1][2] = 0;
    on_axis_Q_matrix[2][0] = 0;
    on_axis_Q_matrix[2][1] = 0;
    on_axis_Q_matrix[2][2] = Q_ss;

    return on_axis_Q_matrix;
}

void Composite::output_all(int layer_number)
{
    std::cout << "Layer number: " << layer_number << std::endl;

	std::cout << "The material used is " << all_parameters.material << std::endl;

	std::cout << "Modulus parameters:" << std::endl
			  << "E_x = " << all_parameters.E_x << "GPa" << std::endl
			  << "E_y = " << all_parameters.E_y << "GPa" << std::endl
			  << "E_s = " << all_parameters.E_s << "GPa" << std::endl
			  << "nu_x = " << all_parameters.nu_x << std::endl
			  << "--------------------------------------------" << std::endl;

    std::cout << "Strength parameters:" << std::endl
			  << "X_t = " << all_parameters.X_t << "MPa" << std::endl
			  << "X_c = " << all_parameters.X_c<< "MPa" << std::endl
			  << "Y_t = " << all_parameters.Y_t << "MPa" << std::endl
			  << "Y_c = " << all_parameters.Y_c << "MPa" << std::endl
			  << "S_c = " << all_parameters.S_c << "MPa" << std::endl
			  << "--------------------------------------------" << std::endl;

    std::cout << "Geometric parameters:" << std::endl;
    std::cout << "theta = " << all_parameters.angle << " degrees" << std::endl;

    std::cout << "ply thickness = " << all_parameters.thickness << "mm" << std::endl;
    std::cout << "--------------------------------------------" << std::endl;

    MatrixOperations::output_vector(off_axis_stress,"off-axis stress","GPa");
    MatrixOperations::output_vector(on_axis_stress,"on-axis stress","GPa");
    MatrixOperations::output_vector(on_axis_strain,"on-axis strain","");

    //calculate_on_axis_S_matrix();
    //calculate_on_axis_Q_matrix();

    //ModulusTransform::calculate_modulus_transform_on_to_off(all_parameters.angle,on_axis_Q_matrix);
    //ComplianceTransform::calculate_compliance_transform_on_to_off(all_parameters.angle,on_axis_S_matrix);

   // MatrixOperations::output_matrix(on_axis_S_matrix, "on-axis S", "1/GPa");
   // MatrixOperations::output_matrix(on_axis_Q_matrix, "on-axis Q", "GPa");

    //MatrixOperations::output_matrix(off_axis_S_matrix, "off-axis S", "1/GPa");
    //MatrixOperations::output_matrix(off_axis_Q_matrix, "off-axis Q", "GPa");
}

void Composite::calculate_stresses_and_strains (std::vector<double> off_axis_strain)
{

    std::vector<std::vector<double>> stress_transform_off_to_on = StressTransformation::calculate_stress_transform_off_to_on(all_parameters.angle);
    std::vector<std::vector<double>> strain_transform_off_to_on = StrainTransformation::calculate_strain_transform_off_to_on(all_parameters.angle);

    on_axis_strain = MatrixOperations::vector_mult(strain_transform_off_to_on, off_axis_strain);
    off_axis_stress = MatrixOperations::vector_mult(off_axis_Q_matrix, off_axis_strain);
    on_axis_stress = MatrixOperations::vector_mult(on_axis_Q_matrix, off_axis_strain);

}


