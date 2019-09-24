/*
 * composite.cpp
 *
 *  Created on: Sep 15, 2019
 *      Author: abtinameri
 */

#include "composite.h"
#include "utilities.h"

Composite::Composite(AllParameters parameters_input)
:
composite_type(composite_type),
all_parameters(parameters_input)
{
    S_matrix.resize(3,std::vector<double>(3));
    Q_matrix.resize(3,std::vector<double>(3));
}

void Composite::calculate_S_Q_matrices()
{
	double S_xx, S_xy, S_yy, S_ss;
	double Q_xx, Q_xy, Q_yy, Q_ss;

	S_xx = 1./all_parameters.E_x;
	S_xy = -all_parameters.nu_x/all_parameters.E_x;
	S_yy = 1./all_parameters.E_y;
	S_ss = 1./all_parameters.E_s;

	double m = 1./(1.-(all_parameters.nu_x * all_parameters.nu_x * all_parameters.E_y/all_parameters.E_x));

	Q_xx = m * all_parameters.E_x;
	Q_xy = m * all_parameters.nu_x * all_parameters.E_y;
	Q_yy = m * all_parameters.E_y;
	Q_ss = all_parameters.E_s;

	S_matrix[0][0] = S_xx;
    S_matrix[0][1] = S_xy;
    S_matrix[0][2] = 0;
    S_matrix[1][0] = S_xy;
    S_matrix[1][1] = S_yy;
    S_matrix[1][2] = 0;
    S_matrix[2][0] = 0;
    S_matrix[2][1] = 0;
    S_matrix[2][2] = S_ss;

    Q_matrix[0][0] = Q_xx;
    Q_matrix[0][1] = Q_xy;
    Q_matrix[0][2] = 0;
    Q_matrix[1][0] = Q_xy;
    Q_matrix[1][1] = Q_yy;
    Q_matrix[1][2] = 0;
    Q_matrix[2][0] = 0;
    Q_matrix[2][1] = 0;
    Q_matrix[2][2] = Q_ss;

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

}

void Composite::output_all()
{
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
    std::cout << "n = " << all_parameters.n << std::endl;
    std::cout << "theta = [";
    for (int i = 0; i < all_parameters.n; ++i)
    {
    	std::cout << all_parameters.angle[i] << " ";
    }
    std::cout << "]" << " degrees" <<std::endl;

    std::cout << "ply thickness = [";
        for (int i = 0; i < all_parameters.n; ++i)
        {
        	std::cout << all_parameters.thickness[i] << " ";
        }
    std::cout << "]" << "mm" <<std::endl;

    std::cout << "z_c = " << all_parameters.z_c << "mm" << std::endl;
    std::cout << "--------------------------------------------" << std::endl;

    calculate_S_Q_matrices();

    MatrixOperations::output_matrix(S_matrix, "S", "GPa");
    MatrixOperations::output_matrix(Q_matrix, "Q", "1/GPa");







}


