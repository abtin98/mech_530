/*
 * composite.h
 *
 *  Created on: Sep 15, 2019
 *      Author: abtinameri
 */

#ifndef COMPOSITE_H_
#define COMPOSITE_H_


#include "properties.h"
#include <array>

class Composite
{
public:
	Composite();
	Composite(AllParameters parameters_input);
    using Matrix = std::vector<std::vector<double>>;
	Matrix calculate_on_axis_S_matrix();
	Matrix calculate_on_axis_Q_matrix();
	//Matrix calculate_off_axis_S_matrix(Matrix on_axis_S_matrix, );
	//Matrix calculate_off_axis_Q_matrix(Matrix on_axis_Q_matrix);

	void output_all(int layer_number);

    AllParameters all_parameters;
    Matrix on_axis_S_matrix;
    Matrix on_axis_Q_matrix;
    Matrix off_axis_S_matrix;
    Matrix off_axis_Q_matrix;

    std::vector<double> off_axis_stress;
    std::vector<double> on_axis_stress;
    std::vector<double> on_axis_strain;

    void calculate_stresses_and_strains (std::vector<double> off_axis_strain);


private:

};


#endif /* COMPOSITE_H_ */
