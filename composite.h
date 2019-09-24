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
	Composite() = delete;
	Composite(AllParameters parameters_input);
	void calculate_S_Q_matrices();
	void output_all();


private:
	AllParameters all_parameters;
	std::vector<std::vector<double>> S_matrix;
	std::vector<std::vector<double>> Q_matrix;
};


#endif /* COMPOSITE_H_ */
