/*
 * laminate.h
 *
 *  Created on: Sep 24, 2019
 *      Author: abtinameri
 */

#ifndef LAMINATE_H_
#define LAMINATE_H_

#include "composite.h"
#include "properties.h"

class Laminate
{
public:
    Laminate();
    Laminate(std::vector<AllParameters> all_parameters_vector, double z_c, std::vector<double> applied_stress_vector);
    double z_c;
    std::vector<double> applied_stress_vector;
    std::vector<Composite> layer;
    std::vector<AllParameters> all_parameters_input;
    void output_all();

    using Matrix = std::vector<std::vector<double>>;
    Matrix A;
    Matrix a;

    std::vector<double> off_axis_strain;

};

#endif /* LAMINATE_H_ */
