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
    Laminate(std::vector<AllParameters> all_parameters_vector, double z_c, std::vector<double> applied_stress_vector, std::vector<double> applied_moment_vector);
    double z_c;
    std::vector<double> applied_stress_vector;
    std::vector<double> applied_moment_vector;
    std::vector<Composite> layer;
    std::vector<AllParameters> all_parameters_input;
    void output_all();

    double perform_failure_analysis(bool is_positive);

    using Matrix = std::vector<std::vector<double>>;
    Matrix A;
    Matrix a;
    Matrix D;
    Matrix d;

    std::vector<double> off_axis_in_plane_strain;
    std::vector<double> curvature;
    std::vector<double> ply_z_coordinate;

    double find_min_safety_factor(bool isPositive);

};

#endif /* LAMINATE_H_ */
