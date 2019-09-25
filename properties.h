/*
 * properties.h
 *
 *  Created on: Sep 15, 2019
 *      Author: abtinameri
 */

#ifndef PROPERTIES_H_
#define PROPERTIES_H_


#include <vector>
#include <iostream>

class GeometryParameters
{
public:
    GeometryParameters() = delete;
    GeometryParameters(
                       double thickness,
                       double angle
                       );
    double thickness;
	double angle;
};

class ModulusParameters
{
public:
    ModulusParameters() = delete;
    ModulusParameters(std::string material);
	double E_x;
	double E_y;
	double E_s;
	double nu_x;
};

class StrengthParameters
{
public:
    StrengthParameters() = delete;
    StrengthParameters(std::string material);
	double X_t;
	double X_c;
	double Y_t;
	double Y_c;
	double S_c;
};

class AllParameters : public GeometryParameters, public ModulusParameters, public StrengthParameters
{
public:
    AllParameters();
    std::string material;
    AllParameters(
                  double thickness,
                  double angle,
                  std::string material);
};


#endif /* PROPERTIES_H_ */
