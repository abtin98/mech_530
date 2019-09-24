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
    GeometryParameters(int n,
                       std::vector<double> thickness,
                       std::vector<double> angle,
                       double z_c);
	int n;
	std::vector<double> thickness;
	std::vector<double> angle;
	double z_c;
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
    AllParameters() = delete;
    std::string material;
    AllParameters(int n,
                  std::vector<double> thickness,
                  std::vector<double> angle,
                  double z_c,
                  std::string material);
};


#endif /* PROPERTIES_H_ */
