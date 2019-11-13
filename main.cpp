#include <iostream>
#include "properties.h"
#include "composite.h"
#include "utilities.h"
#include "laminate.h"

int main(int argc, char **argv) {



	//std::vector<double> angle = {0,0,-20,20,90,10,-10,-10,10,90,20,-20,0,0};
	//std::vector<double> thickness = {0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125};
	//double z_c = 0.0;
    { //example 1

    std::cout << "Example 1:" << std::endl;
	double thickness = 0.125;
	std::vector<double> applied_stress_vector = {5.000, 0.0, 0.0}; //kN/mm
	std::vector<double> applied_moment_vector = {-1250, 0.0, 0.0}; //kN
	std::vector<double> angle = {0.0, 0.0, 90.0, 30.0, -30.0, -30.0, 30.0, 90.0, 0.0, 0.0};
	int n = angle.size();
	std::vector<AllParameters> all_parameters_vector;

	for (int i = 0; i < n; ++i)
	{
	    AllParameters all_parameters(thickness,angle[i],"AS4_PEEK");
	    all_parameters_vector.push_back(all_parameters);
	}

	double z_c = 5.0;
    Laminate laminate (all_parameters_vector, z_c, applied_stress_vector, applied_moment_vector);
    laminate.output_all();

    }
}
