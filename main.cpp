#include <iostream>
#include "properties.h"
#include "composite.h"

int main(int argc, char **argv) {


	int n = 14;
	std::vector<double> angle = {0,0,-20,20,90,10,-10,-10,10,90,20,-20,0,0};
	std::vector<double> thickness = {0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125};
	double z_c = 0.0;

	AllParameters all_parameters(n,angle,thickness,z_c,"AS_H3501");


	Composite composite(all_parameters,"AS/H3501");

	composite.output_all();

}
