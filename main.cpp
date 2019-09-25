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
	std::vector<double> angle = {0.0 , 90.0 , 20.0 , -20.0, 45.0, -45.0, -45.0, 45.0, -20.0, 20.0, 90.0, 0.0};
	int n = angle.size();
	std::vector<AllParameters> all_parameters_vector;

	for (int i = 0; i < n; ++i)
	{
	    AllParameters all_parameters(thickness,angle[i],"E_glass_epoxy");
	    all_parameters_vector.push_back(all_parameters);
	}

	double z_c = 0.0;
    Laminate laminate (all_parameters_vector, z_c);
    laminate.output_all();



    }


    { //Example 2
        std::cout << "Example 2:" << std::endl;
        double thickness = 0.125;
        std::vector<double> angle = {30.0};


        int n = angle.size();
        std::vector<AllParameters> all_parameters_vector;

        for (int i = 0; i < n; ++i)
        {
            AllParameters all_parameters(thickness,angle[i],"E_glass_epoxy");
            all_parameters_vector.push_back(all_parameters);
        }

        double z_c = 0.0;

        Laminate laminate (all_parameters_vector, z_c);

        laminate.output_all();

        std::vector<double> off_axis_stress = {4.0, -0.860, -0.990};
        std::vector<double> off_axis_strain = MatrixOperations::vector_mult (laminate.layer[0].off_axis_S_matrix,off_axis_stress);

        std::vector<std::vector<double>> stress_transform_off_to_on = StressTransformation::calculate_stress_transform_off_to_on(angle[0]);
        std::vector<std::vector<double>> strain_transform_off_to_on = StrainTransformation::calculate_strain_transform_off_to_on(angle[0]);

        std::vector<double> on_axis_stress = MatrixOperations::vector_mult(stress_transform_off_to_on, off_axis_stress);
        std::vector<double> on_axis_strain = MatrixOperations::vector_mult(strain_transform_off_to_on, off_axis_strain);

        MatrixOperations::output_vector(off_axis_stress,"off-axis stress","GPa");
        MatrixOperations::output_vector(off_axis_strain,"off-axis strain","");
        MatrixOperations::output_vector(on_axis_stress,"on-axis stress","GPa");
        MatrixOperations::output_vector(on_axis_strain,"on-axis strain","");

    }

}
