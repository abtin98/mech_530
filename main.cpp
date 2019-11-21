#include <iostream>
#include "properties.h"
#include "composite.h"
#include "utilities.h"
#include "laminate.h"

int main(int argc, char **argv) {

    { //example 1

    std::cout << "EXAMPLE 1:" << std::endl;
	double thickness = 0.125;
	double optimal_theta;
	double p = 5000000.0; //N/m^2
	double D = 0.08; //mm
	double N1, N2;
	N1 = p * D /4.0;
	N2 = p * D /2.0;
    double max_safety_factor = 0;


	for (double theta = 54; theta < 55; ++theta)
	{
	    std::vector<double> applied_stress_vector = {N1/1000.0, N2/1000.0, 0.0}; //kN/mm
	    std::vector<double> applied_moment_vector = {0.0, 0.0, 0.0}; //kN
	    std::vector<double> angle = {+theta, -theta,
	                                 +theta, -theta,
	                                 +theta, -theta,
	                                 +theta, -theta,
	                                 -theta, +theta,
	                                 -theta, +theta,
	                                 -theta, +theta,
	                                 -theta, +theta};
	    int n = angle.size();
	    std::vector<AllParameters> all_parameters_vector;

	    for (int i = 0; i < n; ++i)
	    {
	        AllParameters all_parameters(thickness,angle[i],"T300_N520");
	        all_parameters_vector.push_back(all_parameters);
	    }

	    double z_c = 5.0;
	    Laminate laminate (all_parameters_vector, z_c, applied_stress_vector, applied_moment_vector);
	    //laminate.output_all();
	    double safety_factor = laminate.find_min_safety_factor(true);

	    if (safety_factor > max_safety_factor)
	    {
	        optimal_theta = theta;
	        max_safety_factor = safety_factor;
	    }
	    //std::cout << laminate.find_min_safety_factor(false) << std::endl;

	}

	std::cout << "optimal angle is " << optimal_theta << " and the min safety factor is "
	          << max_safety_factor <<std::endl;

    }
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    {//example 2
        double theta = 0.0;
        std::cout << "EXAMPLE 2:" << std::endl;
        double thickness = 0.125;
        std::vector<double> applied_stress_vector = {5.000, 0.0, 0.0}; //kN/mm
        std::vector<double> applied_moment_vector = {-1250, 0.0, 0.0}; //kN
        std::vector<double> angle = {50.0,
                                     20.0,
                                    -70.0,
                                     50.0,
                                     0.0,
                                     0.0,
                                     50.0,
                                     -70.0,
                                     20.0,
                                     50.0};

        int n = angle.size();
        std::vector<AllParameters> all_parameters_vector;

        for (int i = 0; i < n; ++i)
        {
            AllParameters all_parameters(thickness,angle[i],"Kev49_epoxy");
            all_parameters_vector.push_back(all_parameters);
        }
        double z_c = 0.0;
        Laminate laminate (all_parameters_vector, z_c, applied_stress_vector, applied_moment_vector);
        laminate.output_all();
        std::cout << laminate.find_min_safety_factor(false) << std::endl;

    }
//    std::cout << std::endl;
//    std::cout << std::endl;
//    std::cout << std::endl;
//    {//example 3a
//        std::cout << "EXAMPLE 3a:" << std::endl;
//        double  q, C , L, W, N1, M1, M2;
//
//        N1 = 12.0000;
//        q = 92000;
//        C = 140;
//        L = 0.48;
//        W = 0.27;
//
//        M1 = q * W * W/ 8.0 + C;
//        M2 = q * L * L/ 8.0;
//
//        double thickness = 0.125;
//        std::vector<double> applied_stress_vector = {N1, 0.0, 0.0}; //kN/mm
//        std::vector<double> applied_moment_vector = {M1, M2, 0.0}; //kN
//        std::vector<double> angle = {0.0, 90.0, 90.0, 30.0, -30.0, 45.0, -45.0, 70.0, -70.0
//                                     -70.0, 70.0, -45.0, 45.0, -30.0, 30.0, 90.0, 90.0, 0.0};
//        int n = angle.size();
//        std::vector<AllParameters> all_parameters_vector;
//
//        for (int i = 0; i < n; ++i)
//        {
//            AllParameters all_parameters(thickness,angle[i],"AS4_PEEK");
//            all_parameters_vector.push_back(all_parameters);
//        }
//
//        double z_c = 5.0;
//        Laminate laminate (all_parameters_vector, z_c, applied_stress_vector, applied_moment_vector);
//        laminate.output_all();
//        std::cout << "min safety factor is " << laminate.find_min_safety_factor(true) << std::endl;
//    }
//    std::cout << std::endl;
//    std::cout << std::endl;
//    std::cout << std::endl;
//    {//example 3b
//            std::cout << "EXAMPLE 3b:" << std::endl;
//            double  q, C , L, W, N1, M1, M2;
//
//            N1 = -12.0000; // N/m
//            q = 92000; // N/m^2
//            C = 140; // Nm/m
//            L = 0.48; // m
//            W = 0.27; // m
//
//            M1 = q * W * W/ 8.0 + C;
//            M2 = q * L * L/ 8.0;
//
//            double thickness = 0.125;
//            std::vector<double> applied_stress_vector = {N1, 0.0, 0.0};//{5.000, 0.0, 0.0}; //kN/mm
//            std::vector<double> applied_moment_vector = {M1, M2, 0.0}; //{-1250, 0, 0.0}; //kN
//            std::vector<double> angle = {0.0, 90.0, 90.0, 30.0, -30.0, 45.0, -45.0, 70.0, -70.0
//                                        -70.0, 70.0, -45.0, 45.0, -30.0, 30.0, 90.0, 90.0, 0.0};
//            int n = angle.size();
//            std::vector<AllParameters> all_parameters_vector;
//
//            for (int i = 0; i < n; ++i)
//            {
//                AllParameters all_parameters(thickness,angle[i],"AS4_PEEK");
//                all_parameters_vector.push_back(all_parameters);
//            }
//
//            double z_c = 5.0;
//            Laminate laminate (all_parameters_vector, z_c, applied_stress_vector, applied_moment_vector);
//            laminate.output_all();
//              std::cout << "min safety factor is "  << laminate.find_min_safety_factor(false) << std::endl;
//
//      }
}
