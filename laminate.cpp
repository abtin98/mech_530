/*
 * laminate.cpp
 *
 *  Created on: Sep 24, 2019
 *      Author: abtinameri
 */

#include "laminate.h"
#include <assert.h>

Laminate::Laminate(std::vector<AllParameters> all_parameters_vector, double z_c)
:
z_c(z_c),
all_parameters_input(all_parameters_vector)
{
    for (int i = 0; i < all_parameters_vector.size(); ++i)
    {
        Composite composite(all_parameters_input[i]);
        layer.push_back(composite);
    }
}

void Laminate::output_all()
{
    std::cout << "OUTPUT" << std::endl;
    std::cout << "z_c = " << z_c << std::endl;
    std::cout << "Number of layers = " << all_parameters_input.size() << std::endl;
    for (int i = 0; i < all_parameters_input.size(); ++i)
    {
        layer[i].output_all(i);
    }
}


