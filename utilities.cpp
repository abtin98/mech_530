/*
 * utilities.cpp
 *
 *  Created on: Sep 23, 2019
 *      Author: abtinameri
 */


#include "utilities.h"

namespace MatrixOperations
{
void output_matrix(std::vector<std::vector<double>> matrix, std::string name, std::string units)
{
    std::cout << "The " << name << " matrix is (" << units << "):" <<  std::endl;
    for (int i = 0; i < matrix.size() ; ++i)
    {
        for (int j = 0; j < matrix[0].size(); ++j)
        {
            std::cout << matrix[i][j] << " | ";
        }
        std::cout << std::endl;
    }
    std::cout << "--------------------------------------------" << std::endl;
}

}

