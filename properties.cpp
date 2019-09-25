/*
 * properties.cpp
 *
 *  Created on: Sep 23, 2019
 *      Author: abtinameri
 */
#include "properties.h"

GeometryParameters::GeometryParameters(
                                       double thickness,
                                       double angle
                                       )
:
thickness(thickness),
angle(angle)
{}

ModulusParameters::ModulusParameters(std::string material)
{
    if (material == "T300_N520") {
        E_x = 181.0;
        E_y = 10.30;
        E_s = 7.17;
        nu_x = 0.28;
    } else if (material == "E_glass_epoxy") {
        E_x = 38.6;
        E_y = 8.27;
        E_s = 4.14;
        nu_x = 0.26;
    } else if (material == "Kev49_epoxy") {
        E_x = 76.0;
        E_y = 5.50;
        E_s = 2.30;
        nu_x = 0.34;
    } else if (material == "AS_H3501") {
        E_x = 138.0;
        E_y = 8.96;
        E_s = 7.10;
        nu_x = 0.30;
    } else if (material == "AS4_PEEK") {
        E_x = 134.0;
        E_y = 8.90;
        E_s = 5.10;
        nu_x = 0.28;
    } else {
        std::cout << "Invalid material type. Please try again" << std::endl;
    }
}

StrengthParameters::StrengthParameters(std::string material)
{
    if (material == "T300_N520") {
        X_t = 1500.;
        X_c = 1500.;
        Y_t = 40.;
        Y_c = 246.;
        S_c = 68.;
    } else if (material == "E_glass_epoxy") {
        X_t = 1062.;
        X_c = 610.;
        Y_t = 31.;
        Y_c = 118.;
        S_c = 72.;
    } else if (material == "Kev49_epoxy") {
        X_t = 1400.;
        X_c = 235.;
        Y_t = 12.;
        Y_c = 53.;
        S_c = 34.;
    } else if (material == "AS_H3501") {
        X_t = 1447.;
        X_c = 1447.;
        Y_t = 51.7;
        Y_c = 206.;
        S_c = 93.;
    } else if (material == "AS4_PEEK") {
        X_t = 2130.;
        X_c = 1100.;
        Y_t = 80.;
        Y_c = 200.;
        S_c = 160.;
    } else {
        std::cout << "Invalid material type. Please try again." << std::endl;
    }
}

AllParameters::AllParameters(
                             double thickness,
                             double angle,
                             std::string material)
:
        GeometryParameters(thickness,angle),
        ModulusParameters(material),
        StrengthParameters(material),
        material(material)
{

}




