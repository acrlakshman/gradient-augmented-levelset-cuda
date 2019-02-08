/*
This file defines the functions required to calculate the analytical velocity at
any point in space and time for the given 2D domain
This file is for a Vortex type velocity field. For any other field, this file needs
to be modified completely.
*/

#ifndef _VortexVelocity_h
#define _VortexVelocity_h

#include "Constants.h"

namespace galsfunctions
{
    
    double Velx(double x, double y, double t, double T_period){
        double temp;
        temp = pow(sin(pi * x),2) * sin(2 * pi * y) * cos(pi * t/T_period);
        return temp;
    }
    
    double Vely(double x, double y,double t, double T_period){
        double temp;
        temp = -pow(sin(pi * y),2) * sin(2 * pi * x) * cos(pi * t/T_period);
        return temp;
    }
    
    double gradUx(double x, double y, double t, double T_period){
        double temp;
        temp = pi * sin(2 * pi * y) * sin(2 * pi * x) * cos(pi * t/T_period);
        return temp;
    }
    
    double gradUy(double x, double y, double t, double T_period){
        double temp;
        temp = 2 * pi * pow(sin(pi * x),2) * cos(2 * pi * y) * cos(pi * t/T_period);
        return temp;
    }
    
    double gradVx(double x, double y, double t, double T_period){
        double temp;
        temp = -2 * pi * pow(sin(pi * y),2) * cos(2 * pi * x) * cos(pi * t/T_period);
        return temp;
    }
    
    double gradVy(double x, double y, double t, double T_period){
        double temp;
        temp = -pi * sin(2 * pi * y) * sin(2 * pi * x) * cos(pi * t/T_period);
        return temp;
    }
}

#endif
