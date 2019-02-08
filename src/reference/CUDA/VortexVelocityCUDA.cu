/*
This file defines the functions required to calculate the analytical velocity at
any point in space and time for the given 2D domain
This file is for a Vortex type velocity field. For any other field, this file needs
to be modified completely.
*/

#ifndef _VortexVelocityCUDA_cu
#define _VortexVelocityCUDA_cu

__device__ __host__ double Velx(double x, double y, double t, double T_period){
    double temp;
    temp = pow(sin(pi * x),2) * sin(2 * pi * y) * cos(pi * t/T_period);
    return temp;
}

__device__ __host__ double Vely(double x, double y,double t, double T_period){
    double temp;
    temp = -pow(sin(pi * y),2) * sin(2 * pi * x) * cos(pi * t/T_period);
    return temp;
}

__device__ double gradUx(double x, double y, double t, double T_period){
    double temp;
    temp = pi * sin(2 * pi * y) * sin(2 * pi * x) * cos(pi * t/T_period);
    return temp;
}

__device__ double gradUy(double x, double y, double t, double T_period){
    double temp;
    temp = 2 * pi * pow(sin(pi * x),2) * cos(2 * pi * y) * cos(pi * t/T_period);
    return temp;
}

__device__ double gradVx(double x, double y, double t, double T_period){
    double temp;
    temp = -2 * pi * pow(sin(pi * y),2) * cos(2 * pi * x) * cos(pi * t/T_period);
    return temp;
}

__device__ double gradVy(double x, double y, double t, double T_period){
    double temp;
    temp = -pi * sin(2 * pi * y) * sin(2 * pi * x) * cos(pi * t/T_period);
    return temp;
}
#endif
