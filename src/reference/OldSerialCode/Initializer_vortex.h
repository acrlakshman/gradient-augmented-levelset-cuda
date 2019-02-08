//
//  Initializer.h
//
//
//  Created by Raunak Bardia on 10/10/14.
//
//

#ifndef _Initializer_h
#define _Initializer_h

#include <math.h>
#include <stdio.h>
#include <string.h>

const double pi = M_PI;
const double xo = .5;
const double yo = .75;
const double rcircle = .15;

using namespace std;

double Velx(double x, double y, double t, double T_period){
    double temp;
    double signum = cos(pi*t/T_period);
    temp = pow(sin(pi * x),2) * sin(2 * pi * y) * cos(pi * t/T_period);
    return temp;
}

double Vely(double x, double y,double t, double T_period){
    double temp;
    double signum = cos(pi*t/T_period);
    temp = -pow(sin(pi * y),2) * sin(2 * pi * x) * cos(pi * t/T_period);
    return temp;
}

double gradUx(double x, double y, double t, double T_period){
    double temp;
    double signum = cos(pi*t/T_period);
    temp = pi * sin(2 * pi * y) * sin(2 * pi * x) * cos(pi * t/T_period);
    return temp;
}

double gradUy(double x, double y, double t, double T_period){
    double temp;
    double signum = cos(pi*t/T_period);
    temp = 2 * pi * pow(sin(pi * x),2) * cos(2 * pi * y) * cos(pi * t/T_period);
    return temp;
}

double gradVx(double x, double y, double t, double T_period){
    double temp;
    double signum = cos(pi*t/T_period);
    temp = -2 * pi * pow(sin(pi * y),2) * cos(2 * pi * x) * cos(pi * t/T_period);
    return temp;
}

double gradVy(double x, double y, double t, double T_period){
    double temp;
    double signum = cos(pi*t/T_period);
    temp = -pi * sin(2 * pi * y) * sin(2 * pi * x) * cos(pi * t/T_period);
    return temp;
}
double initialize(double x, double y){
    double k;
    k = exp(-(pow((x - xo),2) + pow((y - yo),2))) - exp(-pow(rcircle,2));
    return k;
}

double derivxinit(double x, double y){
    double k;
    k = -2*(x-xo)*exp(-(pow((x - xo),2) + pow((y - yo),2)));
    return k;
}

double derivyinit(double x, double y){
    double k;
    k = -2*(y-yo)*exp(-(pow((x - xo),2) + pow((y - yo),2)));
    return k;
}

double derivxyinit(double x, double y){
    double k;
    k = 4 * (y - yo) * (x - xo) * exp(-(pow((x - xo),2) + pow((y - yo),2)));
    return k;
}


#endif