/*
 * This file is called by the Allocation.h to calculate the level set
 * and its gradient values at specific points on the grid for initialization.
 */

#ifndef _InitializeLevelSet_h
#define _InitializeLevelSet_h

#include "Constants.h"

namespace galsfunctions
{
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
}

#endif
