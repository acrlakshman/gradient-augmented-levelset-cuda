/*
 * This file is called by the Allocation.h to calculate the level set
 * and its gradient values at specific points on the grid for initialization.
 */

const double pi = M_PI;
const double xo = .5;
const double yo = .75;
const double rcircle = .15;

using namespace std;

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
