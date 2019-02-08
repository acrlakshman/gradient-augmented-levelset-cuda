/*
 * Parameters needed from the user to setup the case.
 * This file is only used in the MPI implementation.
 * Ideally, this will be an input argument file, which
 * is only given at the runtime but due to time constraints
 * this aspect of the code was not made modular.
 */

double xlim1 = 0.0;                       //Lower limit on x-axis
double xlim2 = 1.0;                      //Upper limit on x-axis

// for simplicity, nx=2*multiple
unsigned int nx = 128;                         //Number of nodes in x-direction INCLUDING THE EXTREME VALUES

double ylim1 = 0.0;                       //Lower limit on y-axis
double ylim2 = 1.0;                     //Upper limit on y-axis

// for simplicity, ny=2*multiple
unsigned int ny = 128;                        //Number of nodes INCLUDING THE EXTREME VALUES

double dt = 0.5 * 1.0/128.0;                     //Length of time step
double Tfinal = 1.0;                    //Total time period for the simulation
double T_period = 1.0;                  //Period of the velocity field

unsigned int option = 1;                         //Option - if you need animation initialize at 1 else initialize at 2
unsigned int printstep = 16;                      //How frequently do you want to store the images (every nth time step)

char psischeme[] = "SuperConsistent";   //'SuperConsistent' or 'Heuns'
char backtrace_scheme[] = "RK3" ;      //'Euler' or 'RK3'
