//
//  GALS - Before compiling the program, update the section of this program that is in the beginning of main and update all the functions in Initializer.h
//
//
//  Created by Raunak Bardia on 10/22/14.
//
//

// DISCLAIMER:
// Use the indexes carefully
// First index of array represents movement along y-direction because it represents rows
// Second index of array represents movement along x-direction because it represents columns
//
// Implementing GALS for a given initial level set function
// in a specified velocity field for a grid of cells
//
// Given -
// 1. Defining function at t = 0 which implies that phi and psi values are available for all node points at t = 0
// 2. Given velocity for the complete domain at all times
//
// All required data is stored in separate 2D matrices of phi, psix, psiy and psixy
// Boundary Condition grad(velocity).n > 0
//

// THIS IMPLEMENTATION WON'T WORK IF THE GRID IS SMALLER THAN (2 X 2)

#include <iostream>
#include <iomanip>
#include <math.h>
#include <stdio.h>
#include <fstream>
#include <time.h>
#include <string.h>
//#include "Bilinear.h"
#include "Hermite.h"
#include "Initializer_vortex.h"    //UPDATE THIS HEADER FILE WITH THE RELEVANT FUNCTIONS

#include <chrono>
using namespace std;

// y direction is the first index of array, x direction is the second index of array

int main(){
    
    clock_t tStart = clock();
    
    /* UPDATE ALL THE FOLLOWING VALUES */
    double xlim1 = 0.0;                       //Lower limit on x-axis
    double xlim2 = 	1.0;                       //Upper limit on x-axis
    int nx = 32 + 1;                        //Number of nodes in x-direction INCLUDING THE EXTREME VALUES
    
    double ylim1 = 0.0;                       //Lower limit on y-axis
    double ylim2 = 1.0;                     //Upper limit on y-axis
    int ny = 32 + 1;                        //Number of nodes INCLUDING THE EXTREME VALUES
    
    double dt = (1/33.0);                     //Length of time step
    double Tfinal = 1.0;                    //Total time period for the simulation
    int option = 1;                         //Option - if you need animation initialize at 1 else initialize at 2
    int printstep = 8;                      //How frequently do you want to store the images (every nth time step)
    char psischeme[] = "SuperConsistent";   //'SuperConsistent' or 'Heuns'
    char backtrace_scheme[] = "RK3" ;      //'Euler' or 'RK3'
    double T_period = 1.0;                  //Period of the velocity field
    
    //MAKE SURE THAT YOU HAVE ENOUGH MEMORY SPACE IF YOU ARE STORING A LOT OF TIME STEP VALUES BECAUSE IT STORES ACROSS GRID POINTS FOR EACH PRINTSTEP
    
    /* USER UPDATE OVER */
    
    int n = Tfinal/dt; //Number of time steps
    // If only the initial and final profiles are needed
    if(option != 1)
        printstep = n;
    
    
    /*
     
     Let the following represent one cell
     
     2      3
     *      *
     
     
     *      *
     0      1
     
     the value of the loop for this cell varies from 0 -> 3 but the (i, j) coordinates that represent these points in an array change as (i,j), (i,j+1), (i+1,j), (i+1, j+1)
     
     Hence, tempindexes take care of these changes
     
     */
    
    double xo, yo, rootpsix, rootpsiy; // {xo,yo} - first node in each cell; rootpsi* - temporary variables that stores updated psi*
    
    // Temporary variables to store the values of phi, psix, psiy and psixy at all the nodes of the cell being considered for the iteration
    double phi[4];
    double psix[4];
    double psiy[4];
    double psixy[4];
    int tempindex_x, tempindex_y;   // These variable select each node one by one in a cell
    int cornerx, cornery;   //To store indexes of corners
    int corner; // To save which node of the cell is the grid corner if a corner is encountered
    // Float constants for RK-3
    
    double c1 = (1/6.0);
    double c2 = (1/6.0);
    double c3 = (2/3.0);    //RK-3 constants
    
    // Defining node points
    double dx = (xlim2 - xlim1)/(nx - 1);
    double x[nx];
    for(int i = 0; i < nx; i++){
        x[i] = xlim1 + dx * i;
    }
    double dy = (ylim2 - ylim1)/(ny - 1);
    double y[ny];
    for(int i = 0; i < ny; i++){
        y[i] = ylim1 + dy * i;
    }
    // Node point definition ends
    
    int temp1, temp2;       // These temporary variables are updated when the advection point lies between the x and y limits of the cell
    int tracker[ny][nx];    // Tracker, tracks every node point and updates that particular index to 1 if that point has been updated
    
    // master matrices of phi, psix, psiy, psixy which are updated at each time step & temporary matrices which is updated constantly throughout the cell iterations
    double mphi[ny][nx], mpsix[ny][nx], mpsiy[ny][nx], mpsixy[ny][nx], tempphi[ny][nx], temppsix[ny][nx], temppsiy[ny][nx], temppsixy[ny][nx];
    
    double u[ny][nx], v[ny][nx];    //Storing the velocity at node points only for those time steps when printstep is initiated
    double d1, d2, d3, d4;  //Temporary variable to compute mixed derivatives
    double xadv, yadv;  //Advection point
    double np = nx * ny; // Number of node points
    
    // ux, vy - X velocity & Y velocity; ux1, vy1, ux2, vy2 - X & Y velocity of the xadv & yadv at separate RK steps
    double ux, vy, ux1, vy1, ux2, vy2;
    
    //FOR SUPERCONSISTENT METHOD & HEUNS
    double xr1,yr1,xr2,yr2,gradux1,graduy1,gradvx1,gradvy1,gradux2,graduy2,gradvx2,gradvy2,graduxnode,graduynode,gradvxnode,gradvynode,graduxroot,graduyroot,gradvxroot,gradvyroot;
    double gradxr1x,gradyr1x,gradxr2x,gradyr2x,gradxr1y,gradyr1y,gradxr2y,gradyr2y,gradxrx,gradxry,gradyrx,gradyry, graduxnew, graduynew, gradvxnew, gradvynew, temp;
    
    double t; // Time variable
    
   /* 
    ofstream myfile, myfile1, myfile2, myfile3, details, uu, vv;
    // Removing existing files with these names if any
    int i = remove("phi.txt");
    i = remove("psix.txt");
    i = remove("psiy.txt");
    i = remove("psixy.txt");
    i = remove("details.txt");
    i = remove("Velocity_x.txt");
    i = remove("Velocity_y.txt");
    
    // Opening new files
    myfile.open("phi.txt", ios::out);
    myfile1.open("psix.txt", ios::out);
    myfile2.open("psiy.txt", ios::out);
    myfile3.open("psixy.txt", ios::out);
    details.open("details.txt", ios::out);
    uu.open("Velocity_x.txt", ios::out);
    vv.open("Velocity_y.txt",ios::out);
    */
    // Initializing at t = 0
    for(int i = 0; i < ny; i++){
        
        // Initializing the first element of each row separately to avoid trailing commas in the file
        // Storing these values in the master variable as well as temporary variable matrices
        
        mphi[i][0] = initialize(x[0], y[i]);
        tempphi[i][0] = mphi[i][0];
    //    myfile << std::fixed << std::setprecision(10) << mphi[i][0];
        mpsix[i][0] = derivxinit(x[0], y[i]);
        temppsix[i][0] = mpsix[i][0];
      //  myfile1 << std::fixed << std::setprecision(10) << mpsix[i][0];
        mpsiy[i][0] = derivyinit(x[0], y[i]);
        temppsiy[i][0] = mpsiy[i][0];
        //myfile2 << std::fixed << std::setprecision(10) << mpsiy[i][0];
        mpsixy[i][0] = derivxyinit(x[0], y[i]);
        temppsixy[i][0] = mpsixy[i][0];
      //  myfile3 << std::fixed << std::setprecision(10) << mpsixy[i][0];
        
        u[i][0] = Velx(x[0],y[i],0, T_period);
      //  uu << std::fixed << std::setprecision(10) << u[i][0];
        v[i][0] = Vely(x[0],y[i],0, T_period);
      //  vv << std::fixed << std::setprecision(10) << v[i][0];
        
        for(int j = 1; j < nx; j++){
            mphi[i][j] = initialize(x[j], y[i]);
            tempphi[i][j] = mphi[i][j];
  //          myfile << ",";
   //         myfile << std::fixed << std::setprecision(10) << mphi[i][j];
            mpsix[i][j] = derivxinit(x[j], y[i]);
            temppsix[i][j] = mpsix[i][j];
  //          myfile1 << ",";
  //          myfile1 << std::fixed << std::setprecision(10) << mpsix[i][j];
            mpsiy[i][j] = derivyinit(x[j], y[i]);
            temppsiy[i][j] = mpsiy[i][j];
  //          myfile2 << ",";
  //          myfile2 << std::fixed << std::setprecision(10) << mpsiy[i][j];
            mpsixy[i][j] = derivxyinit(x[j], y[i]);
            temppsixy[i][j] = mpsixy[i][j];
  //          myfile3 << ",";
   //         myfile3 << std::fixed << std::setprecision(10) << mpsixy[i][j];
            
            u[i][j] = Velx(x[j],y[i],0, T_period);
    //        uu << "," << std::fixed << std::setprecision(10) << u[i][j];
            v[i][j] = Vely(x[j],y[i],0, T_period);
      //      vv << "," << std::fixed << std::setprecision(10) << v[i][j];
        }
      /*  myfile << '\n';
        myfile1 << '\n';
        myfile2 << '\n';
        myfile3 << '\n';
        uu << '\n';
        vv << '\n';*/
    }
    // Giving an extra empty line for distinction between next time step
   /* myfile << '\n';
    myfile1 << '\n';
    myfile2 << '\n';
    myfile3 << '\n';
    uu << '\n';
    vv << '\n';
    
    details<< nx << "," << ny << "," << std::fixed << std::setprecision(10) << dx << "," << dy << "," << xlim1 << "," << xlim2 << "," << ylim1 << "," << ylim2 << "," << n << "," << dt << "," << printstep;
    details.close();*/
    float t_calcpts=0.0,t_update=0.0,t_mixed=0.0,t_copy=0.0;
    auto tbegin = chrono::high_resolution_clock::now();
    for(t = 0; t < n; t++){
        
        // Initializing tracker to zero at every time step
        for(int p = 0; p < ny; p ++){
            for(int q = 0; q<nx; q++){
                tracker[p][q] = 0;
            }
        }
        
        for(int i = 0; i < ny-1; i++){  // loop for y - rows in the array
            for(int j = 0; j < nx-1; j++){  // loop for x - columns in the array
                
                // Storing the four values for four nodes of each cell
                phi[0] = mphi[i][j];
                psix[0] = mpsix[i][j];
                psiy[0] = mpsiy[i][j];
                psixy[0] = mpsixy[i][j];
                phi[1] = mphi[i][j+1];
                psix[1] = mpsix[i][j+1];
                psiy[1] = mpsiy[i][j+1];
                psixy[1] = mpsixy[i][j+1];
                phi[2] = mphi[i+1][j];
                psix[2] = mpsix[i+1][j];
                psiy[2] = mpsiy[i+1][j];
                psixy[2] = mpsixy[i+1][j];
                phi[3] = mphi[i+1][j+1];
                psix[3] = mpsix[i+1][j+1];
                psiy[3] = mpsiy[i+1][j+1];
                psixy[3] = mpsixy[i+1][j+1];
                // Node value assignment ends
                
                // Storing the coordinates of first node of the working cell
                xo = x[j];
                yo = y[i];
                corner = -1;
                
    	auto tbegin_function = chrono::high_resolution_clock::now();
                for(int k = 0; k < 4; k++){ //4 nodes on the cell
                    
                    int count = 0;  // Additional variable to make sure that the advection point is available for the relevant node point
                    temp1 = 0;
                    temp2 = 0;
                    tempindex_y = i + (-2 * pow(k,3) + 9 * pow(k,2) - 7 * k)/6;   // Temporary Indexes polynomial
                    tempindex_x = j + (2 * pow(k,3) - 9 * pow(k,2) + 10 * k)/3;
                    
                    if(tracker[tempindex_y][tempindex_x] == 0){ // If the node has already been updated we need not update it again
                        
                        ux = Velx(x[tempindex_x], y[tempindex_y], (t + 1) * dt, T_period);
                        vy = Vely(x[tempindex_x], y[tempindex_y], (t + 1) * dt, T_period);
                        
                        // Advected Points - 1 Step
                        xr1 = x[tempindex_x] - ux * dt;
                        yr1 = y[tempindex_y] - vy * dt;

                        if(strcmp("RK3",backtrace_scheme) == 0){
                            
                            ux1 = Velx(xr1, yr1, t * dt, T_period);
                            vy1 = Vely(xr1, yr1, t * dt, T_period);
                            
                            // 2 Step
                            xr2 = x[tempindex_x] - ((dt)/2.0 * (0.5 * ux + 0.5 * ux1));
                            yr2 = y[tempindex_y] - ((dt)/2.0 * (0.5 * vy + 0.5 * vy1));
                            
                            ux2 = Velx(xr2, yr2,(t * dt + dt/2.0), T_period);
                            vy2 = Vely(xr2, yr2,(t * dt + dt/2.0), T_period);
                            
                            // 3 Step
                            xadv = x[tempindex_x] - (dt * (c1 * ux + c2 * ux1 + c3 * ux2));
                            yadv = y[tempindex_y] - (dt * (c1 * vy + c2 * vy1 + c3 * vy2));
                        }
                        else{
                            xadv = xr1;
                            yadv = yr1;
                        }
    	auto tend_function = chrono::high_resolution_clock::now();
        t_calcpts += chrono::duration_cast<chrono::nanoseconds>(tend_function-tbegin_function).count();
                        
                        if((xadv >= x[j]) && (xadv <= x[j+1]))
                            temp1 = 1;
                        if((yadv >= y[i]) && (yadv <= y[i+1]))
                            temp2 = 1;
                        
                        //count variable makes sure that the values are updated if they can be through this cell's hermite interpolant itself
                        
                        // If the advection point lies within the cell then we don't need to change the advection point
                        if(temp1 && temp2){ count = 1; }
                        
                        // If only the x coordinate is within the cell limits then only update the point if it is on y-boundary
                        if(temp1 && !temp2){
                            if((yadv < ylim1) || (yadv > ylim2)){
                                count = 2;
                            }
                        }
                        
                        // If only y coordinate is within the cell limits then only update the point if it is on x-boundary
                        if(!temp1 && temp2){
                            if((xadv < xlim1) || (xadv > xlim2)){
                                count = 3;
                            }
                        }
                        
                        // If both coordinates dont lie in the cell limits, then the value is not changed for phi and first derivatives
                        if(!temp1 && !temp2){
                            if(((tempindex_y == 0) || (tempindex_y == ny-1)) && ((tempindex_x == 0) || (tempindex_x == nx-1))){
                                count = 4;
                            }
                        }
                        if(((tempindex_y == 0) || (tempindex_y==ny-1)) && ((tempindex_x == 0) || (tempindex_x==nx-1))){
                            corner = k; // This variable is to check if we have a grid corner in the cell because the psixy value for corners is calculated differently
                        }
                        
                        //This condition and variable has been defined so that these statements need to be written only once otherwise all these statements would have been repeated in all the above if's
    	tbegin_function = chrono::high_resolution_clock::now();
                        if(count != 0)
                        {
                            switch(count){
                                case 1:{
                                    tempphi[tempindex_y][tempindex_x] = hp(phi, psix, psiy, psixy, xadv, yadv, xo, yo, dx, dy);
                                    rootpsix = hermx(phi,psix,psiy,psixy,xadv,yadv,xo,yo,dx, dy);
                                    rootpsiy = hermy(phi,psix,psiy,psixy,xadv,yadv,xo,yo,dx, dy);
                                    
                                    if(strcmp("Heuns",psischeme) == 0){
                                        
                                        // The following derivatives are the velocity gradients required for the third equation in the mehtod of characteristics for GRADIENTS
                                        graduxnode = gradUx(x[tempindex_x], y[tempindex_y], (t + 1) * dt, T_period);
                                        graduynode = gradUy(x[tempindex_x], y[tempindex_y], (t + 1) * dt, T_period);
                                        gradvxnode = gradVx(x[tempindex_x], y[tempindex_y], (t + 1) * dt, T_period);
                                        gradvynode = gradVy(x[tempindex_x], y[tempindex_y], (t + 1) * dt, T_period);
                                        
                                        graduxroot = gradUx(xadv, yadv, t * dt, T_period);
                                        graduyroot = gradUy(xadv, yadv, t * dt, T_period);
                                        gradvxroot = gradVx(xadv, yadv, t * dt, T_period);
                                        gradvyroot = gradVy(xadv, yadv, t * dt, T_period);
                                        
                                        // Heun's RK-2 method - first step - for x & y
                                        temppsix[tempindex_y][tempindex_x] = rootpsix - dt * (graduxroot * rootpsix + gradvxroot * rootpsiy);
                                        temppsiy[tempindex_y][tempindex_x] = rootpsiy - dt * (graduyroot * rootpsix + gradvyroot * rootpsiy);
                                        
                                        // Heun's RK-2 method - second step - for x & y
                                        temppsix[tempindex_y][tempindex_x] = rootpsix - 0.5 * dt * (graduxroot * rootpsix + gradvxroot * rootpsiy + graduxnode * temppsix[tempindex_y][tempindex_x] + gradvxnode * temppsiy[tempindex_y][tempindex_x]);
                                        temppsiy[tempindex_y][tempindex_x] = rootpsiy - 0.5 * dt * (graduyroot * rootpsix + gradvyroot * rootpsiy + graduynode * temppsix[tempindex_y][tempindex_x] + gradvynode * temppsiy[tempindex_y][tempindex_x]);
                                    }
                                    
                                    else if(strcmp("SuperConsistent",psischeme) == 0){
                                        
                                        graduxnode = gradUx(x[tempindex_x], y[tempindex_y], (t + 1) * dt, T_period);
                                        graduynode = gradUy(x[tempindex_x], y[tempindex_y], (t + 1) * dt, T_period);
                                        gradvxnode = gradVx(x[tempindex_x], y[tempindex_y], (t + 1) * dt, T_period);
                                        gradvynode = gradVy(x[tempindex_x], y[tempindex_y], (t + 1) * dt, T_period);
                                        
                                        gradxr1x = 1 - dt * graduxnode;
                                        gradxr1y = 0 - dt * graduynode;
                                        gradyr1x = 0 - dt * gradvxnode;
                                        gradyr1y = 1 - dt * gradvynode;
                                        
                                        if(strcmp("RK3",backtrace_scheme) == 0){
                                            
                                            gradux1 = gradUx(xr1, yr1, t * dt, T_period);
                                            graduy1 = gradUy(xr1, yr1, t * dt, T_period);
                                            gradvx1 = gradVx(xr1, yr1, t * dt, T_period);
                                            gradvy1 = gradVy(xr1, yr1, t * dt, T_period);
                                            
                                            gradux2 = gradUx(xr2, yr2, (t * dt + dt/2.0), T_period);
                                            graduy2 = gradUy(xr2, yr2, (t * dt + dt/2.0), T_period);
                                            gradvx2 = gradVx(xr2, yr2, (t * dt + dt/2.0), T_period);
                                            gradvy2 = gradVy(xr2, yr2, (t * dt + dt/2.0), T_period);
                                            
                                            gradxr2x = 1 - dt/4.0 * (graduxnode + gradxr1x * gradux1 + gradyr1x * graduy1);
                                            gradxr2y = 0 - dt/4.0 * (graduynode + gradxr1y * gradux1 + gradyr1y * graduy1);
                                            gradyr2x = 0 - dt/4.0 * (gradvxnode + gradxr1x * gradvx1 + gradyr1x * gradvy1);
                                            gradyr2y = 1 - dt/4.0 * (gradvynode + gradxr1y * gradvx1 + gradyr1y * gradvy1);
                                            
                                            gradxrx = 1 - dt/6.0 * (graduxnode + gradxr1x * gradux1 + gradyr1x * graduy1 + 4 * (gradxr2x * gradux2 + gradyr2x * graduy2));
                                            gradxry = 0 - dt/6.0 * (graduynode + gradxr1y * gradux1 + gradyr1y * graduy1 + 4 * (gradxr2y * gradux2 + gradyr2y * graduy2));
                                            gradyrx = 0 - dt/6.0 * (gradvxnode + gradxr1x * gradvx1 + gradyr1x * gradvy1 + 4 * (gradxr2x * gradvx2 + gradyr2x * gradvy2));
                                            gradyry = 1 - dt/6.0 * (gradvynode + gradxr1y * gradvx1 + gradyr1y * gradvy1 + 4 * (gradxr2y * gradvx2 + gradyr2y * gradvy2));
                                        }
                                        else{
                                            gradxrx = gradxr1x;
                                            gradxry = gradxr1y;
                                            gradyrx = gradyr1x;
                                            gradyry = gradyr1y;
                                        }
                                        
                                        
                                        // Super consistent method -> psi(n+1) = grad(xroot,yroot) * psiroot(n)
                                        temppsix[tempindex_y][tempindex_x] = (gradxrx * rootpsix + gradyrx * rootpsiy);
                                        temppsiy[tempindex_y][tempindex_x] = (gradxry * rootpsix + gradyry * rootpsiy);
                                    }
                                    tracker[tempindex_y][tempindex_x] = 1;
                                    break;
                                }
                                case 2:{
                                    rootpsix = hermx(phi,psix,psiy,psixy,xadv, y[tempindex_y],xo,yo,dx, dy);
                                    rootpsiy = temppsiy[tempindex_y][tempindex_x];
                                    tempphi[tempindex_y][tempindex_x] = hp(phi, psix, psiy, psixy, xadv, y[tempindex_y], xo, yo, dx, dy) - dt * Vely(x[tempindex_x], y[tempindex_y], t * dt, T_period) * rootpsiy;
                                    // The following derivatives are the velocity gradients required for the third equation in the mehtod of characteristics for GRADIENTS
                                    
                                    graduxnode = gradUx(x[tempindex_x], y[tempindex_y], t * dt, T_period);
                                    graduynode = gradUy(x[tempindex_x], y[tempindex_y], t * dt, T_period);
                                    gradvxnode = gradVx(x[tempindex_x], y[tempindex_y], t * dt, T_period);
                                    gradvynode = gradVy(x[tempindex_x], y[tempindex_y], t * dt, T_period);
                                    
                                    graduxnew = gradUx(x[tempindex_x], y[tempindex_y], (t+1) * dt, T_period);
                                    graduynew = gradUy(x[tempindex_x], y[tempindex_y], (t+1) * dt, T_period);
                                    gradvxnew = gradVx(x[tempindex_x], y[tempindex_y], (t+1) * dt, T_period);
                                    gradvynew = gradVy(x[tempindex_x], y[tempindex_y], (t+1) * dt, T_period);
                                    
                                    temp = rootpsix - dt * (gradvxnode * rootpsiy + Vely(x[tempindex_x], y[tempindex_y], t * dt, T_period) * temppsixy[tempindex_y][tempindex_x] + graduxnode * temppsix[tempindex_y][tempindex_x]);
                                    
                                    //HEUNS
                                    temppsix[tempindex_y][tempindex_x] = rootpsix - (dt/2.0) * (gradvxnode * rootpsiy + Vely(x[tempindex_x], y[tempindex_y], t * dt, T_period) * temppsixy[tempindex_y][tempindex_x] + graduxnode * temppsix[tempindex_y][tempindex_x] + gradvxnew * rootpsiy + Vely(x[tempindex_x], y[tempindex_y], (t+1) * dt, T_period) * temppsixy[tempindex_y][tempindex_x] + graduxnew * temp);
                                    tracker[tempindex_y][tempindex_x] = 2;
                                    break;
                                }
                                case 3:{
                                    
                                    rootpsix = temppsix[tempindex_y][tempindex_x];
                                    rootpsiy = hermy(phi,psix,psiy,psixy,x[tempindex_x],yadv,xo,yo,dx, dy);
                                    tempphi[tempindex_y][tempindex_x] = hp(phi, psix, psiy, psixy, x[tempindex_x],yadv, xo, yo, dx, dy) - dt * Velx(x[tempindex_x], y[tempindex_y], t * dt, T_period) * rootpsix;
                                    // The following derivatives are the velocity gradients required for the third equation in the mehtod of characteristics for GRADIENTS
                                    
                                    graduxnode = gradUx(x[tempindex_x], y[tempindex_y], t * dt, T_period);
                                    graduynode = gradUy(x[tempindex_x], y[tempindex_y], t * dt, T_period);
                                    gradvxnode = gradVx(x[tempindex_x], y[tempindex_y], t * dt, T_period);
                                    gradvynode = gradVy(x[tempindex_x], y[tempindex_y], t * dt, T_period);
                                    
                                    graduxnew = gradUx(x[tempindex_x], y[tempindex_y], (t+1) * dt, T_period);
                                    graduynew = gradUy(x[tempindex_x], y[tempindex_y], (t+1) * dt, T_period);
                                    gradvxnew = gradVx(x[tempindex_x], y[tempindex_y], (t+1) * dt, T_period);
                                    gradvynew = gradVy(x[tempindex_x], y[tempindex_y], (t+1) * dt, T_period);
                                    
                                    temp = rootpsiy - dt * (graduynode * rootpsix + Velx(x[tempindex_x], y[tempindex_y], t * dt, T_period) * temppsixy[tempindex_y][tempindex_x] + gradvynode * temppsiy[tempindex_y][tempindex_x]);
                                    
                                    //HEUNS
                                    temppsiy[tempindex_y][tempindex_x] = rootpsiy - (dt/2.0) * (graduynode * rootpsix + Velx(x[tempindex_x], y[tempindex_y], t * dt, T_period) * temppsixy[tempindex_y][tempindex_x] + gradvynode * temppsiy[tempindex_y][tempindex_x] + graduynew * rootpsix + Velx(x[tempindex_x], y[tempindex_y], (t+1) * dt, T_period) * temppsixy[tempindex_y][tempindex_x] + gradvynew * temp);
                                    tracker[tempindex_y][tempindex_x] = 3;
                                    break;
                                }
                                case 4:{
                                    tempphi[tempindex_y][tempindex_x] = tempphi[tempindex_y][tempindex_x] - dt * (Velx(x[tempindex_x], y[tempindex_y], t * dt, T_period) * temppsix[tempindex_y][tempindex_x] + Vely(x[tempindex_x], y[tempindex_y], t * dt, T_period) * temppsiy[tempindex_y][tempindex_x]);
                                    tracker[tempindex_y][tempindex_x] = 4;
                                    break;
                                }
                                default:{
                                    break;
                                }
                            }
                            count = 0;
                        }
    	tend_function = chrono::high_resolution_clock::now();
        t_update += chrono::duration_cast<chrono::nanoseconds>(tend_function-tbegin_function).count();
                    }
                }
            }
        }
        
        // Update the mixed derivatives now for the remaining grid points
        double dxy1, dxy2;;
        
    	auto tbegin_function = chrono::high_resolution_clock::now();
        for(int k = 0; k < ny; k++){
            for(int l = 0; l < nx; l++){
                if ((k == 0 || k == ny - 1) && (l != 0 && l != nx - 1)){
                    temppsixy[k][l] = (temppsiy[k][l+1] - temppsiy[k][l-1])/(2 * dx);
                }
                else
                    if ((k != 0 && k != ny - 1) && (l == 0 || l == nx - 1)){
                        temppsixy[k][l] = (temppsix[k+1][l] - temppsix[k-1][l])/(2 * dy);
                    }
                    else
                        if((k == 0 || k == ny - 1) && (l == 0 || l == nx - 1)){
                            if(k == 0 && l == 0){
                                d1 = (temppsiy[0][1] - temppsiy[0][0])/dx;
                                d2 = (temppsix[1][0] - temppsix[0][0])/dy;
                                d3 = (temppsix[1][1] - temppsix[0][1])/dy;
                                d4 = (temppsiy[1][1] - temppsiy[1][0])/dx;
                                temppsixy[k][l] = 0.75 * (d1 + d2) - 0.25 * (d3 + d4);
                            }
                            else if(k == 0 && l == nx-1){
                                d1 = (temppsiy[0][nx-1] - temppsiy[0][nx-2])/dx;
                                d2 = (temppsix[1][nx-2] - temppsix[0][nx-2])/dy;
                                d3 = (temppsix[1][nx-1] - temppsix[0][nx-1])/dy;
                                d4 = (temppsiy[1][nx-1] - temppsiy[1][nx-2])/dx;
                                temppsixy[k][l] = 0.75 * (d1 + d3) - 0.25 * (d2 + d4);
                                
                            }
                            else if(k == ny-1 && l == 0){
                                d1 = (temppsiy[ny-2][1] - temppsiy[ny-2][0])/dx;
                                d2 = (temppsix[ny-1][0] - temppsix[ny-2][0])/dy;
                                d3 = (temppsix[ny-1][1] - temppsix[ny-2][1])/dy;
                                d4 = (temppsiy[ny-1][1] - temppsiy[ny-1][0])/dx;
                                temppsixy[k][l] = 0.75 * (d2 + d4) - 0.25 * (d3 + d1);
                                
                            }
                            else if(k == ny-1 && l == nx-1){
                                d1 = (temppsiy[ny-2][nx-1] - temppsiy[ny-2][nx-2])/dx;
                                d2 = (temppsix[ny-1][nx-2] - temppsix[ny-2][nx-2])/dy;
                                d3 = (temppsix[ny-1][nx-1] - temppsix[ny-2][nx-1])/dy;
                                d4 = (temppsiy[ny-1][nx-1] - temppsiy[ny-1][nx-2])/dx;
                                temppsixy[k][l] = 0.75 * (d3 + d4) - 0.25 * (d1 + d2);
                            }
                            
                        }
                        else{
                            dxy1 = (temppsiy[k][l+1] - temppsiy[k][l-1])/(2 * dx);
                            dxy2 = (temppsix[k+1][l] - temppsix[k-1][l])/(2 * dy);
                            temppsixy[k][l] = (dxy1 + dxy2)/2.0;
                        }
            }
        }
    	auto tend_function = chrono::high_resolution_clock::now();
        t_mixed += chrono::duration_cast<chrono::nanoseconds>(tend_function-tbegin_function).count();
        
        // Feeding values back to the master matrix
        
    	tbegin_function = chrono::high_resolution_clock::now();
        for(int k = 0; k < ny; k++){
            for(int l = 0; l < nx; l++){
                mphi[k][l] = tempphi[k][l];
                mpsix[k][l] = temppsix[k][l];
                mpsiy[k][l] = temppsiy[k][l];
                mpsixy[k][l] = temppsixy[k][l];
            }
        }
    	tend_function = chrono::high_resolution_clock::now();
        t_copy += chrono::duration_cast<chrono::nanoseconds>(tend_function-tbegin_function).count();
       /* 
        // Feeding phi, psix, psiy and psixy values in their respective files
        int check = t+1;
        if(check % printstep == 0){
            for(int i = 0; i < ny; i++){
                
                myfile << std::fixed << std::setprecision(10) <<  mphi[i][0];
                myfile1 << std::fixed << std::setprecision(10) << mpsix[i][0];
                myfile2 << std::fixed << std::setprecision(10) << mpsiy[i][0];
                myfile3 << std::fixed << std::setprecision(10) << mpsixy[i][0];
                
                u[i][0] = Velx(x[0],y[i],(t+1)*dt, T_period);
                uu << std::fixed << std::setprecision(10) << u[i][0];
                v[i][0] = Vely(x[0],y[i],(t+1)*dt, T_period);
                vv << std::fixed << std::setprecision(10) << v[i][0];
                
                for(int j = 1; j < nx; j++){
                    myfile << ",";
                    myfile << std::fixed << std::setprecision(10) << mphi[i][j];
                    myfile1 << ",";
                    myfile1 << std::fixed << std::setprecision(10) << mpsix[i][j];
                    myfile2 << ",";
                    myfile2 << std::fixed << std::setprecision(10) << mpsiy[i][j];
                    myfile3 << ",";
                    myfile3 << std::fixed << std::setprecision(10) << mpsixy[i][j];
                    
                    u[i][j] = Velx(x[j],y[i],(t+1)*dt, T_period);
                    uu << "," << std::fixed << std::setprecision(10) << u[i][j];
                    v[i][j] = Vely(x[j],y[i],(t+1)*dt, T_period);
                    vv << "," << std::fixed << std::setprecision(10) << v[i][j];
                }
                myfile << '\n';
                myfile1 << '\n';
                myfile2 << '\n';
                myfile3 << '\n';
                uu << '\n';
                vv << '\n';
            }
            myfile << '\n';
            myfile1 << '\n';
            myfile2 << '\n';
            myfile3 << '\n';
            uu << '\n';
            vv << '\n';
        }
        
        cout<< t+1;
        cout<< " Time Step Completed" <<'\n';
        */
    }
    
/*    myfile.close();
    myfile1.close();
    myfile2.close();
    myfile3.close();
    uu.close();
    vv.close();*/
    auto tend = chrono::high_resolution_clock::now();
    float duration = chrono::duration_cast<chrono::nanoseconds>(tend-tbegin).count();
    duration = duration * pow(10.0,-6);
    t_calcpts = t_calcpts * pow(10.0,-6);
    t_update = t_update * pow(10.0,-6);
    t_mixed = t_mixed * pow(10.0,-6);
    t_copy = t_copy * pow(10.0,-6);
    cout << "Time taken for calculation of advection points = " << t_calcpts << '\n';
    cout << "Time taken for hermite Update = " << t_update << '\n';
    cout << "Time taken for copying the matrix = " << t_copy << '\n';
    cout << "Time taken for calculation of mixed derivatives = " << t_mixed << '\n';
    cout << "Total Duration of the Time Loop = " << duration << endl;
    return 0;
}
