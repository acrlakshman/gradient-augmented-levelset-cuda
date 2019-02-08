//  GALS - Before compiling the program, update the section of this program that is in the beginning of main and update all the initializer files
//
//  Created by Raunak Bardia, Chia-Wei Kuo, and Arpit Agarwal on December 20, 2017.
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

/* * * * * * * * * * * * * *  MPI IMPLEMENTATION * * * * * * * * * * * * * * */

// THIS IMPLEMENTATION WON'T WORK IF THE GRID IS SMALLER THAN (2 X 2)
#include <iostream>
#include <iomanip>
#include <math.h>
#include <stdio.h>
#include <fstream>
#include <sys/time.h>
#include <time.h>
#include <string.h>
#include <vector>
#include <tuple>
#include "VariableDefinitions.h"
#include "Hermite.h"
#include "InitializeLevelSet.h"    //UPDATE THIS HEADER FILE WITH THE RELEVANT FUNCTIONS
#include "VortexVelocity.h"    //UPDATE THIS HEADER FILE WITH THE RELEVANT FUNCTIONS
#include "TimeSteppingMethods.h"
#include "Allocation.h"
#include "AdvectionPointCalcs.h"
#include "Constants.h"

#include <mpi.h>  	 // used for MPI
#include "defineMPI.h"   // used for MPI

using namespace std;
using namespace galsfunctions;
using namespace mpi;

// y direction is the first index of array, x direction is the second index of array

int main(int argc, char **argv){
    
    /* UPDATE ALL THE FOLLOWING VALUES */
#include "Setting.h"
    //MAKE SURE THAT YOU HAVE ENOUGH MEMORY SPACE IF YOU ARE STORING A LOT OF TIME STEP VALUES BECAUSE IT STORES ACROSS GRID POINTS FOR EACH PRINTSTEP
    /* USER UPDATE OVER */
    
    struct  timeval start1;       struct  timeval end1;          unsigned  long diff1;
    
    
    /***** used for MPI only ******/
    //
    context ctx(&argc, &argv);
    //
    /***** used for MPI only ******/
    
    
    unsigned int n = Tfinal/dt; //Number of time steps
    if(option != 1)
        printstep = n;
    
    // Node Locations
    double dx = (xlim2 - xlim1)/(nx - 1);
    double dy = (ylim2 - ylim1)/(ny - 1);
    vectorarray x,y;
    gridnodes(x,y,xlim1,ylim1,dx,dy,nx,ny);
    
    // Initializing at t = 0
    gridarray mphi, mpsix, mpsiy, mpsixy; // level set matrix
    allocate_levelset_matrices(mphi,mpsix,mpsiy,mpsixy,x,y,nx,ny); //Initializing level set matrices
    
    gridarray tempphi(mphi), temppsix(mpsix), temppsiy(mpsiy), temppsixy(mpsixy);	// Making level set matrix copies for time loop use
    
    
    // Removing existing files with these names if any
    remove("phi.txt");
    remove("psix.txt");
    remove("psiy.txt");
    remove("psixy.txt");
    remove("details.txt");
    remove("Velocity_x.txt");
    remove("Velocity_y.txt");
    fileprint(mphi,mpsix,mpsiy,mpsixy,nx,ny,x,y,0.0,T_period);
    
    /*
     * Let the following represent one cell
     *
     * 2      3
     *      *
     *
     *
     *      *
     * 0      1
     *
     * the value of the loop for this cell varies from 0 -> 3 but the (i, j) coordinates that represent these points in an array change as (i,j), (i,j+1), (i+1,j), (i+1, j+1)
     *
     * Hence, tempindexes take care of these changes
     *
     */
    
    ofstream details;
    details.open("details.txt", ios::out | ios::app);
    details<< nx << "," << ny << "," << std::fixed << std::setprecision(10) << dx << "," << dy << "," << xlim1 << "," << xlim2 << "," << ylim1 << "," << ylim2 << "," << n << "," << dt << "," << printstep;
    details.close();
    
    // TIME STEPPING LOOP
    // If only the initial and final profiles are needed
    gettimeofday(&start1,NULL);  // start the timer
    for(unsigned int t = 0; t < n; t++){
        
        intgridarray tracker;
        tracker.resize(ny,intvectorarray(nx,0));
        
        gridarray xadv, yadv;
        xadv.resize(ny,vectorarray(nx,0.0));
        yadv.resize(ny,vectorarray(nx,0.0));
        
        
        // used for MPI only, here I use 2 nodes
        unsigned int n_Nodes = 2;
        
        // initialize the memory contiguous matrices
        double **A1, **B1, **A3, **B3, **C3;
        unsigned int nnx = nx/n_Nodes;
        
        A1 = alloc_2d_int(nnx,ny);
        B1 = alloc_2d_int(nnx,ny);
        A3 = alloc_2d_int(nnx,ny);
        B3 = alloc_2d_int(nnx,ny);
        C3 = alloc_2d_int(nnx,ny);
        // end of initialization
        
// 	advection_point_MPI(x, y, xadv, yadv, t, dt, T_period,0,nx, backtrace_scheme);
        
        // node-1
        if(ctx.rank() == 0){
            constexpr int source_rank = 1;
            MPI_Status status;
            
            
            advection_point_MPI(x, y, xadv, yadv, t, dt, T_period,0,nnx, backtrace_scheme);
            
            // receive the data from node-2
            MPI_Recv(&(A1[0][0]), ny*nnx, MPI_FLOAT, source_rank, 0, MPI_COMM_WORLD, &status);
            MPI_Recv(&(B1[0][0]), ny*nnx, MPI_FLOAT, source_rank, 0, MPI_COMM_WORLD, &status);
            
            
            // populate the full data to xadv and yadv
            for(unsigned int i = 0; i < ny; i++){  // loop for y - rows in the array
                for(unsigned int j = nnx; j < nx; j++){  // loop for x - columns in the array
                    xadv[i][j]=A1[j-nnx][i];  //the matrix received from node-2 is the transpose matrix
                    yadv[i][j]=B1[j-nnx][i];
                }
            }
            // done for MPI task
            
            
            intgridarray cellx, celly;
            cellx.resize(ny,intvectorarray(nx,0));
            celly.resize(ny,intvectorarray(nx,0));
            
            
            // do the task in node-1
            find_advection_point_location_MPI(x, y, xadv, yadv, cellx, celly, tracker, xlim1, xlim2, ylim1, ylim2,0,nnx);
            
            // receive the data from node-2
            MPI_Recv(&(A3[0][0]), ny*nnx, MPI_FLOAT, source_rank, 0, MPI_COMM_WORLD, &status);
            MPI_Recv(&(B3[0][0]), ny*nnx, MPI_FLOAT, source_rank, 0, MPI_COMM_WORLD, &status);
            MPI_Recv(&(C3[0][0]), ny*nnx, MPI_FLOAT, source_rank, 0, MPI_COMM_WORLD, &status);
            
            
            // populate the full data to xadv and yadv
            for(unsigned int i = 0; i < ny; i++){  // loop for y - rows in the array
                for(unsigned int j = nnx; j < nx; j++){  // loop for x - columns in the array
                    cellx[i][j]=A3[j-nnx][i];
                    celly[i][j]=B3[j-nnx][i];
                    tracker[i][j]=C3[j-nnx][i];
                }
            }
            // done for MPI task
            
            
            update_levelset_data(x, y, xadv, yadv, cellx, celly, tracker, t, dt,
                    tempphi, temppsix, temppsiy, temppsixy, mphi, mpsix, mpsiy, mpsixy,
                    psischeme,backtrace_scheme,T_period);
            
            //---------------------------------------------------------------------------------------------------------
            // Update the mixed derivatives now for the remaining grid points
            update_mixed_derivatives(temppsix, temppsiy, temppsixy, nx, ny, dx, dy);
            
            //---------------------------------------------------------------------------------------------------------
            // Feeding values back to the master matrix
            
            mphi = tempphi;
            mpsix = temppsix;
            mpsiy = temppsiy;
            mpsixy = temppsixy;
            
            //---------------------------------------------------------------------------------------------------------
            // Feeding phi, psix, psiy and psixy values in their respective files
            if((t+1) % printstep == 0)
                fileprint(mphi,mpsix,mpsiy,mpsixy,nx,ny,x,y,(t+1)*dt,T_period);
            
            cout<< t+1;
            cout<< " Time Step Completed" <<'\n';
            
            //---------------------------------------------------------------------------------------------------------
            xadv.clear();
            yadv.clear();
            tracker.clear();
            cellx.clear();
            celly.clear();
        }  // end of node-1
        
        
        // task for node-2
        else {
            
            constexpr int dest_rank = 0;  // We send a message to Task 0
            
            
            advection_point_MPI(x, y, xadv, yadv, t, dt, T_period, nnx,nx, backtrace_scheme);
            
            // populate data, here I transpose the matrix because for c++, matrix is read in row-wise order
            for(unsigned int i = 0; i < ny; i++){  // loop for y - rows in the array
                for(unsigned int j = nnx; j < nx; j++){  // loop for x - columns in the array
                    A1[j-nnx][i]=xadv[i][j];  // this is the transpose matrix
                    B1[j-nnx][i]=yadv[i][j];
                }
            }
            
            //send data to node-1
            MPI_Send(&(A1[0][0]), ny*nnx, MPI_FLOAT, dest_rank, 0, MPI_COMM_WORLD);
            MPI_Send(&(B1[0][0]), ny*nnx, MPI_FLOAT, dest_rank, 0, MPI_COMM_WORLD);
            
            
            
            intgridarray cellx, celly;
            cellx.resize(ny,intvectorarray(nx,0));
            celly.resize(ny,intvectorarray(nx,0));
            
            find_advection_point_location_MPI(x, y, xadv, yadv, cellx, celly, tracker, xlim1, xlim2, ylim1, ylim2,nnx,nx);
            
            // populate data, here I transpose the matrix because for c++, matrix is read in row-wise order
            for(unsigned int i = 0; i < ny; i++){  // loop for y - rows in the array
                for(unsigned int j = nnx; j < nx; j++){  // loop for x - columns in the array
                    A3[j-nnx][i]=cellx[i][j];
                    B3[j-nnx][i]=celly[i][j];
                    C3[j-nnx][i]=tracker[i][j];
                }
            }
            
            //send data to node-1
            MPI_Send(&(A3[0][0]), ny*nnx, MPI_FLOAT, dest_rank, 0, MPI_COMM_WORLD);
            MPI_Send(&(B3[0][0]), ny*nnx, MPI_FLOAT, dest_rank, 0, MPI_COMM_WORLD);
            MPI_Send(&(C3[0][0]), ny*nnx, MPI_FLOAT, dest_rank, 0, MPI_COMM_WORLD);
            
        }  	    // end of node-2
    } // end of time marching loop
    
    // end of timer, measure time
    gettimeofday(&end1,NULL);
    diff1 = 1000000 * (end1.tv_sec-start1.tv_sec)+ end1.tv_usec-start1.tv_usec;
    cout<<diff1*0.001<<endl;
    
    return 0;
}
