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

/* * * * * * * * * * * * * *  CUDA IMPLEMENTATION * * * * * * * * * * * * * * */
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
#include <cuda.h>
#include <chrono>
#include "Allocation.h"

//Including Kernel
#include "AdvectionPointCalcsCUDA.cu"
#include "VortexVelocityCUDA.cu"

// y direction is the first index of array, x direction is the second index of array

using namespace std;

int main(){
    
    /* UPDATE ALL THE FOLLOWING VALUES */
    double xlim1 = 0.0;                       //Lower limit on x-axis
    double xlim2 = 1.0;                      //Upper limit on x-axis
    unsigned int nx = 128;                         //Number of nodes in x-direction INCLUDING THE EXTREME VALUES
    
    double ylim1 = 0.0;                       //Lower limit on y-axis
    double ylim2 = 1.0;                     //Upper limit on y-axis
    unsigned int ny = 128;                        //Number of nodes INCLUDING THE EXTREME VALUES
    
    double dt = 0.5 * 1.0/128.0;                     //Length of time step
    double Tfinal = 1.0;                    //Total time period for the simulation
    unsigned int option = 1;                         //Option - if you need animation initialize at 1 else initialize at 2
    unsigned int printstep = 16;                      //How frequently do you want to store the images (every nth time step)
    char psischeme[] = "SuperConsistent";   //'SuperConsistent' or 'Heuns'
    char backtrace_scheme[] = "RK3" ;      //'Euler' or 'RK3'
    double T_period = 1.0;                  //Period of the velocity field
    unsigned int TileSize = 16;
    
    //---------------------------------------------------------------------------------------------------------
    //MAKE SURE THAT YOU HAVE ENOUGH MEMORY SPACE IF YOU ARE STORING A LOT OF TIME STEP VALUES BECAUSE IT STORES ACROSS GRID POINTS FOR EACH PRINTSTEP
    
    /* USER UPDATE OVER */
    unsigned long gridmemory = nx * ny * sizeof(double);
    unsigned long gridmemoryint = nx * ny * sizeof(unsigned int);
    unsigned int n = Tfinal/dt; //Number of time steps
    if(option != 1)
        printstep = n;
    
    dim3 dimBlock(TileSize, TileSize);
    dim3 dimGrid(nx/dimBlock.x, ny/dimBlock.y);
    
    // Node Locations
    double dx = (xlim2 - xlim1)/(nx - 1);
    double dy = (ylim2 - ylim1)/(ny - 1);
    double* x = (double*) malloc(nx * sizeof(double));
    double* y = (double*) malloc(ny * sizeof(double));
    gridnodes(x,y,xlim1,ylim1,dx,dy,nx,ny);
    double *devicex, *devicey;
    // allocate device memory for x and y
    cudaMalloc((void**)&devicex,nx * sizeof(double));	// Allocating GPU memory for the x-node values
    cudaMalloc((void**)&devicey,ny * sizeof(double));	// Allocating GPU memory for the y-node values
    // Copy data from host to GPU
    cudaMemcpy(devicex, x, nx * sizeof(double), cudaMemcpyHostToDevice);	// Writing to device memory
    cudaMemcpy(devicey, y, ny * sizeof(double), cudaMemcpyHostToDevice);	// Writing to device memory
    
    // level set matrices
    double* mphi = (double*) malloc(gridmemory);
    double* mpsix = (double*) malloc(gridmemory);
    double* mpsiy = (double*) malloc(gridmemory);
    double* mpsixy = (double*) malloc(gridmemory);
    
    double *masterdphi, *masterdpsix, *masterdpsiy, *masterdpsixy;
    
    // allocate device memory for integer grids
    cudaMalloc((void**)&masterdphi,gridmemory);	// Allocating GPU memory for the x-node values
    cudaMalloc((void**)&masterdpsix,gridmemory);	// Allocating GPU memory for the y-node values
    cudaMalloc((void**)&masterdpsiy,gridmemory);	// Allocating GPU memory for the y-node values
    cudaMalloc((void**)&masterdpsixy,gridmemory);	// Allocating GPU memory for the y-node values
    double *dphi, *dpsix, *dpsiy, *dpsixy;
    // allocate device memory for integer grids
    cudaMalloc((void**)&dphi,gridmemory);	// Allocating GPU memory for the x-node values
    cudaMalloc((void**)&dpsix,gridmemory);	// Allocating GPU memory for the y-node values
    cudaMalloc((void**)&dpsiy,gridmemory);	// Allocating GPU memory for the y-node values
    cudaMalloc((void**)&dpsixy,gridmemory);	// Allocating GPU memory for the y-node values
    
    // Initializing at t = 0
    allocate_levelset_matrices_CUDA<<<dimGrid, dimBlock>>>(masterdphi, masterdpsix, masterdpsiy, masterdpsixy, devicex, devicey, nx, ny); //Initializing level set matrices
    allocate_levelset_matrices_CUDA<<<dimGrid, dimBlock>>>(dphi, dpsix, dpsiy, dpsixy, devicex, devicey, nx, ny); //Initializing level set matrices
    // Initializing at t = 0
    allocate_levelset_matrices(mphi,mpsix,mpsiy,mpsixy,x,y,nx,ny); //Initializing level set matrices
    
    // Removing existing files with these names if any
    remove("phi.txt");
    remove("psix.txt");
    remove("psiy.txt");
    remove("psixy.txt");
    remove("details.txt");
    remove("Velocity_x.txt");
    remove("Velocity_y.txt");
    fileprint(mphi,mpsix,mpsiy,mpsixy,nx,ny,x,y,0.0,T_period);
    ofstream details;
    details.open("details.txt", ios::out | ios::app);
    details<< nx << "," << ny << "," << std::fixed << std::setprecision(10) << dx << "," << dy << "," << xlim1 << "," << xlim2 << "," << ylim1 << "," << ylim2 << "," << n << "," << dt << "," << printstep;
    details.close();
    
    ///*
    // TIME STEPPING LOOP
    // If only the initial and final profiles are needed
// This section will be deleted after a proper profiler is installed on my computer - Raunak Bardia
    cudaEvent_t startEvent, stopEvent;
    cudaEventCreate(&startEvent);
    cudaEventCreate(&stopEvent);
    float t_calcpts=0.0,t_findpts=0.0,t_update=0.0,t_mixed=0.0,t_copy=0.0,t_transfer=0.0,tempt=0.0;
    auto tbegin = chrono::high_resolution_clock::now();
    
//
    for(unsigned int t = 0; t < n; t++){
        
        double *dxadv, *dyadv;
        // allocate device memory for x and y
        cudaMalloc((void**)&dxadv,gridmemory);	// Allocating GPU memory for the x-node values
        cudaMalloc((void**)&dyadv,gridmemory);	// Allocating GPU memory for the y-node values
        
        unsigned int *dcellx, *dcelly,*dtracker;
        // allocate device memory for integer grids
        cudaMalloc((void**)&dcellx,gridmemoryint);	// Allocating GPU memory for the x-node values
        cudaMalloc((void**)&dcelly,gridmemoryint);	// Allocating GPU memory for the y-node values
        cudaMalloc((void**)&dtracker,gridmemoryint);	// Allocating GPU memory for the y-node values
        
        // Find the point from which advection occurs at this time step
        cudaEventRecord(startEvent,0);
        advection_point_cuda<<<dimGrid,dimBlock>>>(devicex,devicey,dxadv,dyadv,nx,t,dt,T_period,TileSize);
        cudaEventRecord(stopEvent,0);
        cudaEventSynchronize(stopEvent);
        cudaEventElapsedTime(&tempt, startEvent, stopEvent);
        t_calcpts += tempt;
        
        // Find the cell in which those advection points lie
        cudaEventRecord(startEvent,0);
        find_advection_point_location_cuda<<<dimGrid,dimBlock>>>(devicex,devicey,dxadv,dyadv,nx,ny,dcellx,dcelly,dtracker,xlim1,xlim2,ylim1,ylim2,TileSize);
        cudaEventRecord(stopEvent,0);
        cudaEventSynchronize(stopEvent);
        cudaEventElapsedTime(&tempt, startEvent, stopEvent);
        t_findpts += tempt;
        
        // Update the level set values
        cudaEventRecord(startEvent,0);
        update_levelset_data_cuda<<<dimGrid,dimBlock>>>(devicex, devicey, dxadv, dyadv, dcellx, dcelly, dtracker, t, dt, dphi, dpsix, dpsiy, dpsixy, masterdphi, masterdpsix, masterdpsiy,psischeme,backtrace_scheme,T_period,nx,ny,TileSize);
        cudaEventRecord(stopEvent,0);
        cudaEventSynchronize(stopEvent);
        cudaEventElapsedTime(&tempt, startEvent, stopEvent);
        t_update += tempt;
        
        // Create another copy to preserve data which gets modified on the fly in next loop
        cudaEventRecord(startEvent,0);
        devicetodevicecopy<<<dimGrid,dimBlock>>>(dphi,dpsix,dpsiy,masterdphi,masterdpsix,masterdpsiy,nx,TileSize);
        cudaEventRecord(stopEvent,0);
        cudaEventSynchronize(stopEvent);
        cudaEventElapsedTime(&tempt, startEvent, stopEvent);
        t_copy += tempt;
        
        // Update the mixed derivatives now for the remaining grid points
        cudaEventRecord(startEvent,0);
        update_mixed_derivatives<<<dimGrid,dimBlock>>>(dpsix, dpsiy, dpsixy, nx, ny, dx, dy,TileSize);
        cudaEventRecord(stopEvent,0);
        cudaEventSynchronize(stopEvent);
        cudaEventElapsedTime(&tempt, startEvent, stopEvent);
        t_mixed += tempt;
        
        cudaDeviceSynchronize();
        
        //---------------------------------------------------------------------------------------------------------
        // Feeding phi, psix, psiy and psixy values in their respective files
        if((t+1) % printstep == 0)
        {
            cudaEventRecord(startEvent,0);
            cudaMemcpy(mphi, masterdphi, gridmemory, cudaMemcpyDeviceToHost);       // Writing back to host memory
            cudaMemcpy(mpsix, masterdpsix, gridmemory, cudaMemcpyDeviceToHost);       // Writing back to host memory
            cudaMemcpy(mpsiy, masterdpsiy, gridmemoryint, cudaMemcpyDeviceToHost);      // Writing back to host memory
            cudaEventRecord(stopEvent,0);
            cudaEventSynchronize(stopEvent);
            cudaEventElapsedTime(&tempt, startEvent, stopEvent);
            t_transfer += tempt;
            //cudaMemcpy(mpsixy, masterdpsixy, gridmemoryint, cudaMemcpyDeviceToHost);  // Writing back to host memory
            fileprint(mphi,mpsix,mpsiy,mpsixy,nx,ny,x,y,(t+1)*dt,T_period);
        }
        cout<< t+1;
        cout<< " Time Step Completed" <<'\n';
        
        //---------------------------------------------------------------------------------------------------------
        cudaFree(dxadv);
        cudaFree(dyadv);
        cudaFree(dcellx);
        cudaFree(dcelly);
        cudaFree(dtracker);
    }  // end of time marching loop
    //*/
    auto tend = chrono::high_resolution_clock::now();
    float duration = chrono::duration_cast<chrono::nanoseconds>(tend-tbegin).count();
    duration = duration * pow(10.0,-6);
    cout << "Time taken for calculation of advection points = " << t_calcpts << '\n';
    cout << "Time taken for finding locatio of advection points in the grid = " << t_findpts << '\n';
    cout << "Time taken for hermite Update = " << t_update << '\n';
    cout << "Time taken for copying the matrix = " << t_copy << '\n';
    cout << "Time taken for calculation of mixed derivatives = " << t_mixed << '\n';
    cout << "Time taken for transfer of level set data = " << t_transfer << '\n';
    cout << "Total Duration of the Time Loop = " << duration << endl;
    
    return 0;
}
