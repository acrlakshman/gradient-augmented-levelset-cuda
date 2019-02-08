/*
This file is an implementation of all the functions that work on a single node
in the 2D grid and update the level set.
*/

#include "HermiteCUDA.cu"
#include "TimeSteppingMethodsCUDA.cu"

#ifndef _AdvectionPointCalcsCUDA_h
#define _AdvectionPointCalcsCIDA_h

__global__ void advection_point_cuda(double *x, double *y, double *xadv, double *yadv, unsigned int nx,
        unsigned int t, double dt, double T_period, unsigned int TileSize)
{
    double c1 = (1/6.0);
    double c2 = (1/6.0);
    double c3 = (2/3.0);    //RK-3 constants
    
    unsigned int bx = blockIdx.x;
    unsigned int by = blockIdx.y;
    
    unsigned int tx = threadIdx.x;
    unsigned int ty = threadIdx.y;
    
    unsigned int index_x = bx * TileSize + tx;
    unsigned int index_y = by * TileSize + ty;
    
    double ux = Velx(x[index_x], y[index_y], (t + 1) * dt, T_period);
    double vy = Vely(x[index_x], y[index_y], (t + 1) * dt, T_period);
    
    unsigned int indexToWrite = index_y * nx + index_x;
    
    // Advected Points - 1 Step
    xadv[indexToWrite] = x[index_x] - ux * dt;
    yadv[indexToWrite] = y[index_y] - vy * dt;
    
    double ux1 = Velx(xadv[indexToWrite], yadv[indexToWrite], t * dt, T_period);
    double vy1 = Vely(xadv[indexToWrite], yadv[indexToWrite], t * dt, T_period);
    
    // 2 Step
    xadv[indexToWrite] = x[index_x] - ((dt)/2.0 * (0.5 * ux + 0.5 * ux1));
    yadv[indexToWrite] = y[index_y] - ((dt)/2.0 * (0.5 * vy + 0.5 * vy1));
    
    double ux2 = Velx(xadv[indexToWrite], yadv[indexToWrite],(t * dt + dt/2.0), T_period);
    double vy2 = Vely(xadv[indexToWrite], yadv[indexToWrite],(t * dt + dt/2.0), T_period);
    
    // 3 Step
    xadv[indexToWrite] = x[index_x] - (dt * (c1 * ux + c2 * ux1 + c3 * ux2));
    yadv[indexToWrite] = y[index_y] - (dt * (c1 * vy + c2 * vy1 + c3 * vy2));
}

__device__ unsigned int locationAlgo(double *x, double xadv, unsigned int nx)
{
    unsigned int location = 0;
    while (x[location] < xadv && location < nx)
        location++;
    if(location == 0)
        return location;
    else
        return location-1;
}

__global__ void find_advection_point_location_cuda(double *x, double *y, double *xadv, double *yadv, unsigned int nx, unsigned int ny,
        unsigned int *cellx, unsigned int *celly, unsigned int *tracker, double xlim1, double xlim2, double ylim1,
        double ylim2, unsigned int TileSize)
{
    unsigned int bx = blockIdx.x;
    unsigned int by = blockIdx.y;
    
    unsigned int tx = threadIdx.x;
    unsigned int ty = threadIdx.y;
    
    unsigned int index_x = bx * TileSize + tx;
    unsigned int index_y = by * TileSize + ty;
    
    unsigned int indexToWrite = index_y * nx + index_x;
    
    bool xoutofbounds = false;
    bool youtofbounds = false;
    
    if(!((xadv[indexToWrite] > xlim1) && (xadv[indexToWrite] < xlim2)))
        xoutofbounds = true;
    if(!((yadv[indexToWrite] > ylim1) && (yadv[indexToWrite] < ylim2)))
        youtofbounds = true;
    
    if(!xoutofbounds && !youtofbounds)
    {
        tracker[indexToWrite] = 1;
        cellx[indexToWrite] = locationAlgo(x,xadv[indexToWrite],nx);
        celly[indexToWrite] = locationAlgo(y,yadv[indexToWrite],ny);
    }
    else
        if(!xoutofbounds && youtofbounds)
        {
            tracker[indexToWrite] = 2;
            cellx[indexToWrite] = locationAlgo(x,xadv[indexToWrite],nx);
            if(yadv[indexToWrite] <= ylim1)
                celly[indexToWrite] = 0;
            else
                if(yadv[indexToWrite] >= ylim2)
                    celly[indexToWrite] = ny-2;
        }
        else
            if(xoutofbounds && !youtofbounds)
            {
                tracker[indexToWrite] = 3;
                celly[indexToWrite] = locationAlgo(y,yadv[indexToWrite],ny);
                if(xadv[indexToWrite] <= xlim1)
                    cellx[indexToWrite] = 0;
                else
                    if(xadv[indexToWrite] >= xlim2)
                        cellx[indexToWrite] = nx-2;
            }
            else
                if(xoutofbounds && youtofbounds)
                    tracker[indexToWrite] = 4;
}

__global__ void update_levelset_data_cuda(double *x, double *y, double *xadv, double *yadv, unsigned int *cellx,
        unsigned int *celly, unsigned int *tracker, unsigned int t, double dt, double *tempphi, double *temppsix,
        double *temppsiy, double *temppsixy, double *mphi, double *mpsix, double *mpsiy, char psischeme[],
        char backtrace_scheme[], double T_period, unsigned int nx, unsigned int ny, unsigned int TileSize)
{
    unsigned int bx = blockIdx.x;
    unsigned int by = blockIdx.y;
    
    unsigned int tx = threadIdx.x;
    unsigned int ty = threadIdx.y;
    
    unsigned int index_x = bx * TileSize + tx;
    unsigned int index_y = by * TileSize + ty;
    
    unsigned int indexToWrite = index_y * nx + index_x;
    
    double dx = x[2] - x[1];
    double dy = y[2] - y[1];
    
    double phi[4], psix[4], psiy[4], psixy[4];
    
    unsigned int cellindex_x = cellx[indexToWrite];
    unsigned int cellindex_y = celly[indexToWrite];
    
    unsigned int cellindex = cellindex_x + cellindex_y * nx;
    
    // Storing the four values for four nodes of each cell
    phi[0] = mphi[cellindex];            psix[0] = mpsix[cellindex];          psiy[0] = mpsiy[cellindex];          psixy[0] = temppsixy[cellindex];
    phi[1] = mphi[cellindex + 1];          psix[1] = mpsix[cellindex + 1];        psiy[1] = mpsiy[cellindex + 1];        psixy[1] = temppsixy[cellindex + 1];
    phi[2] = mphi[cellindex + nx];          psix[2] = mpsix[cellindex + nx];        psiy[2] = mpsiy[cellindex + nx];        psixy[2] = temppsixy[cellindex + nx];
    phi[3] = mphi[cellindex + nx + 1];        psix[3] = mpsix[cellindex + nx + 1];      psiy[3] = mpsiy[cellindex + nx + 1];      psixy[3] = temppsixy[cellindex + nx + 1];
    // Node value assignment ends
    
    // Storing the coordinates of first node of the working cell
    double xo = x[cellindex_x], yo = y[cellindex_y];
    
    switch(tracker[indexToWrite]){
        case 1:{
            tempphi[indexToWrite] = hp(phi, psix, psiy, psixy, xadv[indexToWrite], yadv[indexToWrite], xo, yo, dx, dy);
            double rootpsix = hermx(phi,psix,psiy,psixy,xadv[indexToWrite],yadv[indexToWrite],xo,yo,dx, dy);
            double rootpsiy = hermy(phi,psix,psiy,psixy,xadv[indexToWrite],yadv[indexToWrite],xo,yo,dx, dy);
            
            //Commenting out the options of psi scheme
            //SuperConsistent is being used now irrespective of specification given in
            //GALS_Advection.cu
            //The if/else functionality may be added later
            //if(strcmp("Heuns",psischeme) == 0)
            //Heuns_internal(x[index_x],y[index_y],xadv[indexToWrite],yadv[indexToWrite],rootpsix,rootpsiy,t,dt,T_period,temppsix,temppsiy,indexToWrite);
            //else if(strcmp("SuperConsistent",psischeme) == 0)
            SuperConsistentScheme(x[index_x],y[index_y],rootpsix,rootpsiy,t,dt,T_period,backtrace_scheme,temppsix,temppsiy,indexToWrite);
            
            break;
        } // end of case 1
        
        case 2:{
            double rootpsix = hermx(phi,psix,psiy,psixy,xadv[indexToWrite], y[index_y],xo,yo,dx, dy);
            double rootpsiy = temppsiy[indexToWrite];
            tempphi[indexToWrite] = hp(phi, psix, psiy, psixy, xadv[indexToWrite], y[index_y], xo, yo, dx, dy) - dt * Vely(x[index_x], y[index_y], t * dt, T_period) * rootpsiy;
            temppsix[indexToWrite] = Heuns_X(x[index_x],y[index_y],rootpsix,rootpsiy,temppsixy[indexToWrite],temppsix[indexToWrite],t,dt,T_period);
            break;
        }   // end of case 2
        
        case 3:{
            double rootpsix = temppsix[indexToWrite];
            double rootpsiy = hermy(phi,psix,psiy,psixy,x[index_x],yadv[indexToWrite],xo,yo,dx, dy);
            tempphi[indexToWrite] = hp(phi, psix, psiy, psixy, x[index_x],yadv[indexToWrite], xo, yo, dx, dy) - dt * Velx(x[index_x], y[index_y], t * dt, T_period) * rootpsix;
            temppsiy[indexToWrite] = Heuns_Y(x[index_x],y[index_y],rootpsix,rootpsiy,temppsixy[indexToWrite],temppsiy[indexToWrite],t,dt,T_period);
            break;
        }   // end of case 3
        
        case 4:{
            tempphi[indexToWrite] = tempphi[indexToWrite] - dt * (Velx(x[index_x], y[index_y], t * dt, T_period) * temppsix[indexToWrite] + Vely(x[index_x], y[index_y], t * dt, T_period) * temppsiy[indexToWrite]);
            break;
        } //end of case4
        
        default:{break;}
    }   // end of switch loop
}

__global__ void update_mixed_derivatives(double *temppsix, double *temppsiy, double *temppsixy,
        unsigned int nx, unsigned int ny, double dx, double dy, unsigned int TileSize)
{
    unsigned int bx = blockIdx.x;
    unsigned int by = blockIdx.y;
    
    unsigned int tx = threadIdx.x;
    unsigned int ty = threadIdx.y;
    
    unsigned int index_x = bx * TileSize + tx;
    unsigned int index_y = by * TileSize + ty;
    
    unsigned int indexToWrite = index_y * nx + index_x;
    
    if ((index_y == 0 || index_y == ny - 1) && (index_x != 0 && index_x != nx - 1))
        temppsixy[indexToWrite] = (temppsiy[indexToWrite+1] - temppsiy[indexToWrite-1])/(2 * dx);
    else
        if ((index_y != 0 && index_y != ny - 1) && (index_x == 0 || index_x == nx - 1))
            temppsixy[indexToWrite] = (temppsix[indexToWrite + nx] - temppsix[indexToWrite - nx])/(2 * dy);
        else
            if((index_y == 0 || index_y == ny - 1) && (index_x == 0 || index_x == nx - 1)){
                if(index_y == 0 && index_x == 0){
                    double d1 = (temppsiy[1] - temppsiy[0])/dx;
                    double d2 = (temppsix[nx] - temppsix[0])/dy;
                    double d3 = (temppsix[nx+1] - temppsix[1])/dy;
                    double d4 = (temppsiy[nx+1] - temppsiy[nx])/dx;
                    temppsixy[indexToWrite] = 0.75 * (d1 + d2) - 0.25 * (d3 + d4);
                }
                else if(index_y == 0 && index_x == nx-1){
                    double d1 = (temppsiy[nx-1] - temppsiy[nx-2])/dx;
                    double d2 = (temppsix[nx+nx-2] - temppsix[nx-2])/dy;
                    double d3 = (temppsix[nx+nx-1] - temppsix[nx-1])/dy;
                    double d4 = (temppsiy[nx+nx-1] - temppsiy[nx+nx-2])/dx;
                    temppsixy[indexToWrite] = 0.75 * (d1 + d3) - 0.25 * (d2 + d4);
                    
                }
                else if(index_y == ny-1 && index_x == 0){
                    double d1 = (temppsiy[nx *(ny-2) + 1] - temppsiy[nx *(ny-2)])/dx;
                    double d2 = (temppsix[nx *(ny-1)] - temppsix[nx *(ny-2)])/dy;
                    double d3 = (temppsix[nx *(ny-1)] - temppsix[nx *(ny-2) + 1])/dy;
                    double d4 = (temppsiy[nx *(ny-1) + 1] - temppsiy[nx *(ny-1)])/dx;
                    temppsixy[indexToWrite] = 0.75 * (d2 + d4) - 0.25 * (d3 + d1);
                    
                }
                else if(index_y == ny-1 && index_x == nx-1){
                    double d1 = (temppsiy[nx *(ny-2) + nx - 1] - temppsiy[nx *(ny-2) + nx - 2])/dx;
                    double d2 = (temppsix[nx *(ny-1) + nx - 2] - temppsix[nx *(ny-2) + nx - 2])/dy;
                    double d3 = (temppsix[nx *(ny-1) + nx - 1] - temppsix[nx *(ny-2) + nx - 1])/dy;
                    double d4 = (temppsiy[nx *(ny-1) + nx - 1] - temppsiy[nx *(ny-1) + nx - 2])/dx;
                    temppsixy[indexToWrite] = 0.75 * (d3 + d4) - 0.25 * (d1 + d2);
                }
            }
            else{
                double dxy1 = (temppsiy[indexToWrite+1] - temppsiy[indexToWrite-1])/(2 * dx);
                double dxy2 = (temppsix[indexToWrite + nx] - temppsix[indexToWrite - nx])/(2 * dy);
                temppsixy[indexToWrite] = (dxy1 + dxy2)/2.0;
            }
    
}
#endif
