/*
 * This file contains all the functions that are used to initialize the level set, it's gradients and the velocity field on the 2D grid.
 * The last function is used to print all the details correctly in a file for all grid nodes. It generates 6 different txt files that are 
 * read usingthe Readfiles.m for visualization of the interface and levelset.
 */

#ifndef _Allocation_h
#define _Allocation_h

namespace galsfunctions
{
    void allocate_levelset_matrices(gridarray& mphi, gridarray& mpsix, gridarray& mpsiy, gridarray& mpsixy, vectorarray& x, vectorarray& y, unsigned int nx, unsigned int ny)
    {
        // master matrices of phi, psix, psiy, psixy which are updated at each time step & temporary matrices which is updated constantly throughout the cell iterations
        mphi.resize(ny,vectorarray(nx,0.0));
        mpsix.resize(ny,vectorarray(nx,0.0));
        mpsiy.resize(ny,vectorarray(nx,0.0));
        mpsixy.resize(ny,vectorarray(nx,0.0));
        
        for(unsigned int i = 0; i < ny; i++){
            for(unsigned int j = 0; j < nx; j++){
                mphi[i][j] = initialize(x[j], y[i]);
                mpsix[i][j] = derivxinit(x[j], y[i]);
                mpsiy[i][j] = derivyinit(x[j], y[i]);
                mpsixy[i][j] = derivxyinit(x[j], y[i]);
            }
        }
    }
    
    void allocate_levelset_velocity(gridarray& u, gridarray& v, unsigned int nx, unsigned int ny, vectorarray& x, vectorarray& y, double time, double T_period)
    {
        u.resize(ny,vectorarray(nx,0.0));
        v.resize(ny,vectorarray(nx,0.0));
        
        for(unsigned int i = 0; i < ny; i++){
            for(unsigned int j = 0; j < nx; j++){
                
                u[i][j] = Velx(x[j],y[i],time, T_period);
                v[i][j] = Vely(x[j],y[i],time, T_period);
            }
        }
    }
    
    void gridnodes(vectorarray& x, vectorarray& y, double xlim1, double ylim1, double dx, double dy, unsigned int nx, unsigned int ny)
    {
        // Defining node points
        x.resize(nx,0.0);
        y.resize(ny,0.0);
        for(unsigned int i = 0; i < nx; i++)
            x[i] = xlim1 + dx * i;
        for(unsigned int i = 0; i < ny; i++)
            y[i] = ylim1 + dy * i;
        // Node point definition ends
    }
    
    void fileprint(gridarray& mphi, gridarray& mpsix, gridarray& mpsiy, gridarray& mpsixy, unsigned int nx, unsigned int ny,
            vectorarray& x, vectorarray& y, double time, double T_period)
    {
        gridarray u, v;
        allocate_levelset_velocity(u,v,nx,ny,x,y,time,T_period);
        // Opening new files
        ofstream phifile, psixfile, psiyfile, psixyfile, ufile, vfile;
        phifile.open("phi.txt", ios::out | ios::app);
        psixfile.open("psix.txt", ios::out | ios::app);
        psiyfile.open("psiy.txt", ios::out | ios::app);
        psixyfile.open("psixy.txt", ios::out | ios::app);
        ufile.open("Velocity_x.txt", ios::out | ios::app);
        vfile.open("Velocity_y.txt",ios::out | ios::app);
        
        for(unsigned int i = 0; i < ny; i++){
            
            // Initializing the first element of each row separately to avoid trailing commas in the file
            phifile << std::fixed << std::setprecision(10) << mphi[i][0];
            psixfile << std::fixed << std::setprecision(10) << mpsix[i][0];
            psiyfile << std::fixed << std::setprecision(10) << mpsiy[i][0];
            psixyfile << std::fixed << std::setprecision(10) << mpsixy[i][0];
            ufile << std::fixed << std::setprecision(10) << u[i][0];
            vfile << std::fixed << std::setprecision(10) << v[i][0];
            
            for(unsigned int j = 1; j < nx; j++){
                phifile << "," << std::fixed << std::setprecision(10) << mphi[i][j];
                psixfile << "," << std::fixed << std::setprecision(10) << mpsix[i][j];
                psiyfile << "," << std::fixed << std::setprecision(10) << mpsiy[i][j];
                psixyfile << "," << std::fixed << std::setprecision(10) << mpsixy[i][j];
                ufile << "," << std::fixed << std::setprecision(10) << u[i][j];
                vfile << "," << std::fixed << std::setprecision(10) << v[i][j];
            }
            phifile << '\n';
            psixfile << '\n';
            psiyfile << '\n';
            psixyfile << '\n';
            ufile << '\n';
            vfile << '\n';
        }
        
        phifile << '\n';           psixfile << '\n';             psiyfile << '\n';            psixyfile << '\n';             ufile << '\n';         vfile << '\n';
        phifile.close();       psixfile.close();         psiyfile.close();          psixyfile.close();         ufile.close();        vfile.close();
    }
    
    
    // create contiguous-memory 2D matrices, A[i][j], used for MPI
    double **alloc_2d_int(unsigned int rows, unsigned int cols) {
        double *data = (double *)malloc(rows*cols*sizeof(double));
        double **array= (double **)malloc(rows*sizeof(double*));
        for (unsigned int i=0; i<rows; i++)
            array[i] = &(data[cols*i]);
        return array;
    }
}
#endif
