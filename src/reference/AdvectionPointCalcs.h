/*
This file is an implementation of all the functions that work on a single node
in the 2D grid and update the level set.
*/

#ifndef _AdvectionPointCalcs_h
#define _AdvectionPointCalcs_h

namespace galsfunctions
{

   // used for MPI
   void advection_point_MPI(vectorarray& x, vectorarray& y, gridarray& xadv, gridarray &yadv, unsigned int t,
                         double dt, double T_period, unsigned int ny_Node_min, unsigned int ny_Node_max, char backtrace_scheme[])
    {
        double c1 = (1/6.0);
        double c2 = (1/6.0);
        double c3 = (2/3.0);    //RK-3 constants

        unsigned int ny = y.size();
	# pragma omp parallel
	{
        # pragma omp for schedule(static) collapse(2)  
        for(unsigned int i = 0; i < ny; i++){  // loop for y - rows in the array
            for(unsigned int j = ny_Node_min; j < ny_Node_max; j++){  // loop for x - columns in the array

                double ux = Velx(x[j], y[i], (t + 1) * dt, T_period);
                double vy = Vely(x[j], y[i], (t + 1) * dt, T_period);

                // Advected Points - 1 Step
                   xadv[i][j] = x[j] - ux * dt;
                   yadv[i][j] = y[i] - vy * dt;
                
                  if(strcmp("RK3",backtrace_scheme) == 0)
                   {
                    double ux1 = Velx(xadv[i][j], yadv[i][j], t * dt, T_period);
                    double vy1 = Vely(xadv[i][j], yadv[i][j], t * dt, T_period);

                   // 2 Step
                   xadv[i][j] = x[j] - ((dt)/2.0 * (0.5 * ux + 0.5 * ux1));
                   yadv[i][j] = y[i] - ((dt)/2.0 * (0.5 * vy + 0.5 * vy1));
                
                   double ux2 = Velx(xadv[i][j], yadv[i][j],(t * dt + dt/2.0), T_period);
                   double vy2 = Vely(xadv[i][j], yadv[i][j],(t * dt + dt/2.0), T_period);

 		   // 3 Step
                   xadv[i][j] = x[j] - (dt * (c1 * ux + c2 * ux1 + c3 * ux2));
                   yadv[i][j] = y[i] - (dt * (c1 * vy + c2 * vy1 + c3 * vy2));          
		   }
	     }
	 }
     }
	}

    void advection_point(vectorarray& x, vectorarray& y, gridarray& xadv, gridarray &yadv, unsigned int t,
                         double dt, double T_period, char backtrace_scheme[])
    {
        double c1 = (1/6.0);
        double c2 = (1/6.0);
        double c3 = (2/3.0);    //RK-3 constants

        unsigned int nx = x.size();
        unsigned int ny = y.size();

        for(unsigned int i = 0; i < ny - 1; i++){  // loop for y - rows in the array
            for(unsigned int j = 0; j < nx - 1; j++){  // loop for x - columns in the array

                double ux = Velx(x[j], y[i], (t + 1) * dt, T_period);
                double vy = Vely(x[j], y[i], (t + 1) * dt, T_period);

                // Advected Points - 1 Step
                xadv[i][j] = x[j] - ux * dt;
                yadv[i][j] = y[i] - vy * dt;

                if(strcmp("RK3",backtrace_scheme) == 0)
                {
                    double ux1 = Velx(xadv[i][j], yadv[i][j], t * dt, T_period);
                    double vy1 = Vely(xadv[i][j], yadv[i][j], t * dt, T_period);

                    // 2 Step
                    xadv[i][j] = x[j] - ((dt)/2.0 * (0.5 * ux + 0.5 * ux1));
                    yadv[i][j] = y[i] - ((dt)/2.0 * (0.5 * vy + 0.5 * vy1));

                    double ux2 = Velx(xadv[i][j], yadv[i][j],(t * dt + dt/2.0), T_period);
                    double vy2 = Vely(xadv[i][j], yadv[i][j],(t * dt + dt/2.0), T_period);

                    // 3 Step
                    xadv[i][j] = x[j] - (dt * (c1 * ux + c2 * ux1 + c3 * ux2));
                    yadv[i][j] = y[i] - (dt * (c1 * vy + c2 * vy1 + c3 * vy2));
                }
            }
        }
    }

    void find_advection_point_location (vectorarray &x, vectorarray &y, gridarray &xadv, gridarray &yadv,
            intgridarray &cellx, intgridarray &celly, intgridarray &tracker,
	    double xlim1, double xlim2, double ylim1, double ylim2)
    {
        unsigned int nx = x.size();
        unsigned int ny = y.size();

        for(unsigned int i = 0; i < ny; i++){
            for(unsigned int j = 0; j < nx; j++){
                std::vector<double>::iterator locate;
                bool xoutofbounds = false;
                bool youtofbounds = false;

                if(!((xadv[i][j] > xlim1) && (xadv[i][j] < xlim2)))
                    xoutofbounds = true;
                if(!((yadv[i][j] > ylim1) && (yadv[i][j] < ylim2)))
                    youtofbounds = true;

                if(!xoutofbounds && !youtofbounds)
                {
                    tracker[i][j] = 1;
                    locate = std::lower_bound(x.begin(), x.end(), xadv[i][j]);
                    cellx[i][j] = locate - x.begin() - 1;
                    locate = std::lower_bound(y.begin(), y.end(), yadv[i][j]);
                    celly[i][j] = locate - y.begin() - 1;
                }
                else
                    if(!xoutofbounds && youtofbounds)
                    {
                        tracker[i][j] = 2;
                        if(yadv[i][j] <= ylim1)
                        {
                            locate = std::lower_bound(x.begin(), x.end(), xadv[i][j]);
                            cellx[i][j] = locate - x.begin() - 1;
                            celly[i][j] = 0;
                        }
                        else
                            if(yadv[i][j] >= ylim2)
                            {
                                locate = std::lower_bound(x.begin(), x.end(), xadv[i][j]);
                                cellx[i][j] = locate - x.begin() - 1;
                                celly[i][j] = ny-2;
                            }
                    }
                    else
                        if(xoutofbounds && !youtofbounds)
                        {
                            tracker[i][j] = 3;
                            if(xadv[i][j] <= xlim1)
                            {
                                locate = std::lower_bound(y.begin(), y.end(), yadv[i][j]);
                                celly[i][j] = locate - y.begin() - 1;
                                cellx[i][j] = 0;
                            }
                            else
                                if(xadv[i][j] >= xlim2)
                                {
                                    locate = std::lower_bound(y.begin(), y.end(), yadv[i][j]);
                                    celly[i][j] = locate - y.begin() - 1;
                                    cellx[i][j] = nx-2;
                                }
                        }
                        else
                            if(xoutofbounds && youtofbounds)
                                tracker[i][j] = 4;
            }
        }
    }

    void find_advection_point_location_MPI (vectorarray &x, vectorarray &y, gridarray &xadv, gridarray &yadv,
            intgridarray &cellx, intgridarray &celly, intgridarray &tracker,
	    double xlim1, double xlim2, double ylim1, double ylim2,unsigned int ny_Node_min, unsigned int ny_Node_max)
    {
        unsigned int nx = x.size();
        unsigned int ny = y.size();
	# pragma omp parallel
        {
        # pragma omp for schedule(static) collapse(2)   
        for(unsigned int i = 0; i < ny; i++){
            for(unsigned int j = ny_Node_min; j < ny_Node_max; j++){
                std::vector<double>::iterator locate;
                bool xoutofbounds = false;
                bool youtofbounds = false;

                if(!((xadv[i][j] > xlim1) && (xadv[i][j] < xlim2)))
                    xoutofbounds = true;
                if(!((yadv[i][j] > ylim1) && (yadv[i][j] < ylim2)))
                    youtofbounds = true;

                if(!xoutofbounds && !youtofbounds)
                {
                    tracker[i][j] = 1;
                    locate = std::lower_bound(x.begin(), x.end(), xadv[i][j]);
                    cellx[i][j] = locate - x.begin() - 1;
                    locate = std::lower_bound(y.begin(), y.end(), yadv[i][j]);
                    celly[i][j] = locate - y.begin() - 1;
                }
                else
                    if(!xoutofbounds && youtofbounds)
                    {
                        tracker[i][j] = 2;
                        if(yadv[i][j] <= ylim1)
                        {
                            locate = std::lower_bound(x.begin(), x.end(), xadv[i][j]);
                            cellx[i][j] = locate - x.begin() - 1;
                            celly[i][j] = 0;
                        }
                        else
                            if(yadv[i][j] >= ylim2)
                            {
                                locate = std::lower_bound(x.begin(), x.end(), xadv[i][j]);
                                cellx[i][j] = locate - x.begin() - 1;
                                celly[i][j] = ny-2;
                            }
                    }
                    else
                        if(xoutofbounds && !youtofbounds)
                        {
                            tracker[i][j] = 3;
                            if(xadv[i][j] <= xlim1)
                            {
                                locate = std::lower_bound(y.begin(), y.end(), yadv[i][j]);
                                celly[i][j] = locate - y.begin() - 1;
                                cellx[i][j] = 0;
                            }
                            else
                                if(xadv[i][j] >= xlim2)
                                {
                                    locate = std::lower_bound(y.begin(), y.end(), yadv[i][j]);
                                    celly[i][j] = locate - y.begin() - 1;
                                    cellx[i][j] = nx-2;
                                }
                        }
                        else
                            if(xoutofbounds && youtofbounds)
                                tracker[i][j] = 4;
            }
        }
    }
	}

    void update_levelset_data_mpi(vectorarray &x, vectorarray &y, gridarray &xadv, gridarray &yadv,
            intgridarray &cellx, intgridarray &celly, intgridarray &tracker, unsigned int t, double dt,
            gridarray &tempphi, gridarray &temppsix, gridarray &temppsiy, gridarray &temppsixy,
            gridarray &mphi, gridarray &mpsix, gridarray &mpsiy, gridarray &mpsixy, char psischeme[],
            char backtrace_scheme[], unsigned int ny_Node_min, unsigned int ny_Node_max, double T_period){

        double dx = x[2] - x[1];
        unsigned int ny = y.size();
        double dy = y[2] - y[1];

        for(unsigned int i = 0; i < ny; i++){  // loop for y - rows in the array
            for(unsigned int j = ny_Node_min; j < ny_Node_max; j++){  // loop for x - columns in the array

                double phi[4], psix[4], psiy[4], psixy[4];

                unsigned int cellindex_x = cellx[i][j], cellindex_y = celly[i][j];
                // Storing the four values for four nodes of each cell
                phi[0] = mphi[cellindex_y][cellindex_x];            psix[0] = mpsix[cellindex_y][cellindex_x];          psiy[0] = mpsiy[cellindex_y][cellindex_x];          psixy[0] = mpsixy[cellindex_y][cellindex_x];
                phi[1] = mphi[cellindex_y][cellindex_x+1];          psix[1] = mpsix[cellindex_y][cellindex_x+1];        psiy[1] = mpsiy[cellindex_y][cellindex_x+1];        psixy[1] = mpsixy[cellindex_y][cellindex_x+1];
                phi[2] = mphi[cellindex_y+1][cellindex_x];          psix[2] = mpsix[cellindex_y+1][cellindex_x];        psiy[2] = mpsiy[cellindex_y+1][cellindex_x];        psixy[2] = mpsixy[cellindex_y+1][cellindex_x];
                phi[3] = mphi[cellindex_y+1][cellindex_x+1];        psix[3] = mpsix[cellindex_y+1][cellindex_x+1];      psiy[3] = mpsiy[cellindex_y+1][cellindex_x+1];      psixy[3] = mpsixy[cellindex_y+1][cellindex_x+1];
                // Node value assignment ends

                // Storing the coordinates of first node of the working cell
                double xo = x[cellindex_x], yo = y[cellindex_y];

                switch(tracker[i][j]){
                case 1:{
                    tempphi[i][j] = hp(phi, psix, psiy, psixy, xadv[i][j], yadv[i][j], xo, yo, dx, dy);
                    double rootpsix = hermx(phi,psix,psiy,psixy,xadv[i][j],yadv[i][j],xo,yo,dx, dy);
                    double rootpsiy = hermy(phi,psix,psiy,psixy,xadv[i][j],yadv[i][j],xo,yo,dx, dy);

                    if(strcmp("Heuns",psischeme) == 0)
                        std::tie(temppsix[i][j],temppsiy[i][j]) = Heuns_internal(x[j],y[i],xadv[i][j],yadv[i][j],rootpsix,rootpsiy,t,dt,T_period);
                    else if(strcmp("SuperConsistent",psischeme) == 0)
                        std::tie(temppsix[i][j],temppsiy[i][j]) = SuperConsistentScheme(x[j],y[i],rootpsix,rootpsiy,t,dt,T_period,backtrace_scheme);

                    break;
                } // end of case 1

                case 2:{
                    double rootpsix = hermx(phi,psix,psiy,psixy,xadv[i][j], y[i],xo,yo,dx, dy);
                    double rootpsiy = temppsiy[i][j];
                    tempphi[i][j] = hp(phi, psix, psiy, psixy, xadv[i][j], y[i], xo, yo, dx, dy) - dt * Vely(x[j], y[i], t * dt, T_period) * rootpsiy;
                    temppsix[i][j] = Heuns_X(x[j],y[i],rootpsix,rootpsiy,temppsixy[i][j],temppsix[i][j],t,dt,T_period);
                    break;
                }   // end of case 2

                case 3:{
                    double rootpsix = temppsix[i][j];
                    double rootpsiy = hermy(phi,psix,psiy,psixy,x[j],yadv[i][j],xo,yo,dx, dy);
                    tempphi[i][j] = hp(phi, psix, psiy, psixy, x[j],yadv[i][j], xo, yo, dx, dy) - dt * Velx(x[j], y[i], t * dt, T_period) * rootpsix;
                    temppsiy[i][j] = Heuns_Y(x[j],y[i],rootpsix,rootpsiy,temppsixy[i][j],temppsiy[i][j],t,dt,T_period);
                    break;
                }   // end of case 3

                case 4:{
                    tempphi[i][j] = tempphi[i][j] - dt * (Velx(x[j], y[i], t * dt, T_period) * temppsix[i][j] + Vely(x[j], y[i], t * dt, T_period) * temppsiy[i][j]);
                    break;
                } //end of case4

                default:{break;}
                }   // end of switch loop
            }
        }
    }


    void update_levelset_data(vectorarray &x, vectorarray &y, gridarray &xadv, gridarray &yadv,
            intgridarray &cellx, intgridarray &celly, intgridarray &tracker, unsigned int t, double dt,
            gridarray &tempphi, gridarray &temppsix, gridarray &temppsiy, gridarray &temppsixy,
            gridarray &mphi, gridarray &mpsix, gridarray &mpsiy, gridarray &mpsixy, char psischeme[],
            char backtrace_scheme[], double T_period){

        unsigned int nx = x.size();
        double dx = x[2] - x[1];
        unsigned int ny = y.size();
        double dy = y[2] - y[1];
        # pragma omp parallel
        {
        # pragma omp for schedule(static) collapse(2)    

        for(unsigned int i = 0; i < ny; i++){  // loop for y - rows in the array
            for(unsigned int j = 0; j < nx; j++){  // loop for x - columns in the array

                double phi[4], psix[4], psiy[4], psixy[4];

                unsigned int cellindex_x = cellx[i][j], cellindex_y = celly[i][j];
                // Storing the four values for four nodes of each cell
                phi[0] = mphi[cellindex_y][cellindex_x];            psix[0] = mpsix[cellindex_y][cellindex_x];          psiy[0] = mpsiy[cellindex_y][cellindex_x];          psixy[0] = mpsixy[cellindex_y][cellindex_x];
                phi[1] = mphi[cellindex_y][cellindex_x+1];          psix[1] = mpsix[cellindex_y][cellindex_x+1];        psiy[1] = mpsiy[cellindex_y][cellindex_x+1];        psixy[1] = mpsixy[cellindex_y][cellindex_x+1];
                phi[2] = mphi[cellindex_y+1][cellindex_x];          psix[2] = mpsix[cellindex_y+1][cellindex_x];        psiy[2] = mpsiy[cellindex_y+1][cellindex_x];        psixy[2] = mpsixy[cellindex_y+1][cellindex_x];
                phi[3] = mphi[cellindex_y+1][cellindex_x+1];        psix[3] = mpsix[cellindex_y+1][cellindex_x+1];      psiy[3] = mpsiy[cellindex_y+1][cellindex_x+1];      psixy[3] = mpsixy[cellindex_y+1][cellindex_x+1];
                // Node value assignment ends

                // Storing the coordinates of first node of the working cell
                double xo = x[cellindex_x], yo = y[cellindex_y];

                switch(tracker[i][j]){
                case 1:{
                    tempphi[i][j] = hp(phi, psix, psiy, psixy, xadv[i][j], yadv[i][j], xo, yo, dx, dy);
                    double rootpsix = hermx(phi,psix,psiy,psixy,xadv[i][j],yadv[i][j],xo,yo,dx, dy);
                    double rootpsiy = hermy(phi,psix,psiy,psixy,xadv[i][j],yadv[i][j],xo,yo,dx, dy);

                    if(strcmp("Heuns",psischeme) == 0)
                        std::tie(temppsix[i][j],temppsiy[i][j]) = Heuns_internal(x[j],y[i],xadv[i][j],yadv[i][j],rootpsix,rootpsiy,t,dt,T_period);
                    else if(strcmp("SuperConsistent",psischeme) == 0)
                        std::tie(temppsix[i][j],temppsiy[i][j]) = SuperConsistentScheme(x[j],y[i],rootpsix,rootpsiy,t,dt,T_period,backtrace_scheme);

                    break;
                } // end of case 1

                case 2:{
                    double rootpsix = hermx(phi,psix,psiy,psixy,xadv[i][j], y[i],xo,yo,dx, dy);
                    double rootpsiy = temppsiy[i][j];
                    tempphi[i][j] = hp(phi, psix, psiy, psixy, xadv[i][j], y[i], xo, yo, dx, dy) - dt * Vely(x[j], y[i], t * dt, T_period) * rootpsiy;
                    temppsix[i][j] = Heuns_X(x[j],y[i],rootpsix,rootpsiy,temppsixy[i][j],temppsix[i][j],t,dt,T_period);
                    break;
                }   // end of case 2

                case 3:{
                    double rootpsix = temppsix[i][j];
                    double rootpsiy = hermy(phi,psix,psiy,psixy,x[j],yadv[i][j],xo,yo,dx, dy);
                    tempphi[i][j] = hp(phi, psix, psiy, psixy, x[j],yadv[i][j], xo, yo, dx, dy) - dt * Velx(x[j], y[i], t * dt, T_period) * rootpsix;
                    temppsiy[i][j] = Heuns_Y(x[j],y[i],rootpsix,rootpsiy,temppsixy[i][j],temppsiy[i][j],t,dt,T_period);
                    break;
                }   // end of case 3

                case 4:{
                    tempphi[i][j] = tempphi[i][j] - dt * (Velx(x[j], y[i], t * dt, T_period) * temppsix[i][j] + Vely(x[j], y[i], t * dt, T_period) * temppsiy[i][j]);
                    break;
                } //end of case4

                default:{break;}
                }   // end of switch loop
            }
        }
    }
	}


void update_mixed_derivatives(gridarray &temppsix, gridarray &temppsiy, gridarray &temppsixy,
				unsigned int nx, unsigned int ny, double dx, double dy)
{
double dxy1, dxy2;
double d1, d2, d3, d4;  //Temporary variable to compute mixed derivatives
# pragma omp parallel
{
# pragma omp for schedule(static) collapse(2) 
for(unsigned int k = 0; k < ny; k++){
    for(unsigned int l = 0; l < nx; l++){
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
}
}  // end of function


void update_mixed_derivatives_mpi(gridarray &temppsix, gridarray &temppsiy, gridarray &temppsixy,
				unsigned int nx, unsigned int ny, double dx, double dy, unsigned int ny_Node_min, unsigned int ny_Node_max)
{
double dxy1, dxy2;
double d1, d2, d3, d4;  //Temporary variable to compute mixed derivatives

for(unsigned int k = 0; k < ny; k++){
    for(unsigned int l = ny_Node_min; l < ny_Node_max; l++){
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
}  // end of function






}
#endif
