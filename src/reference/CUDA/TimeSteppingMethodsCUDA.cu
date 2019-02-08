/*
All the time integration methods required for the numerical algorithm are implemented in this header file.
*/

__device__ void Heuns_internal(double xnode, double ynode, double xadv, double yadv, double rootpsix, double rootpsiy,
        unsigned int t, double dt, double T_period, double *temppsix, double *temppsiy, unsigned long indexToWrite)
{
    // The following derivatives are the velocity gradients required for the third equation in the mehtod of characteristics for GRADIENTS
    double graduxnode = gradUx(xnode, ynode, (t + 1) * dt, T_period);
    double graduynode = gradUy(xnode, ynode, (t + 1) * dt, T_period);
    double gradvxnode = gradVx(xnode, ynode, (t + 1) * dt, T_period);
    double gradvynode = gradVy(xnode, ynode, (t + 1) * dt, T_period);
    
    double graduxroot = gradUx(xadv, yadv, t * dt, T_period);
    double graduyroot = gradUy(xadv, yadv, t * dt, T_period);
    double gradvxroot = gradVx(xadv, yadv, t * dt, T_period);
    double gradvyroot = gradVy(xadv, yadv, t * dt, T_period);
    
    // Heun's RK-2 method - first step - for x & y
    temppsix[indexToWrite] = rootpsix - dt * (graduxroot * rootpsix + gradvxroot * rootpsiy);
    temppsiy[indexToWrite] = rootpsiy - dt * (graduyroot * rootpsix + gradvyroot * rootpsiy);
    
    // Heun's RK-2 method - second step - for x & y
    temppsix[indexToWrite] = rootpsix - 0.5 * dt * (graduxroot * rootpsix + gradvxroot * rootpsiy + graduxnode * temppsix[indexToWrite] + gradvxnode * temppsiy[indexToWrite]);
    temppsiy[indexToWrite] = rootpsiy - 0.5 * dt * (graduyroot * rootpsix + gradvyroot * rootpsiy + graduynode * temppsix[indexToWrite] + gradvynode * temppsiy[indexToWrite]);
}

__device__ double Heuns_X(double xnode, double ynode, double rootpsix, double rootpsiy, double temppsixy, double temppsix,
        unsigned int t, double dt, double T_period )
{
    // The following derivatives are the velocity gradients required for the third equation in the mehtod of characteristics for GRADIENTS
    double graduxnode = gradUx(xnode, ynode, t * dt, T_period);
    double gradvxnode = gradVx(xnode, ynode, t * dt, T_period);
    
    double graduxnew = gradUx(xnode, ynode, (t+1) * dt, T_period);
    double gradvxnew = gradVx(xnode, ynode, (t+1) * dt, T_period);
    
    double temp = rootpsix - dt * (gradvxnode * rootpsiy + Vely(xnode, ynode, t * dt, T_period) * temppsixy + graduxnode * temppsix);
    
    //HEUNS
    temp = rootpsix - (dt/2.0) * (gradvxnode * rootpsiy + Vely(xnode, ynode, t * dt, T_period) * temppsixy + graduxnode * temppsix + gradvxnew * rootpsiy + Vely(xnode, ynode, (t+1) * dt, T_period) * temppsixy + graduxnew * temp);
    return temp;
}

__device__ double Heuns_Y(double xnode, double ynode, double rootpsix, double rootpsiy, double temppsixy, double temppsiy,
        unsigned int t, double dt, double T_period)
{
    // The following derivatives are the velocity gradients required for the third equation in the mehtod of characteristics for GRADIENTS
    double graduynode = gradUy(xnode, ynode, t * dt, T_period);
    double gradvynode = gradVy(xnode, ynode, t * dt, T_period);
    
    double graduynew = gradUy(xnode, ynode, (t+1) * dt, T_period);
    double gradvynew = gradVy(xnode, ynode, (t+1) * dt, T_period);
    
    double temp = rootpsiy - dt * (graduynode * rootpsix + Velx(xnode, ynode, t * dt, T_period) * temppsixy + gradvynode * temppsiy);
    
    //HEUNS
    temp = rootpsiy - (dt/2.0) * (graduynode * rootpsix + Velx(xnode, ynode, t * dt, T_period) * temppsixy + gradvynode * temppsiy + graduynew * rootpsix + Velx(xnode, ynode, (t+1) * dt, T_period) * temppsixy + gradvynew * temp);
    return temp;
}

__device__ void SuperConsistentScheme(double xnode, double ynode, double rootpsix, double rootpsiy,
        unsigned int t, double dt, double T_period, char backtrace_scheme[],double *temppsix, double *temppsiy, unsigned long indexToWrite)
{
    
    double graduxnode = gradUx(xnode, ynode, (t + 1) * dt, T_period);
    double graduynode = gradUy(xnode, ynode, (t + 1) * dt, T_period);
    double gradvxnode = gradVx(xnode, ynode, (t + 1) * dt, T_period);
    double gradvynode = gradVy(xnode, ynode, (t + 1) * dt, T_period);
    
    double gradxr1x = 1.0 - dt * graduxnode;
    double gradxr1y = 0.0 - dt * graduynode;
    double gradyr1x = 0.0 - dt * gradvxnode;
    double gradyr1y = 1.0 - dt * gradvynode;
    
    //Commenting out the options of backtrace scheme
    //RK3 is being used now irrespective of specification given in
    //GALS_Advection.cu
    //The if/else functionality may be added later
    //if(strcmp("RK3",backtrace_scheme) == 0){
    
    double ux = Velx(xnode, ynode, (t + 1) * dt, T_period);
    double vy = Vely(xnode, ynode, (t + 1) * dt, T_period);
    double xr = xnode - ux * dt;
    double yr = ynode - vy * dt;
    
    double gradux1 = gradUx(xr, yr, t * dt, T_period);
    double graduy1 = gradUy(xr, yr, t * dt, T_period);
    double gradvx1 = gradVx(xr, yr, t * dt, T_period);
    double gradvy1 = gradVy(xr, yr, t * dt, T_period);
    
    double ux1 = Velx(xr, yr, t * dt, T_period);
    double vy1 = Vely(xr, yr, t * dt, T_period);
    
    // 2 Step
    xr = xnode - ((dt)/2.0 * (0.5 * ux + 0.5 * ux1));
    yr = ynode - ((dt)/2.0 * (0.5 * vy + 0.5 * vy1));
    
    double gradux2 = gradUx(xr, yr, (t * dt + dt/2.0), T_period);
    double graduy2 = gradUy(xr, yr, (t * dt + dt/2.0), T_period);
    double gradvx2 = gradVx(xr, yr, (t * dt + dt/2.0), T_period);
    double gradvy2 = gradVy(xr, yr, (t * dt + dt/2.0), T_period);
    
    double gradxr2x = 1.0 - dt/4.0 * (graduxnode + gradxr1x * gradux1 + gradyr1x * graduy1);
    double gradxr2y = 0.0 - dt/4.0 * (graduynode + gradxr1y * gradux1 + gradyr1y * graduy1);
    double gradyr2x = 0.0 - dt/4.0 * (gradvxnode + gradxr1x * gradvx1 + gradyr1x * gradvy1);
    double gradyr2y = 1.0 - dt/4.0 * (gradvynode + gradxr1y * gradvx1 + gradyr1y * gradvy1);
    
    double gradxrx = 1.0 - dt/6.0 * (graduxnode + gradxr1x * gradux1 + gradyr1x * graduy1 + 4 * (gradxr2x * gradux2 + gradyr2x * graduy2));
    double gradxry = 0.0 - dt/6.0 * (graduynode + gradxr1y * gradux1 + gradyr1y * graduy1 + 4 * (gradxr2y * gradux2 + gradyr2y * graduy2));
    double gradyrx = 0.0 - dt/6.0 * (gradvxnode + gradxr1x * gradvx1 + gradyr1x * gradvy1 + 4 * (gradxr2x * gradvx2 + gradyr2x * gradvy2));
    double gradyry = 1.0 - dt/6.0 * (gradvynode + gradxr1y * gradvx1 + gradyr1y * gradvy1 + 4 * (gradxr2y * gradvx2 + gradyr2y * gradvy2));
    
    // Super consistent method -> psi(n+1) = grad(xroot,yroot) * psiroot(n)
    temppsix[indexToWrite] = (gradxrx * rootpsix + gradyrx * rootpsiy);
    temppsiy[indexToWrite] = (gradxry * rootpsix + gradyry * rootpsiy);
    //}
    /*
     * /else {
     * /    double gradxrx = gradxr1x;
     * /    double gradxry = gradxr1y;
     * /    double gradyrx = gradyr1x;
     * /    double gradyry = gradyr1y;
     *
     * /    // Super consistent method -> psi(n+1) = grad(xroot,yroot) * psiroot(n)
     * /    temppsix[indexToWrite] = (gradxrx * rootpsix + gradyrx * rootpsiy);
     * /    temppsiy[indexToWrite] = (gradxry * rootpsix + gradyry * rootpsiy);
     * /}  // end of if(strcmp("RK3",backtrace_scheme)
     */
}
