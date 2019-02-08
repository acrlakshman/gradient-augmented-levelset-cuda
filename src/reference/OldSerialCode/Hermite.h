//
//  Hermite.h
//  
//
//  Created by Raunak Bardia on 10/6/14.
//
//

#ifndef _Hermite_h
#define _Hermite_h

#include <math.h>
#include <stdio.h>
#include <string.h>

using namespace std;

double basepolynomial(double x, double y, int alphax, int alphay, int vx, int vy, double dx, double dy, double xo, double yo){
    
    double bpx, bpy;
    double etax = (x - xo)/dx;
    double etay = (y - yo)/dy;
    
    switch(2 * alphax + vx + 1){   // Switching between base polynomials for x
        case 1:
            bpx = 1 - 3 * pow(etax,2) + 2 * pow(etax,3);
            break;
        case 2:
            bpx = -2 * pow(etax,3) + 3 * pow(etax,2);
            break;
        case 3:
            bpx = pow(etax,3) - 2 * pow(etax,2) + etax;
            break;
        case 4:
            bpx = pow(etax,3) - pow(etax,2);
            break;
    }

    switch(2 * alphay + vy + 1){   // Switching between base polynomials for y
        case 1:
            bpy = 1 - 3 * pow(etay,2) + 2 * pow(etay,3);
            break;
        case 2:
            bpy = -2 * pow(etay,3) + 3 * pow(etay,2);
            break;
        case 3:
            bpy = pow(etay,3) - 2 * pow(etay,2) + etay;
            break;
        case 4:
            bpy = pow(etay,3) - pow(etay,2);
            break;
    }
    
    double result = bpx * bpy;
    return result;
}

double gradbpx(double x, double y, int alphax, int alphay, int vx, int vy, double dx, double dy, double xo, double yo){
    
    double bpx, bpy;
    double etax = (x - xo)/dx;
    double etay = (y - yo)/dy;
    
    switch(2 * alphax + vx + 1){   // Switching between base polynomials for x
        case 1:
            bpx = - 6 * etax + 6 * pow(etax,2);
            break;
        case 2:
            bpx = -6 * pow(etax,2) + 6 * etax;
            break;
        case 3:
            bpx = 3 * pow(etax,2) - 4 * etax + 1;
            break;
        case 4:
            bpx = 3 * pow(etax,2) - 2 * etax;
            break;
    }
    
    switch(2 * alphay + vy + 1){   // Switching between base polynomials for y
        case 1:
            bpy = 1 - 3 * pow(etay,2) + 2 * pow(etay,3);
            break;
        case 2:
            bpy = -2 * pow(etay,3) + 3 * pow(etay,2);
            break;
        case 3:
            bpy = pow(etay,3) - 2 * pow(etay,2) + etay;
            break;
        case 4:
            bpy = pow(etay,3) - pow(etay,2);
            break;
    }
    
    double result = bpx * bpy;
    
    return result;
}

double gradbpy(double x, double y, int alphax, int alphay, int vx, int vy, double dx, double dy, double xo, double yo){
    
    double bpx, bpy;
    double etax = (x - xo)/dx;
    double etay = (y - yo)/dy;
    
    switch(2 * alphax + vx + 1){   // Switching between base polynomials for x
        case 1:
            bpx = 1 - 3 * pow(etax,2) + 2 * pow(etax,3);
            break;
        case 2:
            bpx = -2 * pow(etax,3) + 3 * pow(etax,2);
            break;
        case 3:
            bpx = pow(etax,3) - 2 * pow(etax,2) + etax;
            break;
        case 4:
            bpx = pow(etax,3) - pow(etax,2);
            break;
    }
    
    switch(2 * alphay + vy + 1){   // Switching between base polynomials for y
        case 1:
            bpy = - 6 * etay + 6 * pow(etay,2);
            break;
        case 2:
            bpy = -6 * pow(etay,2) + 6 * etay;
            break;
        case 3:
            bpy = 3 * pow(etay,2) - 4 * etay + 1;
            break;
        case 4:
            bpy = 3 * pow(etay,2) - 2 * etay;
            break;
    }
    
    double result = bpx * bpy;
    
    return result;
}

double hp(double phi[4], double psix[4], double psiy[4], double psixy[4], double x, double y, double xo, double yo, double dx, double dy){
    
    double H = 0, delta[4],d;
    int alphax, alphay, vx, vy;
    double bp;  //Base Polynomial Section
    for(int i = 0; i < 4; i++){
        switch (i){
            case 0:{
                alphax = 0;
                alphay = 0;
                memcpy(delta, phi, sizeof(delta));
                break;
            }
            case 1:{
                alphax = 0;
                alphay = 1;
                memcpy(delta, psiy, sizeof(delta));
                break;
            }
            case 2:{
                alphax = 1;
                alphay = 0;
                memcpy(delta, psix, sizeof(delta));
                break;
            }
            case 3:{
                alphax = 1;
                alphay = 1;
                memcpy(delta, psixy, sizeof(delta));
                break;
            }
        }
        for(int j = 0; j < 4; j++){
            switch (j){
                case 0:{
                    vx = 0;
                    vy = 0;
                    d = delta[0];
                    break;
                }
                case 1:{
                    vx = 0;
                    vy = 1;
                    d = delta[2];
                    break;
                }
                case 2:{
                    vx = 1;
                    vy = 0;
                    d = delta[1];
                    break;
                }
                case 3:{
                    vx = 1;
                    vy = 1;
                    d = delta[3];
                    break;
                }
            }
            bp = basepolynomial(x, y, alphax, alphay, vx, vy, dx, dy, xo, yo);
            H = H + pow(dx, alphax) * pow(dy, alphay) * bp * d;
        }
    }
    return H;
}

double hermx(double* phi, double* psix, double* psiy, double* psixy, double xadv, double yadv, double xo, double yo, double dx, double dy){
    
    double gradientx = 0, delta[4], d;
    int alphax, alphay, vx, vy;
    double bp;  //Base Polynomial Section
    for(int i = 0; i < 4; i++){
        switch (i){
            case 0:{
                alphax = 0;
                alphay = 0;
                memcpy(delta, phi, sizeof(delta));
                break;
            }
            case 1:{
                alphax = 0;
                alphay = 1;
                memcpy(delta, psiy, sizeof(delta));
                break;
            }
            case 2:{
                alphax = 1;
                alphay = 0;
                memcpy(delta, psix, sizeof(delta));
                break;
            }
            case 3:{
                alphax = 1;
                alphay = 1;
                memcpy(delta, psixy, sizeof(delta));
                break;
            }
        }
        for(int j = 0; j < 4; j++){
            switch (j){
                case 0:{
                    vx = 0;
                    vy = 0;
                    d = delta[0];
                    break;
                }
                case 1:{
                    vx = 0;
                    vy = 1;
                    d = delta[2];
                    break;
                }
                case 2:{
                    vx = 1;
                    vy = 0;
                    d = delta[1];
                    break;
                }
                case 3:{
                    vx = 1;
                    vy = 1;
                    d = delta[3];
                    break;
                }
            }
            bp = gradbpx(xadv, yadv, alphax, alphay, vx, vy, dx, dy, xo, yo);
            gradientx = gradientx + pow(dx, alphax) * pow(dy, alphay) * bp * d * (1/dx);
        }
    }
    return gradientx;
}

double hermy(double* phi, double* psix, double* psiy, double* psixy, double xadv, double yadv, double xo, double yo, double dx, double dy){
    
    double gradienty = 0, d, delta[4];
    int alphax, alphay, vx, vy;
    double bp;  //Base Polynomial Section
    for(int i = 0; i < 4; i++){
        switch (i){
            case 0:{
                alphax = 0;
                alphay = 0;
                memcpy(delta, phi, sizeof(delta));
                break;
            }
            case 1:{
                alphax = 0;
                alphay = 1;
                memcpy(delta, psiy, sizeof(delta));
                break;
            }
            case 2:{
                alphax = 1;
                alphay = 0;
                memcpy(delta, psix, sizeof(delta));
                break;
            }
            case 3:{
                alphax = 1;
                alphay = 1;
                memcpy(delta, psixy, sizeof(delta));
                break;
            }
        }
        for(int j = 0; j < 4; j++){
            switch (j){
                case 0:{
                    vx = 0;
                    vy = 0;
                    d = delta[0];
                    break;
                }
                case 1:{
                    vx = 0;
                    vy = 1;
                    d = delta[2];
                    break;
                }
                case 2:{
                    vx = 1;
                    vy = 0;
                    d = delta[1];
                    break;
                }
                case 3:{
                    vx = 1;
                    vy = 1;
                    d = delta[3];
                    break;
                }
            }
            bp = gradbpy(xadv, yadv, alphax, alphay, vx, vy, dx, dy, xo, yo);
            gradienty = gradienty + pow(dx, alphax) * pow(dy, alphay) * bp * d * (1/dy);
        }
    }
    return gradienty;
}


double hp1D(double phi[2], double psix[2], double x, double xo, double dx){
    
    double H = 0;
    double etax = (x - xo)/dx;
    H = (1 - 3*pow(etax,2) + 2 * pow(etax,3)) * phi[0] + (3 * pow(etax,2) - 2 * pow(etax,3)) * phi[1] + (etax + pow(etax,3) - 2 * pow(etax,2)) * psix[0] + (pow(etax,3) - pow(etax,2)) * psix[1];

    return H;
}

#endif
