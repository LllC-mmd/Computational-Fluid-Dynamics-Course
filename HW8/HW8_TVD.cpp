#include <math.h>
#include <fstream>
#include <iostream>
#include <string>
#include <array>
#include <tuple>
#include <algorithm>
#include "HW8_TVD.hpp"


using namespace std;

int id_map(int id, int id_start, int id_end){
    int id_new;

    if (id < id_start){
        id_new = id_start;
    }
    else if (id > id_end){
        id_new = id_end;
    }
    else {
        id_new = id;
    }
    
    return id_new;
}


Var operator+(const Var& a, const Var& b){
    Var var_sum;

    var_sum.v1 = a.v1 + b.v1;
    var_sum.v2 = a.v2 + b.v2;
    var_sum.v3 = a.v3 + b.v3;

    return var_sum;
}

Var operator-(const Var& a, const Var& b){
    Var var_minus;

    var_minus.v1 = a.v1 - b.v1;
    var_minus.v2 = a.v2 - b.v2;
    var_minus.v3 = a.v3 - b.v3;

    return var_minus;
}

Var operator*(double alpha, const Var&a){
    Var var_mult;

    var_mult.v1 = alpha * a.v1;
    var_mult.v2 = alpha * a.v2;
    var_mult.v3 = alpha * a.v3;

    return var_mult;
}

double get_p(const Var& u){
    double gamma = 1.4;
    double p = (gamma - 1.0) * (u.v3 - 0.5*u.v2*u.v2/u.v1);
    return p;
}

Var get_F(const Var& u){
    Var F;

    double p = get_p(u);
    F.v1 = u.v2;
    F.v2 = u.v2*u.v2/u.v1 + p;
    F.v3 = u.v2*(u.v3+p)/u.v1;

    return F;
}

double get_H(const Var& u){
    double gamma = 1.4;
    double p = get_p(u);
    double a = sqrt(gamma*p/u.v1);
    double H = a*a/(gamma-1.0) + 0.5*u.v2*u.v2/(u.v1*u.v1);
    return H;
}

tuple<double, double, double, double> RoeAvg(const Var& A, const Var& B){
    double rho_roe, u_roe, H_roe, a_roe;
    double Ha, Hb;
    const double gamma = 1.4;

    Ha = get_H(A);
    Hb = get_H(B);
    rho_roe = sqrt(A.v1*B.v1); 
    u_roe = (A.v2/sqrt(A.v1)+B.v2/sqrt(B.v1)) / (sqrt(A.v1)+sqrt(B.v1)); 
    H_roe = (sqrt(A.v1)*Ha+sqrt(B.v1)*Hb) / (sqrt(A.v1)+sqrt(B.v1));
    a_roe = sqrt((gamma-1.0)*(H_roe-0.5*u_roe*u_roe));
    return make_tuple(rho_roe, u_roe, H_roe, a_roe);
}


// For p_ini, p_func, p_gradfunc
// Ref to Toro's book: 4.3 Numerical Solution for Pressure and 4.4 The Complete Solution
double p_ini(Var UL, Var UR, double tol){
    const double gamma = 1.4;
    double uL = UL.v2 / UL.v1;
    double pL = get_p(UL);
    double aL = sqrt(gamma*pL/UL.v1);
    double AL = 2.0/(gamma+1.0)/UL.v1;
    double BL = (gamma-1.0)/(gamma+1.0)*pL;

    double uR = UR.v2 / UR.v1;
    double pR = get_p(UR);
    double aR = sqrt(gamma*pR/UR.v1);
    double AR = 2.0/(gamma+1.0)/UR.v1;
    double BR = (gamma-1.0)/(gamma+1.0)*pR;

    double pmax = max(pL, pR);
    double pmin = min(pL, pR);
    double Qmax = pmax/pmin;

    double p_pv = max(tol, 0.5*(pL+pR) - 0.125*(uR-uL)*(UL.v1+UR.v1)*(aL+aR));
    if (Qmax<=2.0 && (pmin<=p_pv && p_pv<=pmax)){
        return p_pv;
    }
    else {
        if (p_pv<=pmin){
            double p_TR = pow((aL+aR-0.5*(gamma-1.0)*(uR-uL))/(aL/pow(pL, (gamma-1.0)/2.0/gamma)+aR/pow(pR, (gamma-1.0)/2.0/gamma)), 2.0*gamma/(gamma-1.0));
            return p_TR;
        }
        else {
            double p_TS = (sqrt(AL/(p_pv+BL))*pL+sqrt(AR/(p_pv+BR))*pR-(uR-uL))/(sqrt(AL/(p_pv+BL))+sqrt(AR/(p_pv+BR)));
            return p_TS;
        }
    }
}

double p_func(Var UK, double p){
    double fK;

    const double gamma = 1.4;
    double uL = UK.v2 / UK.v1;
    double pK = get_p(UK);
    
    double AK = 2.0/(gamma+1.0)/UK.v1;
    double BK = (gamma-1.0)/(gamma+1.0)*pK;
    double aK = sqrt(gamma*pK/UK.v1);

    if (p>pK){
        fK = (p-pK)*sqrt(AK/(p+BK));
    }
    else{
        fK = 2*aK/(gamma-1.0)*(pow(p/pK, (gamma-1.0)/2/gamma)-1.0);
    }

    return fK;
}

double p_gradfunc(Var UK, double p){
    double gradK;

    const double gamma = 1.4;
    double pK = get_p(UK);
    double AK = 2.0/(gamma+1.0)/UK.v1;
    double BK = (gamma-1.0)/(gamma+1.0)*pK;
    double aK = sqrt(gamma*pK/UK.v1);
   
    if (p>pK){
        gradK = sqrt(AK/(p+BK))*(1.0-(p-pK)/2.0/(BK+p));
    }
    else{
        gradK = 1.0/UK.v1/aK*pow(p/pK, -(gamma+1.0)/2.0/gamma);
    }

    return gradK;
}


Var RP_solver(Var UL, Var UR, double x, double t, double tol, double maxiter){
    Var RP_sol;
    double p_star, u_star;
    double p_tmp, grad_p;
    int c;
    double SL, SR;
    double SL_tail, SL_head, SR_tail, SR_head;
    bool Lshock, Rshock;
    const double gamma = 1.4;

    double simVar = x/t;
    double uL = UL.v2 / UL.v1;
    double pL = get_p(UL);
    double aL = sqrt(gamma*pL/UL.v1);
    double uR = UR.v2 / UR.v1;
    double pR = get_p(UR);
    double aR = sqrt(gamma*pR/UR.v1);

    // Middle Region 
    // ------Solve p_star
    p_tmp = p_ini(UL, UR, tol);
    c = 0;
    while (true){
        grad_p = p_gradfunc(UL, p_tmp) + p_gradfunc(UR, p_tmp);
        p_star = p_tmp - (p_func(UL, p_tmp)+p_func(UR, p_tmp)+uR-uL)/grad_p;
        if (fabs(2.0*(p_star-p_tmp)/(p_star+p_tmp)) < tol){
            break;
        }
        else if (c>=maxiter){
            cout << "Max iteration exceed, last update: " << fabs(2.0*(p_star-p_tmp)/(p_star+p_tmp)) << endl;
            break;
        }
        else{
            p_tmp = p_star;
        }
        c += 1;
    }
    // ------Solve u*
    u_star = 0.5*(uL+uR) + 0.5*(p_func(UR, p_star)-p_func(UL, p_star));
    
    // Left Region
    if (p_star>pL){
        // ------Left Shock Wave
        SL = uL - aL*(sqrt((gamma+1.0)/2.0/gamma*p_star/pL+(gamma-1.0)/2.0/gamma));
        Lshock = true;
    }
    else{
        // ------Left Rarefaction Wave
        SL_head = uL - aL;
        SL_tail = u_star - aL*pow(p_star/pL, (gamma-1.0)/2.0/gamma);
        Lshock = false;
    }

    // Right Region
    if (p_star>pR){
        // ------Right Shock Wave
        SR = uR + aR*(sqrt((gamma+1.0)/2.0/gamma*p_star/pR+(gamma-1.0)/2.0/gamma));
        Rshock = true;
    }
    else{
        // ------Right Rarefaction Wave
        SR_head = uR + aR;
        SR_tail = u_star + aR*pow(p_star/pR, (gamma-1.0)/2.0/gamma);
        Rshock = false;
    }

    if (Lshock && Rshock){
        if (simVar < SL){
            RP_sol = UL;
        }
        else if (simVar < u_star){
            RP_sol.v1 = UL.v1*((p_star/pL+(gamma-1.0)/(gamma+1.0))/((gamma-1.0)/(gamma+1.0)*p_star/pL+1.0));
            RP_sol.v2 = RP_sol.v1*u_star;
            RP_sol.v3 = p_star/(gamma-1.0) + 0.5*RP_sol.v1*(u_star*u_star);
        }
        else if (simVar < SR){
            RP_sol.v1 = UR.v1*((p_star/pR+(gamma-1.0)/(gamma+1.0))/((gamma-1.0)/(gamma+1.0)*p_star/pR+1.0));
            RP_sol.v2 = RP_sol.v1*u_star;
            RP_sol.v3 = p_star/(gamma-1.0) + 0.5*RP_sol.v1*(u_star*u_star);
        }
        else {
            RP_sol = UR;
        }
    }
    else if (Lshock && (!Rshock)){
        if (simVar < SL){
            RP_sol = UL;
        }
        else if (simVar < u_star){
            RP_sol.v1 = UL.v1*((p_star/pL+(gamma-1.0)/(gamma+1.0))/((gamma-1.0)/(gamma+1.0)*p_star/pL+1.0));
            RP_sol.v2 = RP_sol.v1*u_star;
            RP_sol.v3 = p_star/(gamma-1.0) + 0.5*RP_sol.v1*(u_star*u_star);
        }
        else if (simVar < SR_tail){
            RP_sol.v1 = UR.v1*pow(p_star/pR, 1.0/gamma);
            RP_sol.v2 = RP_sol.v1*u_star;
            RP_sol.v3 = p_star/(gamma-1.0) + 0.5*RP_sol.v1*(u_star*u_star);
        }
        else if (simVar < SR_head){
            RP_sol.v1 = UR.v1*pow((2.0/(gamma+1.0)-(gamma-1.0)/(gamma+1.0)/aR*(uR-simVar)), 2.0/(gamma-1.0));
            RP_sol.v2 = RP_sol.v1*2.0/(gamma+1.0)*(-aR + (gamma-1.0)/2.0*uR + simVar);
            RP_sol.v3 = pR*pow(2.0/(gamma+1.0)-(gamma-1.0)/(gamma+1.0)/aR*(uR-simVar), 2*gamma/(gamma-1))/(gamma-1.0)
                        + 0.5*RP_sol.v2*RP_sol.v2/RP_sol.v1;
        }
        else {
            RP_sol = UR;
        }
    }
    else if ((!Lshock) && Rshock){
        if (simVar < SL_head){
            RP_sol = UL;
        }
        else if (simVar < SL_tail){
            RP_sol.v1 = UL.v1*pow((2.0/(gamma+1.0)+(gamma-1.0)/(gamma+1.0)/aL*(uL-simVar)), 2.0/(gamma-1.0));
            RP_sol.v2 = RP_sol.v1*2.0/(gamma+1.0)*(aL + (gamma-1.0)/2.0*uL + simVar);
            RP_sol.v3 = pL*pow(2.0/(gamma+1.0)+(gamma-1.0)/(gamma+1.0)/aL*(uL-simVar), 2*gamma/(gamma-1))/(gamma-1.0)
                        + 0.5*RP_sol.v2*RP_sol.v2/RP_sol.v1;
        }
        else if (simVar < u_star){
            RP_sol.v1 = UL.v1*pow(p_star/pL, 1.0/gamma);
            RP_sol.v2 = RP_sol.v1*u_star;
            RP_sol.v3 = p_star/(gamma-1.0) + 0.5*RP_sol.v1*(u_star*u_star);
        }
        else if (simVar < SR){
            RP_sol.v1 = UR.v1*((p_star/pR+(gamma-1.0)/(gamma+1.0))/((gamma-1.0)/(gamma+1.0)*p_star/pR+1.0));
            RP_sol.v2 = RP_sol.v1*u_star;
            RP_sol.v3 = p_star/(gamma-1.0) + 0.5*RP_sol.v1*(u_star*u_star);
        }
        else {
            RP_sol = UR;
        }
    }
    else {
        if (simVar < SL_head){
            RP_sol = UL;
        }
        else if (simVar < SL_tail){
            RP_sol.v1 = UL.v1*pow((2.0/(gamma+1.0)+(gamma-1.0)/(gamma+1.0)/aL*(uL-simVar)), 2.0/(gamma-1.0));
            RP_sol.v2 = RP_sol.v1*2.0/(gamma+1.0)*(aL + (gamma-1.0)/2.0*uL + simVar);
            RP_sol.v3 = pL*pow(2.0/(gamma+1.0)+(gamma-1.0)/(gamma+1.0)/aL*(uL-simVar), 2*gamma/(gamma-1))/(gamma-1.0)
                        + 0.5*RP_sol.v2*RP_sol.v2/RP_sol.v1;
        }
        else if (simVar < u_star){
            RP_sol.v1 = UL.v1*pow(p_star/pL, 1.0/gamma);
            RP_sol.v2 = RP_sol.v1*u_star;
            RP_sol.v3 = p_star/(gamma-1.0) + 0.5*RP_sol.v1*(u_star*u_star);
        }
        else if (simVar < SR_tail){
            RP_sol.v1 = UR.v1*pow(p_star/pR, 1.0/gamma);
            RP_sol.v2 = RP_sol.v1*u_star;
            RP_sol.v3 = p_star/(gamma-1.0) + 0.5*RP_sol.v1*(u_star*u_star);
        }
        else if (simVar < SR_head){
            RP_sol.v1 = UR.v1*pow((2.0/(gamma+1.0)-(gamma-1.0)/(gamma+1.0)/aR*(uR-simVar)), 2.0/(gamma-1.0));
            RP_sol.v2 = RP_sol.v1*2.0/(gamma+1.0)*(-aR + (gamma-1.0)/2.0*uR + simVar);
            RP_sol.v3 = pR*pow(2.0/(gamma+1.0)-(gamma-1.0)/(gamma+1.0)/aR*(uR-simVar), 2*gamma/(gamma-1))/(gamma-1.0)
                        + 0.5*RP_sol.v2*RP_sol.v2/RP_sol.v1;
        }
        else {
            RP_sol = UR;
        }
    }
    return RP_sol;
}


double upwindSlope(Var U_ip1, Var U_i, double rho_roe, double u_roe, double a_roe, int index){
    double drho, dp, du;
    double alpha;

    if (index==1){
        dp = get_p(U_ip1)-get_p(U_i);
        du = U_ip1.v2/U_ip1.v1 - U_i.v2/U_i.v1;
        alpha = 1.0/(2.0*a_roe*a_roe)*(dp - rho_roe*a_roe*du);
    }
    else if (index==2){
        drho = U_ip1.v1 - U_i.v1;
        dp = get_p(U_ip1)-get_p(U_i);
        alpha = drho - dp/(a_roe*a_roe);        
    }
    else{
        dp = get_p(U_ip1)-get_p(U_i);
        du = U_ip1.v2/U_ip1.v1 - U_i.v2/U_i.v1;
        alpha = 1.0/(2.0*a_roe*a_roe)*(dp + rho_roe*a_roe*du);
    }
    return alpha;
}

double upwindSlopeRatio(double slope, double slopeUpwind){
    double theta;
    const double epsilon = 0.001;

    if (slope*slope + slopeUpwind*slopeUpwind == 0){
        theta = 0;
    }
    else if (slope==0){
        theta = fabs(slopeUpwind)/epsilon;
    }
    else{
        theta = slopeUpwind/slope;
    }
    return theta;
}

double minimodLimiter(double theta){
    double phi;

    if (theta<=0){
        phi = 0;
    }
    else if (theta<=1){
        phi = theta;
    }
    else{
        phi = 1;
    }
    return phi;
}

double SuperBeeLimiter(double theta){
    double phi;

    if (theta<=0){
        phi = 0;
    }
    else if (theta<=0.5){
        phi = 2*theta;
    }
    else if (theta<=1.0){
        phi = 1;
    }
    else if (theta<=2.0){
        phi = theta;
    }
    else{
        phi = 2.0;
    }
    return phi;
}


// 2nd upwind scheme with TVD property
void fluxCalc_TVD(Vector<Var>& sol, Vector<Var>& FluxX, double dx, double dt, int nx, const string fluxLimiter){
    // for c >= 0, im1, i0, ip1 are chosen as stencil points
    // for c < 0, i0, ip1, ip2 are chosen as stencil points
    int im1, i0; 
    int iUpwind_1, iUpwind_2, iUpwind_3;
    double rho_roe, u_roe, H_roe, a_roe;
    double lambda_1, lambda_2, lambda_3;
    double c1, c2, c3;
    double drho, dp, du;
    double alpha_1, alpha_2, alpha_3;
    double alphaUpw_1, alphaUpw_2, alphaUpw_3;
    double theta_1, theta_2, theta_3;
    double phi_1, phi_2, phi_3;

    Var r1, r2, r3;
    r1.v1 = 1.0;
    r2.v1 = 1.0;
    r3.v1 = 1.0;
    
    for (int i = 0; i < nx+1; i++){
        im1 = id_map(i-1, 0, nx-1);
        i0 = id_map(i, 0, nx-1);

        // Roe averaging related calculation
        tie(rho_roe, u_roe, H_roe, a_roe) = RoeAvg(sol[im1], sol[i0]);
        lambda_1 = u_roe - a_roe;
        lambda_2 = u_roe;
        lambda_3 = u_roe + a_roe;
        c1 = lambda_1*dt/dx;
        c2 = lambda_2*dt/dx;
        c3 = lambda_3*dt/dx;
        drho = sol[i0].v1 - sol[im1].v1;
        dp = get_p(sol[i0])-get_p(sol[im1]);
        du = sol[i0].v2/sol[i0].v1 - sol[im1].v2/sol[im1].v1;
        alpha_1 = 1.0/(2.0*a_roe*a_roe)*(dp - rho_roe*a_roe*du);
        alpha_2 = drho - dp/(a_roe*a_roe);
        alpha_3 = 1.0/(2.0*a_roe*a_roe)*(dp + rho_roe*a_roe*du);

        r1.v2 = u_roe - a_roe;
        r1.v3 = H_roe - u_roe*a_roe;
        r2.v2 = u_roe;
        r2.v3 = u_roe*u_roe/2.0;
        r3.v2 = u_roe + a_roe;
        r3.v3 = H_roe + u_roe*a_roe;

        // Determine the upwind direction and its slope
        if (c1>=0){
            iUpwind_1 = id_map(i-2, 0, nx-1);
            alphaUpw_1 = upwindSlope(sol[im1], sol[iUpwind_1], rho_roe, u_roe, a_roe, 1);
        }
        else{
            iUpwind_1 = id_map(i+1, 0, nx-1);
            alphaUpw_1 = upwindSlope(sol[iUpwind_1], sol[i0], rho_roe, u_roe, a_roe, 1);
        }
        if (c2>=0){
            iUpwind_2 = id_map(i-2, 0, nx-1);
            alphaUpw_2 = upwindSlope(sol[im1], sol[iUpwind_2], rho_roe, u_roe, a_roe, 2);
        }
        else{
            iUpwind_2 = id_map(i+1, 0, nx-1);
            alphaUpw_2 = upwindSlope(sol[iUpwind_2], sol[i0], rho_roe, u_roe, a_roe, 2);
        }
        if (c3>=0){
            iUpwind_3 = id_map(i-2, 0, nx-1);
            alphaUpw_3 = upwindSlope(sol[im1], sol[iUpwind_3], rho_roe, u_roe, a_roe, 3);
        }
        else{
            iUpwind_3 = id_map(i+1, 0, nx-1);
            alphaUpw_3 = upwindSlope(sol[iUpwind_3], sol[i0], rho_roe, u_roe, a_roe, 3);
        }

        theta_1 = upwindSlopeRatio(alpha_1, alphaUpw_1);
        theta_2 = upwindSlopeRatio(alpha_2, alphaUpw_2);
        theta_3 = upwindSlopeRatio(alpha_3, alphaUpw_3);
        
        // Flux limiting
        if (fluxLimiter=="Minimod"){
            phi_1 = minimodLimiter(theta_1);
            phi_2 = minimodLimiter(theta_2);
            phi_3 = minimodLimiter(theta_3);
        }
        else if (fluxLimiter=="SuperBee"){
            phi_1 = SuperBeeLimiter(theta_1);
            phi_2 = SuperBeeLimiter(theta_2);
            phi_3 = SuperBeeLimiter(theta_3);
        }
        else{
            cout << "Unknown Limter" << endl;
        }

        FluxX[i] =  0.5*(get_F(sol[im1])+get_F(sol[i0])) 
                    - 0.5*(fabs(lambda_1)-fabs(lambda_1)*(1-c1)*phi_1)*alpha_1*r1
                    - 0.5*(fabs(lambda_2)-fabs(lambda_2)*(1-c2)*phi_2)*alpha_2*r2
                    - 0.5*(fabs(lambda_3)-fabs(lambda_3)*(1-c3)*phi_3)*alpha_3*r3;
    }
}



void step(Vector<Var>& sol, Vector<Var>& FluxX, double dx, double dt, int nx, const string fluxLimiter){
    Vector<Var> sol_tmp(nx);

    fluxCalc_TVD(sol, FluxX, dx, dt, nx, fluxLimiter);
    for (int i=0; i<nx; i++){
        sol_tmp[i] = sol[i] - dt/dx * (FluxX[i+1] - FluxX[i]);
    }

    fluxCalc_TVD(sol_tmp, FluxX, dx, dt, nx, fluxLimiter);
    for (int i=0; i<nx; i++){
        sol_tmp[i] = 3.0/4.0*sol[i] + 1.0/4.0*(sol_tmp[i] - dt/dx*(FluxX[i+1]-FluxX[i]));
    }

    fluxCalc_TVD(sol_tmp, FluxX, dx, dt, nx, fluxLimiter);
    for (int i=0; i<nx; i++){
        sol[i] = 1.0/3.0*sol[i] + 2.0/3.0*(sol_tmp[i] - dt/dx*(FluxX[i+1]-FluxX[i]));
    }
}


/* g++ -Og -o HW8_TVD HW8_TVD.cpp HW8_TVD.hpp */
/* ./HW8_TVD 0.2 200 0.001 -0.5 0.5 0.0 0.01 SuperBee 2ndUpwind */
/* ./HW8_TVD 0.2 200 0.001 -0.5 0.5 0.0 0.01 Minimod 2ndUpwind */
/* ./HW8_TVD 0.2 200 0.001 -0.5 0.5 0.0 0.01 SuperBee ExactRP */
/* ./HW8_TVD 0.2 200 0.001 -0.5 0.5 0.0 0.01 SuperBee 2ndUpwind */
/* ./HW8_TVD 0.2 400 0.0005 -0.5 0.5 0.0 0.005 SuperBee 2ndUpwind */
/* ./HW8_TVD 0.2 400 0.0005 -0.5 0.5 0.0 0.005 Minimod 2ndUpwind */
/* ./HW8_TVD 0.2 400 0.0005 -0.5 0.5 0.0 0.005 SuperBee ExactRP */
/* ./HW8_TVD 0.2 1000 0.0002 -0.5 0.5 0.0 0.002 SuperBee 2ndUpwind */
/* ./HW8_TVD 0.2 1000 0.0002 -0.5 0.5 0.0 0.002 Minimod 2ndUpwind */
/* ./HW8_TVD 0.2 1000 0.0002 -0.5 0.5 0.0 0.002 SuperBee ExactRP */
/* ./HW8_TVD 0.2 2000 0.0001 -0.5 0.5 0.0 0.001 SuperBee 2ndUpwind */
/* ./HW8_TVD 0.2 2000 0.0001 -0.5 0.5 0.0 0.001 Minimod 2ndUpwind */
/* ./HW8_TVD 0.2 2000 0.0001 -0.5 0.5 0.0 0.001 SuperBee ExactRP */
int main(int argc, char *argv[]) {
    double elapsedTime; 
    int stepCounter, reportStep;
    double Xstart, Xend, x0;
    double dx, dt, x_i;
    double px;
    double rhoL, uL, pL;
    double rhoR, uR, pR;
    Var UL;
    Var UR;
    ofstream outFile;

    const double gamma = 1.4;
    const double tol = 1e-6;
    const double maxiter = 1e5;
    
    elapsedTime =  atof(argv[1]);
    reportStep = atoi(argv[2]);
    dt = atof(argv[3]);
    Xstart = atof(argv[4]);
    Xend = atof(argv[5]);
    x0 = atof(argv[6]);
    dx = stod(argv[7]);
    const string fluxLimter = argv[8];
    const string solverScheme = argv[9];

    cout << "dx: " << dx << "\tdt: " << dt << endl;
    cout << "ReportTime Interval: " << reportStep*dt << endl;
    cout << "Solution Domain: " << "[" << Xstart << ", " << Xend << "]" << endl;

    int nx = int((Xend - Xstart) / dx);

    Vector<Var> U(nx);
    // place Ui in the middle of F_i and F_i+1
    Vector<Var> FluxX(nx+1);

    rhoL = 1.0;
    uL = 0.0;
    pL = 1.0;
    rhoR = 0.125;
    uR = 0.0;
    pR = 0.1;

    UL.v1 = rhoL;
    UL.v2 = rhoL*uL;
    UL.v3 = pL/(gamma-1.0) + 0.5*rhoL*(uL*uL);
    UR.v1 = rhoR;
    UR.v2 = rhoR*uR;
    UR.v3 = pR/(gamma-1.0) + 0.5*rhoR*(uR*uR);

    // Initialization
    for (int i = 0; i < nx; i++){
        px = Xstart + (i+0.5) * dx - x0;
        if (px < 0){
            U[i] = UL;
        }
        else {
            U[i] = UR;
        }
    }

    stepCounter = 0;
    // Exact Solution of Riemann Problem
    if (solverScheme=="ExactRP"){
        for (double t = 0; t < elapsedTime; t+=dt){
            for (int i = 0; i < nx; i++){
                px = Xstart + (i+0.5) * dx - x0;
                // Riemann Problem Solver
                U[i] = RP_solver(UL, UR, px, t, tol, maxiter);
            }
            stepCounter += 1;
            // Report result
            if (stepCounter % reportStep==0){
                // ------report rho
                string name = "Exact_Rho_" + to_string(double(stepCounter*dt)) + ".txt";
                cout << "Time: " << t+dt << ", " << name << " recorded" << endl;
                outFile.open(name, ios::out);
                for (int i = 0; i < nx; i++){
                    outFile << U[i].v1 << " ";
                }
                outFile.close();
                // ------report U
                name = "Exact_U_" + to_string(double(stepCounter*dt)) + ".txt";
                cout << "Time: " << t+dt << ", " << name << " recorded" << endl;
                outFile.open(name, ios::out);
                for (int i = 0; i < nx; i++){
                    outFile << U[i].v2/U[i].v1 << " ";
                }
                outFile.close();
                // ------report P
                name = "Exact_P_" + to_string(double(stepCounter*dt)) + ".txt";
                cout << "Time: " << t+dt << ", " << name << " recorded" << endl;
                outFile.open(name, ios::out);
                for (int i = 0; i < nx; i++){
                    outFile << get_p(U[i]) << " ";
                }
                outFile.close();
            }
        }
    }
    // 2nd upwind solver
    else {
        for (double t = 0; t < elapsedTime; t+=dt){
            step(U, FluxX, dx, dt, nx, fluxLimter);
            stepCounter += 1;
            // Report result
            if (stepCounter % reportStep==0){
                // ------report rho
                string name = fluxLimter +"_Rho_" + to_string(double(stepCounter*dt)) + ".txt";
                cout << "Time: " << t+dt << ", " << name << " recorded" << endl;
                outFile.open(name, ios::out);
                for (int i = 0; i < nx; i++){
                    outFile << U[i].v1 << " ";
                }
                outFile.close();
                // ------report U
                name = fluxLimter + "_U_" + to_string(double(stepCounter*dt)) + ".txt";
                cout << "Time: " << t+dt << ", " << name << " recorded" << endl;
                outFile.open(name, ios::out);
                for (int i = 0; i < nx; i++){
                    outFile << U[i].v2/U[i].v1 << " ";
                }
                outFile.close();
                // ------report P
                name = fluxLimter + "_P_" + to_string(double(stepCounter*dt)) + ".txt";
                cout << "Time: " << t+dt << ", " << name << " recorded" << endl;
                outFile.open(name, ios::out);
                for (int i = 0; i < nx; i++){
                    outFile << get_p(U[i]) << " ";
                }
                outFile.close();
            }
        }
    }

    return 0;
}