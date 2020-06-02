#include <math.h>
#include <fstream>
#include <iostream>
#include <string>
#include <array>
#include <tuple>
#include <map>
#include <algorithm>
#include "HW9_WENO.hpp"


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

Var operator+(const double alpha, const Var& a){
    Var var_sum;

    var_sum.v1 = a.v1 + alpha;
    var_sum.v2 = a.v2 + alpha;
    var_sum.v3 = a.v3 + alpha;

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

Var operator*(const Var&a, const Var&b){
    Var var_mult;

    var_mult.v1 = a.v1 * b.v1;
    var_mult.v2 = a.v2 * b.v2;
    var_mult.v3 = a.v3 * b.v3;

    return var_mult;
}


Var operator/(const Var&a, const Var&b){
    Var var_div;

    var_div.v1 = a.v1 / b.v1;
    var_div.v2 = a.v2 / b.v2;
    var_div.v3 = a.v3 / b.v3;

    return var_div;
}


Var operator/(const double alpha, const Var&b){
    Var var_div;

    var_div.v1 = alpha / b.v1;
    var_div.v2 = alpha / b.v2;
    var_div.v3 = alpha / b.v3;

    return var_div;
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

// smooth indicator calculation
Var get_beta(Var v1, Var v2, Var v3){
    Var beta;
    beta = 13.0/12.0*(v1-2.0*v2+v3)*(v1-2.0*v2+v3) + 1.0/4.0*(3.0*v1-4.0*v2+v3)*(3.0*v1-4.0*v2+v3);
    return beta;
}

// variable reconstruction with 5th WENO scheme
void varRec_WENO(Vector<Var>& sol, Vector<Var>& varL, Vector<Var>& varR, int nx, const map<string, double> coeff_table){
    int im3, im2, im1, i0, ip1, ip2;
    Var v0, v1, v2;
    Var w0, w1, w2;
    Var alpha0, alpha1, alpha2;
    Var beta0, beta1, beta2;
    
    const double epsilon = 1e-6;

    for (int i=0; i<nx+1; i++){
        im3 = id_map(i-3, 0, nx-1);
        im2 = id_map(i-2, 0, nx-1);
        im1 = id_map(i-1, 0, nx-1);
        i0 = id_map(i, 0, nx-1);
        ip1 = id_map(i+1, 0, nx-1);
        ip2 = id_map(i+2, 0, nx-1);
        // For the variable on the left side of the i-th edge, we have:
        // 1) stencil for r=2: [i-3, i-2, i-1]
        v2 = coeff_table.at("c20")*sol[im3] + coeff_table.at("c21")*sol[im2] + coeff_table.at("c22")*sol[im1];
        // 2) stencil for r=1: [i-2, i-1, i]
        v1 = coeff_table.at("c10")*sol[im2] + coeff_table.at("c11")*sol[im1] + coeff_table.at("c12")*sol[i0];
        // 3) stencil for r=0: [i-1, i, i+1]
        v0 = coeff_table.at("c00")*sol[im1] + coeff_table.at("c01")*sol[i0] + coeff_table.at("c02")*sol[ip1];
        // calculate the weight
        beta0 = get_beta(sol[im1], sol[i0], sol[ip1]);
        beta1 = get_beta(sol[im2], sol[im1], sol[i0]);
        beta2 = get_beta(sol[im3], sol[im2], sol[im1]);
        alpha0 = 0.3 / ((epsilon+beta0)*(epsilon+beta0));
        alpha1 = 0.6 / ((epsilon+beta1)*(epsilon+beta1));
        alpha2 = 0.1 / ((epsilon+beta2)*(epsilon+beta2));
        w0 = alpha0 / (alpha0+alpha1+alpha2);
        w1 = alpha1 / (alpha0+alpha1+alpha2);
        w2 = alpha2 / (alpha0+alpha1+alpha2);
        // weighted average
        varL[i] = w0*v0 + w1*v1 + w2*v2;

        // For the variable on the right side of the i-th edge, we have:
        // 1) stencil for r=2: [i-2, i-1, i]
        v2 = coeff_table.at("c10")*sol[im2] + coeff_table.at("c11")*sol[im1] + coeff_table.at("c12")*sol[i0];
        // 2) stencil for r=1: [i-1, i, i+1] 
        v1 = coeff_table.at("c00")*sol[im1] + coeff_table.at("c01")*sol[i0] + coeff_table.at("c02")*sol[ip1];
        // 3) stencil for r=0: [i, i+1, i+2]
        v0 = coeff_table.at("c-10")*sol[i0] + coeff_table.at("c-11")*sol[ip1] + coeff_table.at("c-12")*sol[ip2];
        // calculate the weight
        beta0 = get_beta(sol[i0], sol[ip1], sol[ip2]);
        beta1 = get_beta(sol[im1], sol[i0], sol[ip1]);
        beta2 = get_beta(sol[im2], sol[im1], sol[i0]);
        alpha0 = 0.3 / ((epsilon+beta0)*(epsilon+beta0));
        alpha1 = 0.6 / ((epsilon+beta1)*(epsilon+beta1));
        alpha2 = 0.1 / ((epsilon+beta2)*(epsilon+beta2));
        w0 = alpha0 / (alpha0+alpha1+alpha2);
        w1 = alpha1 / (alpha0+alpha1+alpha2);
        w2 = alpha2 / (alpha0+alpha1+alpha2);
        // weighted average
        varR[i] = w0*v0 + w1*v1 + w2*v2;
    }
}


// flux calculation with HLLC approximation for Riemann Problem
void fluxCalc_HLLC(Vector<Var>& varL, Vector<Var>& varR, Vector<Var>& FluxX, int nx){
    double p_tmp, qL, qR, aL, aR;
    double rho_roe, u_roe, H_roe, a_roe;
    double SL, SR, S_star;
    double uL, uR, pL, pR;
    Var sol_star;
    const double gamma = 1.4;

    for (int i = 0; i < nx+1; i++){
        uL = varL[i].v2/varL[i].v1;
        pL = get_p(varL[i]);
        aL = sqrt(gamma*pL/varL[i].v1);
        uR = varR[i].v2/varR[i].v1;
        pR = get_p(varR[i]);
        aR = sqrt(gamma*pR/varR[i].v1);
        
        // Roe-Based Wave Speed Estimates
        tie(rho_roe, u_roe, H_roe, a_roe) = RoeAvg(varL[i], varR[i]);
        SL = u_roe - a_roe;
        SR = u_roe + a_roe;

        if (SL==SR){
            S_star = SL;
        }
        else {
            S_star = (pR - pL + varL[i].v2*(SL-uL) - varR[i].v2*(SR-uR)) / (varL[i].v1*(SL-uL) - varR[i].v1*(SR-uR)); 
        }
        
        if (SL>=0){
            FluxX[i] = get_F(varR[i]);
        }
        else if (S_star>=0){
            sol_star.v1 = varL[i].v1*(SL-uL)/(SL-S_star)*1.0;
            sol_star.v2 = varL[i].v1*(SL-uL)/(SL-S_star)*S_star;
            sol_star.v3 = varL[i].v1*(SL-uL)/(SL-S_star)*(varL[i].v3/varL[i].v1+(S_star-uL)*(S_star+pL/(varL[i].v1*(SL-uL))));
            FluxX[i] =  get_F(varL[i]) + SL*(sol_star-varL[i]);
        }
        else if (SR>=0){
            sol_star.v1 = varR[i].v1*(SR-uR)/(SR-S_star)*1.0;
            sol_star.v2 = varR[i].v1*(SR-uR)/(SR-S_star)*S_star;
            sol_star.v3 = varR[i].v1*(SR-uR)/(SR-S_star)*(varR[i].v3/varR[i].v1+(S_star-uR)*(S_star+pR/(varR[i].v1*(SR-uR))));
            FluxX[i] =  get_F(varR[i]) + SR*(sol_star-varR[i]);
        }
        else {
            FluxX[i] = get_F(varR[i]);
        }
    }
}


void step(Vector<Var>& sol, Vector<Var>& FluxX, double dx, double dt, int nx, const map<string, double> coeff_table){
    Vector<Var> sol_tmp(nx);
    Vector<Var> varL(nx+1);  // reconstructed variable on the left side of the edge
    Vector<Var> varR(nx+1);  // reconstructed variable on the right side of the edge

    varRec_WENO(sol, varL, varR, nx, coeff_table);
    fluxCalc_HLLC(varL, varR, FluxX, nx);
    for (int i=0; i<nx; i++){
        sol_tmp[i] = sol[i] - dt/dx * (FluxX[i+1] - FluxX[i]);
    }

    varRec_WENO(sol_tmp, varL, varR, nx, coeff_table);
    fluxCalc_HLLC(varL, varR, FluxX, nx);
    for (int i=0; i<nx; i++){
        sol_tmp[i] = 3.0/4.0*sol[i] + 1.0/4.0*(sol_tmp[i] - dt/dx*(FluxX[i+1]-FluxX[i]));
    }

    varRec_WENO(sol_tmp, varL, varR, nx, coeff_table);
    fluxCalc_HLLC(varL, varR, FluxX, nx);
    for (int i=0; i<nx; i++){
        sol[i] = 1.0/3.0*sol[i] + 2.0/3.0*(sol_tmp[i] - dt/dx*(FluxX[i+1]-FluxX[i]));
    }
}


/* g++ -Og -o HW9_WENO HW9_WENO.cpp HW9_WENO.hpp */
/* ./HW9_WENO 0.2 200 0.001 -0.5 0.5 0.0 0.01 HLLC */
/* ./HW9_WENO 0.2 200 0.001 -0.5 0.5 0.0 0.01 ExactRP */
/* ./HW9_WENO 0.2 400 0.0005 -0.5 0.5 0.0 0.005 HLLC */
/* ./HW9_WENO 0.2 400 0.0005 -0.5 0.5 0.0 0.005 ExactRP */
/* ./HW9_WENO 0.2 1000 0.0002 -0.5 0.5 0.0 0.002 HLLC */
/* ./HW9_WENO 0.2 1000 0.0002 -0.5 0.5 0.0 0.002 ExactRP */
/* ./HW9_WENO 0.2 2000 0.0001 -0.5 0.5 0.0 0.001 HLLC */
/* ./HW9_WENO 0.2 2000 0.0001 -0.5 0.5 0.0 0.001 ExactRP */
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
    const string solverScheme = argv[8];

    cout << "dx: " << dx << "\tdt: " << dt << endl;
    cout << "ReportTime Interval: " << reportStep*dt << endl;
    cout << "Solution Domain: " << "[" << Xstart << ", " << Xend << "]" << endl;

    int nx = int((Xend - Xstart) / dx);

    const map<string, double> weno_table = {
        {"c-10", 11.0/6.0}, {"c-11", -7.0/6.0}, {"c-12", 1.0/3.0},
        {"c00", 1.0/3.0}, {"c01", 5.0/6.0}, {"c02", -1.0/6.0},
        {"c10", -1.0/6.0}, {"c11", 5.0/6.0}, {"c12", 1.0/3.0},
        {"c20", 1.0/3.0}, {"c21", -7.0/6.0}, {"c22", 11.0/6.0}
    };

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
    // 5th WENO scheme with HLLC flux
    else {
        for (double t = 0; t < elapsedTime; t+=dt){
            step(U, FluxX, dx, dt, nx, weno_table);
            stepCounter += 1;
            // Report result
            if (stepCounter % reportStep==0){
                // ------report rho
                string name = "HLLC_Rho_" + to_string(double(stepCounter*dt)) + ".txt";
                cout << "Time: " << t+dt << ", " << name << " recorded" << endl;
                outFile.open(name, ios::out);
                for (int i = 0; i < nx; i++){
                    outFile << U[i].v1 << " ";
                }
                outFile.close();
                // ------report U
                name = "HLLC_U_" + to_string(double(stepCounter*dt)) + ".txt";
                cout << "Time: " << t+dt << ", " << name << " recorded" << endl;
                outFile.open(name, ios::out);
                for (int i = 0; i < nx; i++){
                    outFile << U[i].v2/U[i].v1 << " ";
                }
                outFile.close();
                // ------report P
                name = "HLLC_P_" + to_string(double(stepCounter*dt)) + ".txt";
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