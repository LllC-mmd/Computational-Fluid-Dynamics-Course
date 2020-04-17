#include <math.h>
#include <fstream>
#include <iostream>
#include <string>
#include <array>
#include <tuple>
#include <algorithm>
#include "HW6_RiemannProblem.hpp"


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

tuple<double, double, double> RoeAvg(const Var& A, const Var& B){
    double u_roe, H_roe, a_roe;
    double Ha, Hb;
    const double gamma = 1.4;

    Ha = get_H(A);
    Hb = get_H(B);
    u_roe = (A.v2/sqrt(A.v1)+B.v2/sqrt(B.v1)) / (sqrt(A.v1)+sqrt(B.v1)); 
    H_roe = (sqrt(A.v1)*Ha+sqrt(B.v1)*Hb) / (sqrt(A.v1)+sqrt(B.v1));
    a_roe = sqrt((gamma-1.0)*(H_roe-0.5*u_roe*u_roe));
    return make_tuple(u_roe, H_roe, a_roe);
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

    double p_pv = max(0.0, 0.5*(pL+pR) - 0.125*(uR-uL)*(UL.v1+UR.v1)*(aL+aR));
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


void fluxCalc_Godunov(Vector<Var>& sol, Vector<Var>& FluxX, double dx, double dt, int nx, double tol, double maxiter){
    int im1, i0;
    for (int i = 0; i < nx+1; i++){
        im1 = id_map(i-1, 0, nx-1);
        i0 = id_map(i, 0, nx-1);
        FluxX[i] = get_F(RP_solver(sol[im1], sol[i0], 0.0, dt, tol, maxiter));
    }
}

void fluxCalc_HLL(Vector<Var>& sol, Vector<Var>& FluxX, double dx, double dt, int nx, double tol, double maxiter){
    int im1, i0;
    double u_roe, H_roe, a_roe;
    double SL, SR;
    for (int i = 0; i < nx+1; i++){
        im1 = id_map(i-1, 0, nx-1);
        i0 = id_map(i, 0, nx-1);
        tie(u_roe, H_roe, a_roe) = RoeAvg(sol[im1], sol[i0]);
        SL = u_roe - a_roe;
        SR = u_roe + a_roe;
        if (SL>=0){
            FluxX[i] = get_F(sol[im1]);
        }
        else if (SR>=0){
            FluxX[i] = 1.0/(SR-SL) * (SR*get_F(sol[im1]) - SL*get_F(sol[i0]) + SL*SR*(sol[i0]-sol[im1]));
        }
        else {
            FluxX[i] = get_F(sol[i0]);
        }
    }
}

void fluxCalc_HLLC(Vector<Var>& sol, Vector<Var>& FluxX, double dx, double dt, int nx, double tol, double maxiter){
    int im1, i0;
    double u_roe, H_roe, a_roe;
    double SL, SR, S_star;
    double uL, uR, pL, pR;
    Var sol_star;

    for (int i = 0; i < nx+1; i++){
        im1 = id_map(i-1, 0, nx-1);
        i0 = id_map(i, 0, nx-1);
        tie(u_roe, H_roe, a_roe) = RoeAvg(sol[im1], sol[i0]);
        SL = u_roe - a_roe;
        SR = u_roe + a_roe;
        uL = sol[im1].v2/sol[im1].v1;
        pL = get_p(sol[im1]);
        uR = sol[i0].v2/sol[i0].v1;
        pR = get_p(sol[i0]);
        S_star = (pR - pL + sol[im1].v2*(SL-uL) - sol[i0].v2*(SR-uR)) / (sol[im1].v1*(SL-uL) - sol[i0].v1*(SR-uR)); 
        if (SL>=0){
            FluxX[i] = get_F(sol[im1]);
        }
        else if (S_star>=0){
            sol_star.v1 = sol[im1].v1*(SL-uL)/(SL-S_star)*1.0;
            sol_star.v2 = sol[im1].v1*(SL-uL)/(SL-S_star)*S_star;
            sol_star.v3 = sol[im1].v1*(SL-uL)/(SL-S_star)*(sol[im1].v3/sol[im1].v1+(S_star-uL)*(S_star+pL/(sol[im1].v1*(SL-uL))));
            FluxX[i] =  get_F(sol[im1]) + SL*(sol_star-sol[im1]);
        }
        else if (SR>=0){
            sol_star.v1 = sol[i0].v1*(SR-uR)/(SR-S_star)*1.0;
            sol_star.v2 = sol[i0].v1*(SR-uR)/(SR-S_star)*S_star;
            sol_star.v3 = sol[i0].v1*(SR-uR)/(SR-S_star)*(sol[i0].v3/sol[i0].v1+(S_star-uR)*(S_star+pR/(sol[i0].v1*(SR-uR))));
            FluxX[i] =  get_F(sol[i0]) + SR*(sol_star-sol[i0]);
        }
        else {
            FluxX[i] = get_F(sol[i0]);
        }
    }
}


void step(Vector<Var>& sol, Vector<Var>& FluxX, double dx, double dt, int nx, double tol, double maxiter, const string scheme){
    if (scheme=="Godunov"){
        fluxCalc_Godunov(sol, FluxX, dx, dt, nx, tol, maxiter);
    }
    else if (scheme=="HLL"){
        fluxCalc_HLL(sol, FluxX, dx, dt, nx, tol, maxiter);
    }
    else if (scheme=="HLLC"){
        fluxCalc_HLLC(sol, FluxX, dx, dt, nx, tol, maxiter);
    }
    else {
        cout << "Unknown scheme" << endl;
    }
    
    for (int i = 0; i < nx; i++){
        sol[i] = sol[i] - dt/dx*(FluxX[i+1] - FluxX[i]);
    }
}

/* g++ -Og -o HW6_RiemannProblem HW6_RiemannProblem.cpp HW6_RiemannProblem.hpp */
/* ./HW6_RiemannProblem 0.3 100 0.001 0.0 1.0 0.3 0.01 HW6Case1 Godunov */
/* ./HW6_RiemannProblem 0.15 50 0.001 0.0 1.0 0.5 0.01 HW6Case2 Godunov */
/* ./HW6_RiemannProblem 0.012 40 0.001 0.0 1.0 0.5 0.01 HW6Case3 Godunov */
/* ./HW6_RiemannProblem 0.035 350 0.0001 0.0 1.0 0.4 0.01 HW6Case4 Godunov */
/* ./HW6_RiemannProblem 0.012 120 0.0001 0.0 1.0 0.8 0.01 HW6Case5 Godunov */
/* ./HW6_RiemannProblem 0.012 40 0.0001 0.0 1.0 0.5 0.01 HW7Case3 HLLC */
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
    const string CaseNo = argv[8];
    const string scheme = argv[9];

    cout << "dx: " << dx << "\tdt: " << dt << endl;
    cout << "ReportTime Interval: " << reportStep*dt << endl;
    cout << "Solution Domain: " << "[" << Xstart << ", " << Xend << "]" << endl;

    int nx = int((Xend - Xstart) / dx);

    Vector<Var> U(nx);
    // place Ui in the middle of F_i and F_i+1
    Vector<Var> FluxX(nx+1);

    // These cases are from Toro's book: 6.4 Numerical Results and Discussion for HW6
    // and 10.8 Numerical Results for HW7
    if (CaseNo=="HW6Case1"){
        rhoL = 1.0;
        uL = 0.75;
        pL = 1.0;
        rhoR = 0.125;
        uR = 0.0;
        pR = 0.1;
    }
    else if (CaseNo=="HW6Case2"){
        rhoL = 1.0;
        uL = -2.0;
        pL = 0.4;
        rhoR = 1.0;
        uR = 2.0;
        pR = 0.4;
    }
    else if (CaseNo=="HW6Case3"){
        rhoL = 1.0;
        uL = 0.0;
        pL = 1000.0;
        rhoR = 1.0;
        uR = 0.0;
        pR = 0.01;
    }
    else if (CaseNo=="HW6Case4"){
        rhoL = 5.99924;
        uL = 19.5975;
        pL = 460.894;
        rhoR = 5.99242;
        uR = -6.19633;
        pR = 46.0950;
    }
    else if (CaseNo=="HW6Case5"){
        rhoL = 1.0;
        uL = -19.59745;
        pL = 1000.0;
        rhoR = 1.0;
        uR = -19.59745;
        pR = 0.01;
    }
    else if (CaseNo=="HW7Case1"){
        rhoL = 1.0;
        uL = 0.0;
        pL = 1.0;
        rhoR = 0.125;
        uR = 0.0;
        pR = 0.1;
    }
    else if (CaseNo=="HW7Case2"){
        rhoL = 1.0;
        uL = -2.0;
        pL = 0.4;
        rhoR = 1.0;
        uR = 2.0;
        pR = 0.4;
    }
    else if (CaseNo=="HW7Case3"){
        rhoL = 1.0;
        uL = 0.0;
        pL = 1000.0;
        rhoR = 1.0;
        uR = 0.0;
        pR = 0.01;
    }
    else if (CaseNo=="HW7Case4"){
        rhoL = 1.0;
        uL = 0.0;
        pL = 0.01;
        rhoR = 1.0;
        uR = 0.0;
        pR = 100.0;
    }
    else if (CaseNo=="HW7Case5"){
        rhoL = 5.99924;
        uL = 19.5975;
        pL = 460.894;
        rhoR = 5.99242;
        uR = -6.19633;
        pR = 46.0950;
    }
    else {
        cout << "Unknown Case" << endl;
    }

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
    // Case1: Exact Solution of Riemann Problem
    if (scheme=="ExactRP"){
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
                string name = CaseNo + "Exact_Rho_" + to_string(double(stepCounter*dt)) + ".txt";
                cout << "Time: " << t+dt << ", " << name << " recorded" << endl;
                outFile.open(name, ios::out);
                for (int i = 0; i < nx; i++){
                    outFile << U[i].v1 << " ";
                }
                outFile.close();
                // ------report U
                name = CaseNo + "Exact_U_" + to_string(double(stepCounter*dt)) + ".txt";
                cout << "Time: " << t+dt << ", " << name << " recorded" << endl;
                outFile.open(name, ios::out);
                for (int i = 0; i < nx; i++){
                    outFile << U[i].v2/U[i].v1 << " ";
                }
                outFile.close();
                // ------report P
                name = CaseNo + "Exact_P_" + to_string(double(stepCounter*dt)) + ".txt";
                cout << "Time: " << t+dt << ", " << name << " recorded" << endl;
                outFile.open(name, ios::out);
                for (int i = 0; i < nx; i++){
                    outFile << get_p(U[i]) << " ";
                }
                outFile.close();
            }
        }
    }
    // Case2: FVM: Godunov, HLL, HLLC
    else {
        for (double t = 0; t < elapsedTime; t+=dt){
            step(U, FluxX, dx, dt, nx, tol, maxiter, scheme);
            stepCounter += 1;
            // Report result
            if (stepCounter % reportStep==0){
                // ------report rho
                string name = CaseNo + scheme +"_Rho_" + to_string(double(stepCounter*dt)) + ".txt";
                cout << "Time: " << t+dt << ", " << name << " recorded" << endl;
                outFile.open(name, ios::out);
                for (int i = 0; i < nx; i++){
                    outFile << U[i].v1 << " ";
                }
                outFile.close();
                // ------report U
                name = CaseNo + scheme + "_U_" + to_string(double(stepCounter*dt)) + ".txt";
                cout << "Time: " << t+dt << ", " << name << " recorded" << endl;
                outFile.open(name, ios::out);
                for (int i = 0; i < nx; i++){
                    outFile << U[i].v2/U[i].v1 << " ";
                }
                outFile.close();
                // ------report P
                name = CaseNo + scheme + "_P_" + to_string(double(stepCounter*dt)) + ".txt";
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

