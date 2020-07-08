#include "fvm.hpp"
#include "global.hpp"
#include <iostream>
#include <math.h>
#include <string>
#include <tuple>


using namespace std;


// rotate a vector according to the rotation matrix
U rotate(U u, double theta){
    U u_r;

    u_r.h = u.h;
    u_r.hu = cos(theta)*u.hu + sin(theta)*u.hv;
    u_r.hv = -sin(theta)*u.hu + cos(theta)*u.hv;
    u_r.eta = u.eta;

    return u_r;
}


// rotate a vector back according to the inverse of the rotation matrix
U rotate_back(U u, double theta){
    U u_back;

    u_back.h = u.h;
    u_back.hu = cos(theta)*u.hu - sin(theta)*u.hv;
    u_back.hv = sin(theta)*u.hu + cos(theta)*u.hv;
    u_back.eta = u.eta;

    return u_back;
}

FluxF rotate_back(FluxF f, double theta){
    FluxF f_back;

    f_back.F1 = f.F1;
    f_back.F2 = cos(theta)*f.F2 - sin(theta)*f.F3;
    f_back.F3 = sin(theta)*f.F2 + cos(theta)*f.F3;

    return f_back;
}


tuple<U, U> grad_Calc(const U& u1, double x1, double y1, const U& u2, double x2, double y2, const U& u3, double x3, double y3){
    double a = 1.0 / ((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1));
    double J11 = a * (y3-y1);
    double J12 = a * (y1-y2);
    double J21 = a * (x1-x3);
    double J22 = a * (x2-x1);

    U gradX = J11 * (u2 - u1) + J12 * (u3 - u1);
    U gradY = J21 * (u2 - u1) + J22 * (u3 - u1);
    
    return make_tuple(gradX, gradY);
}


double Edge::getEdgeSpectralRadius(const U& UL, const U& UR){
    double h_roe, u_roe;
    U UL_rot, UR_rot;
    double lambda;

    UL_rot = rotate(UL, theta_n);
    UR_rot = rotate(UR, theta_n);
    tie(h_roe, u_roe) = RoeAvg(UL_rot, UR_rot);
    lambda = fabs(u_roe)+sqrt(g*h_roe);

    return lambda;
}


double Cell::getCellSpectralRadius(){
    int eId, ncId;
    double lambda = 0.0;

    for (int i = 0; i < edgeId_Cell.size(); i++){
        eId = edgeId_Cell[i];
        ncId = neighborCellId_Cell[i];
        if (ncId >= 1){
            lambda = max(globalEdge[eId-1].getEdgeSpectralRadius(_Uc, globalCell[ncId-1].getValUc()), lambda);
        }
    }
    return lambda;
}

// ********************************* Reconstruction *********************************
// [1] MUSCL reconstruction
void Cell::CellRec_MUSCL(int alpha){
    U gradX, gradY, U_nodal;
    int nId, eId, ncId;
    double xn, yn, zn;
    double eta_max, eta_min, hu_max, hu_min, hv_max, hv_min;
    double eta_limiter = 1.0;
    double hu_limiter = 1.0;
    double hv_limiter = 1.0;

    vector<U> u(4);
    vector<U> u_tmp(3);
    vector<double> x(3);
    vector<double> y(3);
    
    // ---unlimited gradient calculation
    int c = 0;
    for (int i = 0; i < neighborCellId_Cell.size(); i++){
        ncId = neighborCellId_Cell[i];
        eId = edgeId_Cell[i];
        if ((neighborCellId_Cell[i] >= 1) && (globalEdge[eId-1].boundaryFlag==0)){
            u[c] = globalCell[ncId-1].getValU_alpha(alpha);
            x[c] = globalCell[ncId-1].cpx;
            y[c] = globalCell[ncId-1].cpy;
            c += 1;
        }
    }

    // ---switch to the 1st scheme if the cell only has 1 neighbour cells
    if (c==1){
        gradX.hu = 0.0;
        gradX.hv = 0.0;
        gradX.eta = 0.0;
        gradY.hu = 0.0;
        gradY.hv = 0.0;
        gradY.eta = 0.0;
    }
    else {
        // ------the cell has 2 neighbour cells
        if (c==2){
            u[c] = _U_alpha[alpha];
            x[c] = cpx;
            y[c] = cpy;
            u[c+1] = _U_alpha[alpha];
        }
        else{
            u[c] = _U_alpha[alpha];
        }

        if ((x[1]-x[0])*(y[2]-y[0]) - (x[2]-x[0])*(y[1]-y[0]) > 0){
            tie(gradX, gradY) = grad_Calc(u[0], x[0], y[0], u[1], x[1], y[1], u[2], x[2], y[2]);
        }
        else{
            tie(gradX, gradY) = grad_Calc(u[0], x[0], y[0], u[2], x[2], y[2], u[1], x[1], y[1]);
        }

        // ------extrapolation using the unlimited gradients
        for (int i = 0; i < edgeId_Cell.size(); i++){
            eId = edgeId_Cell[i];
            u_tmp[i] = _U_alpha[alpha] + (globalEdge[eId-1].cx-cpx)*gradX + (globalEdge[eId-1].cy-cpy)*gradY;
        }
        eta_max = (*max_element(u.begin(), u.end(), [](U a, U b){
                        return (a.eta < b.eta);})).eta;
        eta_min = (*max_element(u.begin(), u.end(), [](U a, U b){
                        return (a.eta > b.eta);})).eta;
        hu_max = (*max_element(u.begin(), u.end(), [](U a, U b){
                        return (a.hu < b.hu);})).hu;
        hu_min = (*max_element(u.begin(), u.end(), [](U a, U b){
                        return (a.hu > b.hu);})).hu;
        hv_max = (*max_element(u.begin(), u.end(), [](U a, U b){
                        return (a.hv < b.hv);})).hv;
        hv_min = (*max_element(u.begin(), u.end(), [](U a, U b){
                        return (a.hv > b.hv);})).hv;
        // ------gradient limiting
        for (int i = 0; i < edgeId_Cell.size(); i++){
            // ---------eta limiting
            if (u_tmp[i].eta == _U_alpha[alpha].eta){
                eta_limiter = min(1.0, eta_limiter);
            }
            else if (u_tmp[i].eta < _U_alpha[alpha].eta){
                eta_limiter = min(min(1.0, min(eta_min-_U_alpha[alpha].eta, 0.0)/(u_tmp[i].eta-_U_alpha[alpha].eta)), eta_limiter);
            }
            else{
                eta_limiter = min(min(1.0, max(eta_max-_U_alpha[alpha].eta, 0.0)/(u_tmp[i].eta-_U_alpha[alpha].eta)), eta_limiter);
            }
            // ---------hu limiting
            if (u_tmp[i].hu == _U_alpha[alpha].hu){
                hu_limiter = min(1.0, hu_limiter);
            }
            else if (u_tmp[i].hu < _U_alpha[alpha].hu){
                hu_limiter = min(min(1.0, min(hu_min-_U_alpha[alpha].hu, 0.0)/(u_tmp[i].hu-_U_alpha[alpha].hu)), hu_limiter);
            }
            else{
                hu_limiter = min(min(1.0, max(hu_max-_U_alpha[alpha].hu, 0.0)/(u_tmp[i].hu-_U_alpha[alpha].hu)), hu_limiter);
            }
            // ---------hv limiting
            if (u_tmp[i].hv == _U_alpha[alpha].hv){
                hv_limiter = min(1.0, hv_limiter);
            }
            else if (u_tmp[i].hv < _U_alpha[alpha].hv){
                hv_limiter = min(min(1.0, min(hv_min-_U_alpha[alpha].hv, 0.0)/(u_tmp[i].hv-_U_alpha[alpha].hv)), hv_limiter);
            }
            else{
                hv_limiter = min(min(1.0, max(hv_max-_U_alpha[alpha].hv, 0.0)/(u_tmp[i].hv-_U_alpha[alpha].hv)), hv_limiter);
            }
        }
        gradX.eta *= eta_limiter;
        gradX.hu *= hu_limiter;
        gradX.hv *= hv_limiter;
        gradY.eta *= eta_limiter;
        gradY.hu *= hu_limiter;
        gradY.hv *= hv_limiter;
    }

    for (int i = 0; i < nodeId_Cell.size(); i++){
        nId = nodeId_Cell[i];
        xn = globalNode[nId-1].px;
        yn = globalNode[nId-1].py;
        zn = globalNode[nId-1].zb;
        U_nodal = _U_alpha[alpha] + (xn-cpx)*gradX + (yn-cpy)*gradY;
        if (U_nodal.eta - zn < 0){
            gradX = 0.0 * gradX;
            gradY = 0.0 * gradY;
            break;
        };
    }

    for (int i = 0; i < edgeId_Cell.size(); i++){
        eId = edgeId_Cell[i];
        if (globalEdge[eId-1].getOCell() == cId){
            globalEdge[eId-1].UL_GP[0] = _U_alpha[alpha] + (globalEdge[eId-1].cx-cpx)*gradX + (globalEdge[eId-1].cy-cpy)*gradY;
            globalEdge[eId-1].UL_GP[0].h = globalEdge[eId-1].UL_GP[0].eta - globalEdge[eId-1].zm;
        }
        else{
            globalEdge[eId-1].UR_GP[0] = _U_alpha[alpha] + (globalEdge[eId-1].cx-cpx)*gradX + (globalEdge[eId-1].cy-cpy)*gradY;
            globalEdge[eId-1].UR_GP[0].h = globalEdge[eId-1].UR_GP[0].eta - globalEdge[eId-1].zm;
        }
    }
}


void Edge::reconstruct_BC(int alpha, double tol){
    if (boundaryFlag==6){
        for (int i = 0; i < UL_GP.size()/2; i++){
            reconstruct_double_solid(i, alpha);
        }
    }
    else{
        for (int i = 0; i < UL_GP.size(); i++){
            if (boundaryFlag==1){
                reconstruct_free(i);
            }
            else if (boundaryFlag==2){
                reconstruct_solid(i);
            }
            // for super-critical flow, u > sqrt(gh)
            else if (fabs(cos(theta_n)*UL_GP[i].hu+sin(theta_n)*UL_GP[i].hv) > UL_GP[i].h*sqrt(g*UL_GP[i].h)){
                reconstruct_OpenBC_super(i);
            }
            // for sub-critical flow
            else if (boundaryFlag==3){
                reconstruct_OpenBC_h(i);
            }
            else if (boundaryFlag==4){
                reconstruct_OpenBC_u(i);
            }
            else if (boundaryFlag==5){
                reconstruct_OpenBC_hu(i, tol);
            }
            else if (boundaryFlag==6){
                reconstruct_double_solid(i, alpha);
            }
        }
    }
}; 


// Boundary treatment: supercritical flow 
void Edge::reconstruct_OpenBC_super(int i){
    UR_GP[i] = UL_GP[i];
}

// Boundary treatment: free outflow
void Edge::reconstruct_free(int i){
    UR_GP[i] = UL_GP[i];
}

// Boundary treatment: no-slip solid wall
void Edge::reconstruct_solid(int i){
    UR_GP[i].hu = 0.0;
    UR_GP[i].hv = 0.0;
    UR_GP[i].eta = UL_GP[i].eta;
    UR_GP[i].h = max(0.0, UL_GP[i].eta - zm);
}

// Boundary treatment: h is given
void Edge::reconstruct_OpenBC_h(int i){
    U UL_rot, UR_rot;
    double hL, uL_rot, vL_rot;

    _hBC = 5.0;

    UL_rot = rotate(UL_GP[i], theta_n);
    hL = UL_rot.h;
    uL_rot = (hL==0.0) ? 0.0 : UL_rot.hu/UL_rot.h;
    vL_rot = (hL==0.0) ? 0.0 : UL_rot.hv/UL_rot.h;

    UR_rot.h = _hBC;
    UR_rot.hu = _hBC*uL_rot + 2*_hBC*cos(theta_n)*(sqrt(g*hL)-sqrt(g*_hBC));
    UR_rot.hv = _hBC*vL_rot;

    UR_GP[i] = rotate_back(UR_rot, theta_n);
}

// Boundary treatment: u is given
void Edge::reconstruct_OpenBC_u(int i){
    U UR_rot, UL_rot;
    double hL_rot, uL_rot, vL_rot;

    UL_rot = rotate(UL_GP[i], theta_n);

    hL_rot = UL_rot.h;
    uL_rot = (UL_rot.h==0.0) ? 0.0 : UL_rot.hu/UL_rot.h;
    vL_rot = (UL_rot.h==0.0) ? 0.0 : UL_rot.hv/UL_rot.h;
    
    UR_rot.h = pow(uL_rot+2*sqrt(g*hL_rot)-_uBC, 2)/(4*g);
    UR_rot.hu = UR_rot.h*_uBC;
    UR_rot.hv = UR_rot.h*vL_rot;

    UR_GP[i] = rotate_back(UR_rot, theta_n);
}

// Boundary treatment: hu is given
void Edge::reconstruct_OpenBC_hu(int i, double tol){
    U UR_rot, UL_rot;
    double h0;
    double C, h_sqrt, res;
    double hL_rot, uL_rot, vL_rot;

    for (int i = 0; i < UL_GP.size(); i++){
        UL_rot = rotate(UL_GP[i], theta_n);

        hL_rot = UL_rot.h;
        uL_rot = (UL_rot.h==0) ? 0.0 : UL_rot.hu/UL_rot.h;
        vL_rot = (UL_rot.h==0) ? 0.0 : UL_rot.hv/UL_rot.h;

        C = uL_rot + 2*sqrt(g*hL_rot);
        if (C >= 0){
            h0 = pow(C+1e-3, 2)/g;
        }
        else{
            h0 = pow(1e-3, 2)/g;
        }

        // solve hR_rot, uR_rot using Newton-Raphson method:
        // 1) hR_rot * uR_rot = huBC
        // 2) uR_rot + 2*sqrt(g*hR_rot) = uL_rot + 2*sqrt(g*hL_rot)
        h_sqrt = sqrt(h0);
        res = 2*sqrt(g)*pow(h_sqrt, 3) - C*pow(h_sqrt, 2) + _huBC;
        while (res > tol){
            h_sqrt -= res/(6*sqrt(g)-2*C*h_sqrt);
            res = 2*sqrt(g)*pow(h_sqrt, 3) - C*pow(h_sqrt, 2) + _huBC;
        }

        UR_rot.h = pow(h_sqrt, 2);
        UR_rot.hu = _huBC;
        UR_rot.hv = UR_rot.h*vL_rot;

        UR_GP[i] = rotate_back(UR_rot, theta_n);
    }
}

void Edge::reconstruct_double_solid(int i, int alpha){
    int nGP = UR_GP.size() / 2;
    U u_owner = globalCell[cellId_Edge[0]-1].getValU_alpha(alpha);
    U u_neighbour = globalCell[cellId_Edge[1]-1].getValU_alpha(alpha);

    UR_GP[i] = u_owner;
    UR_GP[i].hu = 0.0;
    UR_GP[i].hv = 0.0;
    UR_GP[i].h = max(0.0, UR_GP[i].eta - zm);
    UR_GP[i+nGP] = u_neighbour;
    UR_GP[i+nGP].hu = 0.0;
    UR_GP[i+nGP].hv = 0.0;
    UR_GP[i+nGP].h = max(0.0, UR_GP[i+nGP].eta - zm);
}


void Edge::reconstruct_adjust(const U& u_owner, const U& u_neighbour){
    double uL, vL, uR, vR;
    double hL_temp, hR_temp;

    for (int i = 0; i < UL_GP.size(); i++){
        // ---Switch to the 1st scheme: 
        // ------[1] negative reconstructed water depth
        if (UL_GP[i].h < 0){
            UL_GP[i].eta = u_owner.eta;
            UL_GP[i].h = UL_GP[i].eta - zm;
            UL_GP[i].hu = u_owner.hu;
            UL_GP[i].hv = u_owner.hv;
            // ---Non-negative depth reconstruction
            if (UL_GP[i].h < 0.0){
                UL_GP[i].eta = zm;
                UL_GP[i].h = 0.0;
                UL_GP[i].hu = 0.0;
                UL_GP[i].hv = 0.0;
            }
        }

        if (UR_GP[i].h < 0){
            UR_GP[i].eta = u_neighbour.eta;
            UR_GP[i].h = UR_GP[i].eta - zm;
            UR_GP[i].hu = u_neighbour.hu;
            UR_GP[i].hv = u_neighbour.hv;
            if (UR_GP[i].h < 0.0){
                UR_GP[i].eta = zm;
                UR_GP[i].h = 0.0;
                UR_GP[i].hu = 0.0;
                UR_GP[i].hv = 0.0;
            }
        }

        uL = (UL_GP[i].h==0) ? 0.0 : UL_GP[i].hu/UL_GP[i].h;
        vL = (UL_GP[i].h==0) ? 0.0 : UL_GP[i].hv/UL_GP[i].h;
        uR = (UR_GP[i].h==0) ? 0.0 : UR_GP[i].hu/UR_GP[i].h;
        vR = (UR_GP[i].h==0) ? 0.0 : UR_GP[i].hv/UR_GP[i].h;

        // ---Switch to the 1st scheme if the velocity is too high
        if ((uL > u_max) || (vL > u_max)){
            UL_GP[i].eta = u_owner.eta;
            UL_GP[i].h = UL_GP[i].eta - zm;
            UL_GP[i].hu = u_owner.hu;
            UL_GP[i].hv = u_owner.hv;
            if (UL_GP[i].h <= 0.0){
                UL_GP[i].eta = zm;
                UL_GP[i].h = 0.0;
                UL_GP[i].hu = 0.0;
                UL_GP[i].hv = 0.0;
            }
        }
        if ((uR > u_max) || (vR > u_max)){
            UR_GP[i].eta = u_neighbour.eta;
            UR_GP[i].h = UR_GP[i].eta - zm;
            UR_GP[i].hu = u_neighbour.hu;
            UR_GP[i].hv = u_neighbour.hv;
            if (UR_GP[i].h <= 0.0){
                UR_GP[i].eta = zm;
                UR_GP[i].h = 0.0;
                UR_GP[i].hu = 0.0;
                UR_GP[i].hv = 0.0;
            }
        }
    }
}


// [2] Variational Reconstruction
void Cell::VarRecUpdate(){
    int ncId;

    arma::Col<double> b_eta_tmp;
    arma::Col<double> b_hu_tmp;
    arma::Col<double> b_hv_tmp;
    b_eta_tmp.resize(coef_eta.size());
    b_eta_tmp.fill(0.0);
    b_hu_tmp.resize(coef_hu.size());
    b_hu_tmp.fill(0.0);
    b_hv_tmp.resize(coef_hv.size());
    b_hv_tmp.fill(0.0);
 
    int Nc = (vfv_order+1)*(vfv_order+2)/2 - 1;
    arma::Mat<double> I_mat = arma::eye(Nc,Nc); 

    arma::Mat<double> A_eta_inv = solve(A_eta, I_mat);
    arma::Mat<double> A_hu_inv = solve(A_hu, I_mat);
    arma::Mat<double> A_hv_inv = solve(A_hv, I_mat);

    // block SOR
    coef_eta *= (1 - omega);
    coef_hu *= (1 - omega);
    coef_hv *= (1 - omega);
    // ---sum over edges
    int num_nc = 0;
    for (int j = 0; j < 3; j++){
        if (edgeBCTag_Cell[j] == 0){
            ncId = neighborCellId_Cell[j];
            // ------ update the coef of h
            coef_eta += omega * (A_eta_inv * B_eta[num_nc] * globalCell[ncId-1].coef_eta); 
            b_eta_tmp += b_eta_unit.col(j) * (globalCell[ncId-1].getValU_inner().eta - _u_inner.eta);
            // ------ update the coef of hu
            coef_hu += omega * (A_hu_inv * B_hu[num_nc] * globalCell[ncId-1].coef_hu); 
            b_hu_tmp += b_hu_unit.col(j) * (globalCell[ncId-1].getValU_inner().hu - _u_inner.hu);
            // ------ update the coef of hv
            coef_hv += omega * (A_hv_inv * B_hv[num_nc] * globalCell[ncId-1].coef_hv); 
            b_hv_tmp += b_hv_unit.col(j) * (globalCell[ncId-1].getValU_inner().hv - _u_inner.hv);
            num_nc += 1;
        }
        else if (edgeBCTag_Cell[j] == 2){
            b_hu_tmp += b_hu_unit.col(j) * _u_inner.hu;
            b_hv_tmp += b_hv_unit.col(j) * _u_inner.hv;
        }
    }
    coef_eta += omega * (A_eta_inv * b_eta_tmp); 
    coef_hu += omega * (A_hu_inv * b_hu_tmp); 
    coef_hv += omega * (A_hv_inv * b_hv_tmp); 
}


void Cell::CellRec_VarRec(const U& u_c){
    int eId, xf, yf, xt, yt;
    double alpha;

    for (int i = 0; i < edgeId_Cell.size(); i++){
        eId = edgeId_Cell[i];
        xf = globalNode[globalEdge[eId-1].getFrom()-1].px;
        yf = globalNode[globalEdge[eId-1].getFrom()-1].py;
        xt = globalNode[globalEdge[eId-1].getTo()-1].px;
        yt = globalNode[globalEdge[eId-1].getTo()-1].py;
        // set UL_GP
        if (globalEdge[eId-1].getOCell() == cId){
            for (int j = 0; j < globalEdge[eId-1].UL_GP.size(); j++){
                alpha = LineGQ_table.at(vfv_order)[j][1];
                globalEdge[eId-1].UL_GP[j].hu = rec_hu(u_c, xf, yf, xt, yt, alpha);
                globalEdge[eId-1].UL_GP[j].hv = rec_hv(u_c, xf, yf, xt, yt, alpha);
                globalEdge[eId-1].UL_GP[j].eta = rec_eta(u_c, xf, yf, xt, yt, alpha);
                globalEdge[eId-1].UL_GP[j].h = globalEdge[eId-1].UL_GP[j].eta - globalEdge[eId-1].zm;
            }
        }
        // set UR_GP
        else{
            for (int j = 0; j < globalEdge[eId-1].UR_GP.size(); j++){
                alpha = LineGQ_table.at(vfv_order)[j][1];
                globalEdge[eId-1].UR_GP[j].hu = rec_hu(u_c, xf, yf, xt, yt, alpha);
                globalEdge[eId-1].UR_GP[j].hv = rec_hv(u_c, xf, yf, xt, yt, alpha);
                globalEdge[eId-1].UR_GP[j].eta = rec_eta(u_c, xf, yf, xt, yt, alpha);
                globalEdge[eId-1].UR_GP[j].h = globalEdge[eId-1].UR_GP[j].eta - globalEdge[eId-1].zm;
            }
        }   
    }
}


double Cell::rec_hu(const U& u_c, double x1, double y1, double x2, double y2, double alpha){
    int m, n;

    double xm = 0.5*(x1+x2);
    double ym = 0.5*(y1+y2);
    double x_GP = xm + 0.5*alpha*(x2-x1);
    double y_GP = ym + 0.5*alpha*(y2-y1);
    
    double hu = u_c.hu;
    for (int i = 0; i < coef_hu.size(); i++){
        m = taylor_table.at(vfv_order)[i][0];
        n = taylor_table.at(vfv_order)[i][1];
        hu += coef_hu(i) * (pow((x_GP-cpx)/dx_max, m)*pow((y_GP-cpy)/dy_max, n) - basis_mean(i));
    }
    return hu;
}


double Cell::rec_hv(const U& u_c, double x1, double y1, double x2, double y2, double alpha){
    int m, n;

    double xm = 0.5*(x1+x2);
    double ym = 0.5*(y1+y2);
    double x_GP = xm + 0.5*alpha*(x2-x1);
    double y_GP = ym + 0.5*alpha*(y2-y1);
    
    double hv = u_c.hv;
    for (int i = 0; i < coef_hv.size(); i++){
        m = taylor_table.at(vfv_order)[i][0];
        n = taylor_table.at(vfv_order)[i][1];
        hv += coef_hv(i) * (pow((x_GP-cpx)/dx_max, m)*pow((y_GP-cpy)/dy_max, n) - basis_mean(i));
    }
    return hv;
}


double Cell::rec_eta(const U& u_c, double x1, double y1, double x2, double y2, double alpha){
    int m, n;

    double xm = 0.5*(x1+x2);
    double ym = 0.5*(y1+y2);
    double x_GP = xm + 0.5*alpha*(x2-x1);
    double y_GP = ym + 0.5*alpha*(y2-y1);
    
    double eta = u_c.eta;
    // use the reconstruction coefficients of h for eta's reconstruction
    for (int i = 0; i < coef_eta.size(); i++){
        m = taylor_table.at(vfv_order)[i][0];
        n = taylor_table.at(vfv_order)[i][1];
        eta += coef_eta(i) * (pow((x_GP-cpx)/dx_max, m)*pow((y_GP-cpy)/dy_max, n) - basis_mean(i));
    }
    return eta;
}


// 1st reconstruction for dry cell
void Cell::CellRec_dry(const U& u_c){
    int eId;
    for (int i = 0; i < edgeId_Cell.size(); i++){
        eId = edgeId_Cell[i];
        if (globalEdge[eId-1].getOCell() == cId){
            for (int j = 0; j < globalEdge[eId-1].UL_GP.size(); j++){
                globalEdge[eId-1].UL_GP[j] = u_c;
                globalEdge[eId-1].UL_GP[j].h = globalEdge[eId-1].UL_GP[j].eta - globalEdge[eId-1].zm;
            }
        }
        else{
            for (int j = 0; j < globalEdge[eId-1].UR_GP.size(); j++){
                globalEdge[eId-1].UR_GP[j] = u_c;
                globalEdge[eId-1].UR_GP[j].h = globalEdge[eId-1].UR_GP[j].eta - globalEdge[eId-1].zm;
            }
        }
    }
}


void Cell::CellRec_partial(const U& u_c){
    int eId;
    double z1, z2;
    double h_rec, u_rec, v_rec;
    for (int i = 0; i < edgeId_Cell.size(); i++){
        eId = edgeId_Cell[i];
        z1 = globalEdge[eId-1].z_Edge[0];
        z2 = globalEdge[eId-1].z_Edge[1];

        if (u_c.eta <= z1){
            h_rec = 0.0;
        }
        else if (u_c.eta <= z2){
            h_rec = pow(u_c.eta-z1, 2)/2.0/(z2-z1);
        }
        else{
            h_rec = u_c.eta - (z1+z2)/2.0;
        }
        u_rec = (u_c.h==0.0) ? 0.0 : u_c.hu/u_c.h;
        v_rec = (u_c.h==0.0) ? 0.0 : u_c.hv/u_c.h;

        if (globalEdge[eId-1].getOCell() == cId){
            for (int j = 0; j < globalEdge[eId-1].UL_GP.size(); j++){
                globalEdge[eId-1].UL_GP[j].h = h_rec;
                globalEdge[eId-1].UL_GP[j].hu = h_rec*u_rec;
                globalEdge[eId-1].UL_GP[j].hv = h_rec*v_rec;
                globalEdge[eId-1].UL_GP[j].eta = h_rec + globalEdge[eId-1].zm;
            }
        }
        else{
            for (int j = 0; j < globalEdge[eId-1].UL_GP.size(); j++){
                globalEdge[eId-1].UR_GP[j].h = h_rec;
                globalEdge[eId-1].UR_GP[j].hu = h_rec*u_rec;
                globalEdge[eId-1].UR_GP[j].hv = h_rec*v_rec;
                globalEdge[eId-1].UR_GP[j].eta = h_rec + globalEdge[eId-1].zm;
            }
        }
    }
};


// ********************************* Flux calculation *********************************
tuple<double, double> RoeAvg(const U& A, const U& B){
    double h_roe, u_roe;

    h_roe = 0.5*(A.h+B.h);
    if (pow(A.h, 2) + pow(B.h, 2) > 0.0){
        if (A.h*B.h > 0.0){
            u_roe = (A.hu/sqrt(A.h)+B.hu/sqrt(B.h)) / (sqrt(A.h)+sqrt(B.h));
        }
        else if (A.h > 0.0){
            u_roe = A.hu/A.h;
        }
        else{
            u_roe = B.hu/B.h;
        }
    }
    else{
        u_roe = 0.0;
    }

    return make_tuple(h_roe, u_roe);
}


void Edge::fluxCalc(double flux_dir){
    if (boundaryFlag==0){
        fluxCalc_hllc(flux_dir);
    }
    else if (boundaryFlag<=5){
        fluxCalc_BC(flux_dir);
    }
    else if (boundaryFlag==6){
        fluxCalc_double_solid(flux_dir);
    }
    else{
        cout << "Unknown Boundary Condition" << endl;
    }
}

void Edge::fluxCalc_BC(double flux_dir){
    double uR_rot;
    U UR_rot;

    // determine the direction of the outward vector normal to the boundary
    double d_theta = (flux_dir==1) ? 0.0 : -pi;
    double theta_edgeN = theta_n + d_theta;

    for (int i = 0; i < F_GP.size(); i++){
        UR_rot = rotate(UR_GP[i], theta_edgeN);
        uR_rot = (UR_rot.h==0.0) ? 0.0 : UR_rot.hu/UR_rot.h;
        F_GP[i].F1 = UR_GP[i].h * uR_rot;
        F_GP[i].F2 = UR_GP[i].hu * uR_rot + 0.5*g*(pow(UR_GP[i].h, 2) - pow(zm, 2))*cos(theta_edgeN);
        F_GP[i].F3 = UR_GP[i].hv * uR_rot + 0.5*g*(pow(UR_GP[i].h, 2) - pow(zm, 2))*sin(theta_edgeN);
    }
}

void Edge::fluxCalc_double_solid(double flux_dir){
    int nGP = UR_GP.size() / 2;
    // flux_dir = 1 for owner cell flux; flux_dir = -1 for neighbour cell flux
    int shift = (flux_dir==1) ? 0 : nGP;

    // determine the direction of the outward vector normal to the boundary
    double d_theta = (flux_dir==1) ? 0.0 : -pi;
    double theta_edgeN = theta_n + d_theta;
    
    for (int i = 0; i < F_GP.size(); i++){
        F_GP[i].F1 = 0.0;
        F_GP[i].F2 = 0.5*g*(pow(UR_GP[i+shift].h, 2) - pow(zm, 2))*cos(theta_edgeN);
        F_GP[i].F3 = 0.5*g*(pow(UR_GP[i+shift].h, 2) - pow(zm, 2))*sin(theta_edgeN);
    }
}

// HLLC solver
void Edge::fluxCalc_hllc(double flux_dir){
    double h_roe, u_roe;
    double S_star, SL, SR;
    double uL_rot, vL_rot, uR_rot, vR_rot;
    U UL_local, UL_GP_rot, UR_local,UR_GP_rot;
    FluxF FL_rot, FR_rot, F_hll;

    // determine the direction of the outward vector normal to the boundary
    double d_theta = (flux_dir==1) ? 0.0 : -pi;
    double theta_edgeN = theta_n + d_theta;
	
    for (int i = 0; i < UL_GP.size(); i++){
        if (flux_dir==1){
            UL_local = UL_GP[i];
            UR_local = UR_GP[i];
        }
        else{
            UL_local = UR_GP[i];
            UR_local = UL_GP[i];
        }
        UL_GP_rot = rotate(UL_local, theta_edgeN);
        UR_GP_rot = rotate(UR_local, theta_edgeN);

        // Wave speed estimation
        if (pow(UL_GP_rot.h, 2) + pow(UR_GP_rot.h, 2) > 0){
            if (UL_GP_rot.h * UR_GP_rot.h > 0){
                tie(h_roe, u_roe) = RoeAvg(UL_GP_rot, UR_GP_rot);
                uL_rot = UL_GP_rot.hu / UL_GP_rot.h;
                uR_rot = UR_GP_rot.hu / UR_GP_rot.h;
                SL = min(uL_rot-sqrt(g*UL_GP_rot.h), u_roe-sqrt(g*h_roe));
                SR = max(uR_rot+sqrt(g*UR_GP_rot.h), u_roe+sqrt(g*h_roe));
            }
            else if (UR_GP_rot.h>0.0){
                uR_rot = UR_GP_rot.hu / UR_GP_rot.h;
                SL = uR_rot - 2*sqrt(g*UR_GP_rot.h);
                SR = max(uR_rot+sqrt(g*UR_GP_rot.h), u_roe+sqrt(g*h_roe));
            }
            else{
                uL_rot = UL_GP_rot.hu / UL_GP_rot.h;
                SL = min(uL_rot-sqrt(g*UL_GP_rot.h), u_roe-sqrt(g*h_roe));
                SR = uL_rot + 2*sqrt(g*UL_GP_rot.h);
            }
            S_star = (SL*UR_GP_rot.hu - SL*SR*UR_GP_rot.h - SR*UL_GP_rot.hu + SL*SR*UL_GP_rot.h) / (UR_GP_rot.hu - SR*UR_GP_rot.h - UL_GP_rot.hu + SL*UL_GP_rot.h);
        }
        else{
            SL = 0.0;
            SR = 0.0;
            S_star = 0.0;
        }
        
        // HLL flux calculation
        FL_rot = get_F(UL_GP_rot);
        FR_rot = get_F(UR_GP_rot);

        if (SR != SL){
            F_hll = (SR*FL_rot - SL*FR_rot + SL*SR*toF(toVar(UR_GP_rot-UL_GP_rot))) / (SR - SL);
        }
        else{
            F_hll.F1 = 0.0;
            F_hll.F2 = 0.0;
            F_hll.F3 = 0.0;
        }

        // HLLC flux calculation
        if (SL >= 0){
            F_GP[i] = rotate_back(FL_rot, theta_edgeN);
        }
        else if (S_star >= 0){
            F_GP[i].F1 = F_hll.F1;
            vL_rot = (UL_GP_rot.h == 0) ? 0.0 : UL_GP_rot.hv / UL_GP_rot.h;
            F_GP[i].F2 = cos(theta_edgeN)*F_hll.F2 - vL_rot*sin(theta_edgeN)*F_hll.F1;
            F_GP[i].F3 = sin(theta_edgeN)*F_hll.F2 + vL_rot*cos(theta_edgeN)*F_hll.F1;
        }
        else if (SR > 0){
            F_GP[i].F1 = F_hll.F1;
            vR_rot = (UR_GP_rot.h == 0) ? 0.0 : UR_GP_rot.hv / UR_GP_rot.h;
            F_GP[i].F2 = cos(theta_edgeN)*F_hll.F2 - vR_rot*sin(theta_edgeN)*F_hll.F1;
            F_GP[i].F3 = sin(theta_edgeN)*F_hll.F2 + vR_rot*cos(theta_edgeN)*F_hll.F1;
        }
        else{
            F_GP[i] = rotate_back(FR_rot, theta_edgeN);
        }
    }    
}


Var Cell::getResFlux_MUSCL(const U& u_tmp){
    int eId;
    double l;
    double flux_dir;
    Var res;

    res.v1 = 0.0;
    res.v2 = 0.0;
    res.v3 = 0.0;
    
    for (int i = 0; i < edgeId_Cell.size(); i++){
        eId = edgeId_Cell[i];
        l = globalEdge[eId-1].length;
        // determine the flux direction
        if (globalEdge[eId-1].getOCell()==cId){
            flux_dir = 1.0;
        }
        else{
            flux_dir = -1.0;
        }
        // calculate the flux at the edge
        globalEdge[eId-1].fluxCalc(flux_dir);
        // add to the residual
        // For MUSCL reconstruction, we only have 1 point at the middle of the edge
        res.v1 += -1.0/area * l * globalEdge[eId-1].F_GP[0].F1;
        res.v2 += -1.0/area * l * globalEdge[eId-1].F_GP[0].F2;
        res.v3 += -1.0/area * l * globalEdge[eId-1].F_GP[0].F3;
    }

    return res;
}


Var Cell::getRes_VarRec(const U& u_tmp){
    int eId;
    double w, l;
    double flux_dir;
    
    Var res;
    res.v1 = 0.0;
    res.v2 = 0.0;
    res.v3 = 0.0;

    double u = (u_tmp.h == 0.0) ? 0.0 : u_tmp.hu / u_tmp.h;
    double v = (u_tmp.h == 0.0) ? 0.0 : u_tmp.hv / u_tmp.h;
    
    // convection flux term
    for (int i = 0; i < edgeId_Cell.size(); i++){
        eId = edgeId_Cell[i];
        // determine the flux direction
        if (globalEdge[eId-1].getOCell()==cId){
            flux_dir = 1.0;
        }
        else{
            flux_dir = -1.0;
        }
        // calculate the flux at the edge
        globalEdge[eId-1].fluxCalc(flux_dir);
        // add to the residual
        l = globalEdge[eId-1].length;
        for (int j = 0; j < globalEdge[eId-1].F_GP.size(); j++){
            res.v1 += -1.0/area * l * 0.5 * LineGQ_table.at(vfv_order)[j][0] * globalEdge[eId-1].F_GP[j].F1;
            res.v2 += -1.0/area * l * 0.5 * LineGQ_table.at(vfv_order)[j][0] * globalEdge[eId-1].F_GP[j].F2;
            res.v3 += -1.0/area * l * 0.5 * LineGQ_table.at(vfv_order)[j][0] * globalEdge[eId-1].F_GP[j].F3;
        }
    }

     // source term
    res.v1 += S_p;
    res.v2 += -g*(u_tmp.h+cpz)*ix + n*n*u*sqrt(u*u+v*v)/pow(u_tmp.h, 4.0/3.0);
    res.v3 += -g*(u_tmp.h+cpz)*iy + n*n*v*sqrt(u*u+v*v)/pow(u_tmp.h, 4.0/3.0);

    return res;
}


// ********************************* Wetting and Drying Treatment *********************************
int Edge::isFlood(int cId){
    int nGP, shift;
    int _FloodFlag = 0;

    // inner edge
    if (boundaryFlag==0){
        for (auto i = UL_GP.begin(); i != UL_GP.end(); i++){
            if ((*i).h > h_FLOOD){
                _FloodFlag = 1;
                return _FloodFlag;
            }
        }
        for (auto i = UR_GP.begin(); i != UR_GP.end(); i++){
            if ((*i).h > h_FLOOD){
                _FloodFlag = 1;
                return _FloodFlag;
            }
        }
    }
    // single face boundary
    else if (boundaryFlag<=5){
        for (auto i = UR_GP.begin(); i != UR_GP.end(); i++){
            if ((*i).h > h_FLOOD){
                _FloodFlag = 1;
                return _FloodFlag;
            }
        }
    }
    // double face boundary
    else{
        nGP = UR_GP.size() / 2;
        shift = (cellId_Edge[0]==cId) ? 0 : nGP;
        for (int i = 0; i < nGP; i++){
            if (UR_GP[i+shift].h > h_FLOOD){
                _FloodFlag = 1;
                return _FloodFlag;
            }
        }
    }

    return _FloodFlag;
}

void Cell::FloodEdgeCount(){
    int eId;
    _FloodedEdgeCount = 0;

    for (int i = 0; i < edgeId_Cell.size(); i++){
        eId = edgeId_Cell[i];
        if (globalEdge[eId-1].isFlood(cId)){
            _FloodedEdgeCount += 1;
        }
    }
}

void Cell::WetDryCheck(const U& u_tmp){
    int ncId;

    // determine the type of cell 
    // according to the number of flooded edges and free water surface elevation
    FloodEdgeCount();

    if ((_FloodedEdgeCount == 3) || (u_tmp.eta > z_node[2])){
        label = 1;   // fully submerged cell
    }
    else if ((_FloodedEdgeCount == 0) && (u_tmp.h < h_DRY)){
        label = 2;   // dry cell
    }
    else if (_FloodedEdgeCount >= 1){
        label = 3;   // partially submerged cell
    }
    else{
        label = 4;   // other cell
    }
}


// If the cell is dry, only update the mass
void Cell::update_Dry(const U& u_rhs, int alpha){
    _U_alpha[alpha+1].h = u_rhs.h;
    _U_alpha[alpha+1].hu = _U_alpha[alpha].hu;
    _U_alpha[alpha+1].hv = _U_alpha[alpha].hv;
    _U_alpha[alpha+1] = VFR_set_eta(_U_alpha[alpha+1]);
}

// if the cell is semi-flooded, only update the mass
void Cell::update_SemiFlooded(const U& u_rhs, int alpha){
    if (_U_alpha[alpha+1].h + u_rhs.h <= 0.0){
        _U_alpha[alpha+1].h = 0.0;
        _U_alpha[alpha+1].hu = 0.0;
        _U_alpha[alpha+1].hv = 0.0;
    }
    else{
        _U_alpha[alpha+1].h += u_rhs.h;
        _U_alpha[alpha+1].hu = _U_alpha[alpha].hu;
        _U_alpha[alpha+1].hv = _U_alpha[alpha].hv;   
    }
    _U_alpha[alpha+1] = VFR_set_eta(_U_alpha[alpha+1]);
}

// if all of the cell's edges are dry but the cell water depth >= h_DRY,
// distribute its water to neighbouring semi-flooded cell
void Cell::update_Others(int alpha){
    int ncId;
    double av_water, water_d;

    av_water = _U_alpha[alpha].h * area;

    int semi_count = 0;
   
    for (int i = 0; i < neighborCellId_Cell.size(); i++){
        ncId = neighborCellId_Cell[i];
        if ((ncId >= 1) && (globalCell[ncId-1].label==3)){
            semi_count += 1;
        }
    }

    if (semi_count > 0){
        water_d = av_water/semi_count;
        for (int i = 0; i < neighborCellId_Cell.size(); i++){
            ncId = neighborCellId_Cell[i];
            if ((ncId >= 1) && (globalCell[ncId-1].label==3)){
                globalCell[ncId-1].increValh_alpha(water_d/globalCell[ncId-1].area, alpha+1);
            }
        }
        _U_alpha[alpha+1].h = 0.0;
        _U_alpha[alpha+1].hu = 0.0;
        _U_alpha[alpha+1].hv = 0.0;
        _U_alpha[alpha+1] = VFR_set_eta(_U_alpha[alpha+1]);
    }
    else{
        _U_alpha[alpha+1].h = _U_alpha[alpha].h;
        _U_alpha[alpha+1].hu = _U_alpha[alpha].hu;
        _U_alpha[alpha+1].hv = _U_alpha[alpha].hv;
        _U_alpha[alpha+1] = VFR_set_eta(_U_alpha[alpha+1]);
    }
}


void Cell::update_Others_Uc(){
    int ncId;
    double av_water, water_d;

    av_water = _Uc.h * area;

    int semi_count = 0;
   
    for (int i = 0; i < neighborCellId_Cell.size(); i++){
        ncId = neighborCellId_Cell[i];
        if ((ncId >= 1) && (globalCell[ncId-1].label==3)){
            semi_count += 1;
        }
    }

    if (semi_count > 0){
        water_d = av_water/semi_count;
        for (int i = 0; i < neighborCellId_Cell.size(); i++){
            ncId = neighborCellId_Cell[i];
            if ((ncId >= 1) && (globalCell[ncId-1].label==3)){
                globalCell[ncId-1].increValh_Uc(water_d/globalCell[ncId-1].area);
            }
        }
        _Uc.h = 0.0;
        _Uc.hu = 0.0;
        _Uc.hv = 0.0;
        _Uc = VFR_set_eta(_Uc);
    }
}


// ********************************* Time Stepping *********************************
// [1] RK3 stepping
double Cell::getCellLocalTime(double dt){
    int eId, ncId;
    double l;
    double lambda_l = 0.0;

    for (int i = 0; i < edgeId_Cell.size(); i++){
        eId = edgeId_Cell[i];
        l = globalEdge[eId-1].length;
        ncId = neighborCellId_Cell[i];
        if (ncId >= 1){
            lambda_l = max(l*globalEdge[eId-1].getEdgeSpectralRadius(_Uc, globalCell[ncId-1].getValUc()), lambda_l);
        }
    }

    /*if ((lambda_l > 0.0) && (cfl * area / lambda_l < 1e-6)){
        cout << cId << "\t" << _Uc.h << endl;
    }*/
    
    if (lambda_l == 0.0){
        return dt;
    }
    else{
        return cfl * area / lambda_l;
    }
}

void Cell::update_RK3(double dt, int alpha){
    U u_rhs;
    double u, v, u_p1, v_p1;
    //double q_hat, S_fx, S_fy, Dx, Dy;
    //double qx_grad, qy_grad;

    /*if ((cId==2552)){
        cout << cId << "\t" << alpha << endl;
    }*/

    double a1 = RK3_table[alpha][0];
    double a2 = RK3_table[alpha][1];
    double a3 = RK3_table[alpha][2];

    u_rhs = a1 * _Uc + a2 * _U_alpha[alpha];
    
    // Add net rainfall and slope source term
    u_rhs.h += a3 * dt * S_p;
    u_rhs.hu += -a3 * dt * g*(_U_alpha[alpha].h+cpz)*ix;
    u_rhs.hv += -a3 * dt * g*(_U_alpha[alpha].h+cpz)*iy;

    // fully submerged cell
    if (label==1){
        // Add convection flux term
        u_rhs = u_rhs + a3 * dt * getResFlux_MUSCL(_U_alpha[alpha]);
        // solve semi-implicitly
        if (u_rhs.h <= 0){
            _U_alpha[alpha+1].h = 0.0;
            _U_alpha[alpha+1].hu = 0.0;
            _U_alpha[alpha+1].hv = 0.0;
            _U_alpha[alpha+1] = VFR_set_eta(_U_alpha[alpha+1]);
        }
        else{
            /*u = u_rhs.hu / u_rhs.h;
            v = u_rhs.hv / u_rhs.h;
            q_hat = sqrt(pow(u_rhs.hu, 2)+pow(u_rhs.hv, 2));
            S_fx = g * pow(n, 2) * u * sqrt(pow(u, 2)+pow(v, 2)) / pow(u_rhs.h, 1.0/3.0);
            S_fy = g * pow(n, 2) * v * sqrt(pow(u, 2)+pow(v, 2)) / pow(u_rhs.h, 1.0/3.0);
            Dx = 1.0 + a3 * dt * g * pow(n, 2) / pow(u_rhs.h, 7.0/3.0) * (q_hat + pow(u_rhs.hu, 2)/q_hat);
            Dy = 1.0 + a3 * dt * g * pow(n, 2) / pow(u_rhs.h, 7.0/3.0) * (q_hat + pow(u_rhs.hv, 2)/q_hat);
            if (u_rhs.hu >= 0.0){
                qx_grad = max(-u_rhs.hu/(a3*dt), S_fx/Dx);
            }
            else{
                qx_grad = min(-u_rhs.hu/(a3*dt), S_fx/Dx);
            }
            if (u_rhs.hv>= 0.0){
                qy_grad = max(-u_rhs.hv/(a3*dt), S_fy/Dy);
            }
            else{
                qy_grad = min(-u_rhs.hv/(a3*dt), S_fy/Dy);
            }
            _U_alpha[alpha+1].h = u_rhs.h;
            _U_alpha[alpha+1].hu = u_rhs.hu + a3 * dt * qx_grad;
            _U_alpha[alpha+1].hv = u_rhs.hv + a3 * dt * qy_grad;
            _U_alpha[alpha+1] = VFR_set_eta(_U_alpha[alpha+1]);
            */

            u = u_rhs.hu / u_rhs.h;
            v = u_rhs.hv / u_rhs.h;
            u_p1 = u / (1.0 + a3*dt*g*n*n*sqrt(u*u+v*v)/pow(u_rhs.h, 4.0/3.0));
            v_p1 = v / (1.0 + a3*dt*g*n*n*sqrt(u*u+v*v)/pow(u_rhs.h, 4.0/3.0));

            _U_alpha[alpha+1].h = u_rhs.h;
            _U_alpha[alpha+1].hu = u_rhs.h * u_p1;
            _U_alpha[alpha+1].hv = u_rhs.h * v_p1;
            _U_alpha[alpha+1] = VFR_set_eta(_U_alpha[alpha+1]);
        }
    }
    // dry cell
    else if (label==2){
        update_Dry(u_rhs, alpha);
    }
    // partially submerged cell
    else if (label==3){
        // Add convection flux term
        u_rhs = u_rhs + a3 * dt * getResFlux_MUSCL(_U_alpha[alpha]);
        update_SemiFlooded(u_rhs, alpha);
    }
    // other cell
    else{
        update_Others(alpha);
    }
}

void Cell::update_RK3_Step(){
    _Uc.h = _U_alpha[3].h;
    _Uc.hu = _U_alpha[3].hu;
    _Uc.hv = _U_alpha[3].hv;
    _Uc = VFR_set_eta(_Uc);

    for (int i = 0; i < _U_alpha.size(); i++){
        _U_alpha[i].h = 0.0;
        _U_alpha[i].hu = 0.0;
        _U_alpha[i].hv = 0.0;
    }
}

// [2] SDIRK4 stepping
void Cell::update_SDIRK4_InnerUpdate(){
    if (_u_inner.h + _dU.v1 <= 0.0){
        _u_inner.h = 0.0;
        _u_inner.hu = 0.0;
        _u_inner.hv = 0.0;
    }
    else{
        _u_inner.h += _dU.v1;
        _u_inner.hu += _dU.v2;
        _u_inner.hv += _dU.v3;
    }
    
    _u_inner = VFR_set_eta(_u_inner);

    _dU.v1 = 0.0;
    _dU.v2 = 0.0;
    _dU.v3 = 0.0;
}


void Cell::update_SDIRK4_Step(double dt){
    if (label==1){
        for (int a = 0; a < _U_alpha.size(); a++){
            _Uc = _Uc + dt * SDIRK4_weight[a] * getRes_VarRec(_U_alpha[a]);
        }
    }
    // dry cell
    else if (label==2){
        _Uc.h = _Uc.h + dt * S_p;
    }
    // partially submerged cell
    else if (label==3){
        for (int a = 0; a < _U_alpha.size(); a++){
            if (_Uc.h + dt * SDIRK4_weight[a] * getRes_VarRec(_U_alpha[a]).v1 <= 0.0){
                _Uc.h = 0.0;
                _Uc.hu = 0.0;
                _Uc.hv = 0.0;
            }
            else{
                _Uc.h = _Uc.h + dt * SDIRK4_weight[a] * getRes_VarRec(_U_alpha[a]).v1;
            }
        }
    }
    // other cell
    else{
        update_Others_Uc();
    }
    // Update the free surface level with new h using VFR
    _Uc = VFR_set_eta(_Uc);
}


// LU-SGS: forward sweep
void Cell::update_LUSGS_forward(double dt, double dtau, int alpha){
    Var rhs;
    U U_nc, U_bd;
    Var dU_nc;
    double h_nc, u_nc, v_nc;
    int eId, ncId;
    double l, lambda, cc, theta, nx, ny;
    
    double J_11, J_12, J_13, J_21, J_22, J_23, J_31, J_32, J_33;
    double lambda_l = 0.0;
    J_11 = 0.0;
    rhs.v1 = 0.0;
    rhs.v2 = 0.0;
    rhs.v3 = 0.0;

    /*if (cId == 61){
        cout << cId << "\t" << rhs.v1 << "\t" << rhs.v2 << "\t" << rhs.v3 << endl;
    }*/
    
    rhs = toVar(_Uc - _u_inner)/dt + SDIRK4_table.at(alpha)[alpha] * getRes_VarRec(_u_inner);
    for (int beta = 0; beta < alpha; beta++){
        rhs = rhs + SDIRK4_table.at(alpha)[beta] * _Res_alpha[beta];
    }

    for (int i = 0; i < edgeId_Cell.size(); i++){
        ncId = neighborCellId_Cell[i];
        eId = edgeId_Cell[i];
        l = globalEdge[eId-1].length;
        theta = globalEdge[eId-1].theta_n;
        nx = cos(theta);
        ny = sin(theta);
        // if neighbour cell exists
        if (edgeBCTag_Cell[i] == 0){
            dU_nc = globalCell[ncId-1].getValdU();
            // calculate the SpectralRadius
            U_nc = globalCell[ncId-1].getValUc();
            h_nc = U_nc.h;
            u_nc = (U_nc.h==0.0) ? 0.0 : U_nc.hu / U_nc.h;
            v_nc = (U_nc.h==0.0) ? 0.0 : U_nc.hv / U_nc.h;
            lambda = globalEdge[eId-1].getEdgeSpectralRadius(_Uc, U_nc);
            cc = SDIRK4_table.at(alpha)[alpha]/area*l*0.5;
            // Forward sweep only considers the neighbour cell whose index is smaller
            if (ncId < cId){
                // calculate the Jacobian of F for the neighbour cell
                J_12 = nx;
                J_13 = ny;
                J_21 = (g * h_nc - u_nc * u_nc) * nx - u_nc * v_nc * ny;
                J_22 = 2 * u_nc * nx + v_nc * ny;
                J_23 = u_nc * ny;
                J_31 = -u_nc * v_nc * nx + (g * h_nc - v_nc * v_nc) * ny;
                J_32 = v_nc * nx;
                J_33 = u_nc * nx + 2 * v_nc * ny;
                // update the residual 
                rhs.v1 += - cc * (J_11*dU_nc.v1 + J_12*dU_nc.v2 + J_13*dU_nc.v3) + cc * lambda*dU_nc.v1;
                rhs.v2 += - cc * (J_21*dU_nc.v1 + J_22*dU_nc.v2 + J_23*dU_nc.v3) + cc * lambda*dU_nc.v2;
                rhs.v3 += - cc * (J_31*dU_nc.v1 + J_32*dU_nc.v2 + J_33*dU_nc.v3) + cc * lambda*dU_nc.v3;
            }
            lambda_l += cc * lambda;
        }
        else{
            U_bd.h = 0.0;
            U_bd.hu = 0.0;
            U_bd.hv = 0.0;
            U_bd.eta = 0.0;
            for (int j = 0; j < globalEdge[eId-1].UR_GP.size(); j++){
                U_bd = U_bd + 0.5 * LineGQ_table.at(vfv_order)[j][0] * globalEdge[eId-1].UR_GP[j];
            }
            lambda = globalEdge[eId-1].getEdgeSpectralRadius(_Uc, U_bd);
            cc = SDIRK4_table.at(alpha)[alpha]/area*l*0.5;
            lambda_l += cc * lambda;
        }
    }
    // set the delta U for the cell
    _dU = rhs / (1.0/dtau + 1.0/dt + lambda_l);
}

// LU-SGS: backward sweep
double Cell::update_LUSGS_backward(double dt, double dtau, int alpha){
    Var rhs;
    U U_nc, U_bd;
    Var dU_nc;
    double h_nc, u_nc, v_nc;
    int eId, ncId;
    double l, lambda, cc, theta, nx, ny;

    double J_11, J_12, J_13, J_21, J_22, J_23, J_31, J_32, J_33;
    double lambda_l = 0.0;
    J_11 = 0.0;
    rhs.v1 = 0.0;
    rhs.v2 = 0.0;
    rhs.v3 = 0.0;

    for (int i = 0; i < edgeId_Cell.size(); i++){
        ncId = neighborCellId_Cell[i];
        eId = edgeId_Cell[i];
        l = globalEdge[eId-1].length;
        theta = globalEdge[eId-1].theta_n;
        nx = cos(theta);
        ny = sin(theta);
        // if neighbour cell exists
        if (edgeBCTag_Cell[i] == 0){
            dU_nc = globalCell[ncId-1].getValdU();
            // calculate the SpectralRadius
            U_nc = globalCell[ncId-1].getValUc();
            h_nc = U_nc.h;
            u_nc = (U_nc.h==0.0) ? 0.0 : U_nc.hu / U_nc.h;
            v_nc = (U_nc.h==0.0) ? 0.0 : U_nc.hv / U_nc.h;
            lambda = globalEdge[eId-1].getEdgeSpectralRadius(_Uc, U_nc);
            cc = SDIRK4_table.at(alpha)[alpha]/area*l*0.5;
            // Backward sweep only considers the neighbour cell whose index is larger
            if (ncId > cId){
                // calculate the Jacobian of F for the neighbour cell
                J_12 = nx;
                J_13 = ny;
                J_21 = (g * h_nc - u_nc * u_nc) * nx - u_nc * v_nc * ny;
                J_22 = 2 * u_nc * nx + v_nc * ny;
                J_23 = u_nc * ny;
                J_31 = -u_nc * v_nc * nx + (g * h_nc - v_nc * v_nc) * ny;
                J_32 = v_nc * nx;
                J_33 = u_nc * nx + 2 * v_nc * ny;
                // update the residual 
                rhs.v1 += cc * (J_11*dU_nc.v1 + J_12*dU_nc.v2 + J_13*dU_nc.v3) - cc * lambda*dU_nc.v1;
                rhs.v2 += cc * (J_21*dU_nc.v1 + J_22*dU_nc.v2 + J_23*dU_nc.v3) - cc * lambda*dU_nc.v2;
                rhs.v3 += cc * (J_31*dU_nc.v1 + J_32*dU_nc.v2 + J_33*dU_nc.v3) - cc * lambda*dU_nc.v3;
            }
            lambda_l += cc * lambda;
        }
        else{
            U_bd.h = 0.0;
            U_bd.hu = 0.0;
            U_bd.hv = 0.0;
            U_bd.eta = 0.0;
            for (int j = 0; j < globalEdge[eId-1].UR_GP.size(); j++){
                U_bd = U_bd + 0.5 * LineGQ_table.at(vfv_order)[j][0] * globalEdge[eId-1].UR_GP[j];
            }
            lambda = globalEdge[eId-1].getEdgeSpectralRadius(_Uc, U_bd);
            cc = SDIRK4_table.at(alpha)[alpha]/area*l*0.5;
            lambda_l += cc * lambda;
        }
    }
    _dU = _dU - rhs / (1.0/dtau + 1.0/dt + lambda_l);

    return _dU.v1;
}