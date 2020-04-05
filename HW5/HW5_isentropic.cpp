#include <math.h>
#include <fstream>
#include <iostream>
#include <string>
#include <array>
#include <tuple>
#include <algorithm>
#include "HW5_isentropic.hpp"


using namespace std;

int id_map(int id, int id_start, int id_end, double T){
    int n1, n2, id_new;

    n1 = floor((id - id_start) / T);
    n2 = floor((id - id_end) / T);

    if (n1 < n2){
        id_new = id - n2*T;
    }
    else {
        id_new = id - n1*T;
    }

    return id_new;
}

double sign(double x){
    if (x>=0.0){
        return 1.0;
    }
    else{
        return -1.0;
    }
}

Var operator+(const Var& a, const Var& b){
    Var var_sum;

    var_sum.v1 = a.v1 + b.v1;
    var_sum.v2 = a.v2 + b.v2;
    var_sum.v3 = a.v3 + b.v3;
    var_sum.v4 = a.v4 + b.v4;

    return var_sum;
}

Var operator-(const Var& a, const Var& b){
    Var var_minus;

    var_minus.v1 = a.v1 - b.v1;
    var_minus.v2 = a.v2 - b.v2;
    var_minus.v3 = a.v3 - b.v3;
    var_minus.v4 = a.v4 - b.v4;

    return var_minus;
}

Var operator*(double alpha, const Var&a){
    Var var_mult;

    var_mult.v1 = alpha * a.v1;
    var_mult.v2 = alpha * a.v2;
    var_mult.v3 = alpha * a.v3;
    var_mult.v4 = alpha * a.v4;

    return var_mult;
}


Matrix<double> operator*(const Matrix<double>& A, const Matrix<double>& B){
    Matrix<double> C(A.ny, B.nx);
    for (int i = 0; i < C.ny; i++){
        for (int j = 0; j < C.nx; j++){
            C[i][j] = 0.0;
        }
    }

    for (int i = 0; i < C.ny; i++){
        for (int j = 0; j < C.nx; j++){
            for (int k = 0; k < A.nx; k++){
                C[i][j] += A[i][k]*B[k][j];
            }
        }
    }
    
    return C;
}

Matrix<double> operator*(double alpha, const Matrix<double>& A){
    Matrix<double> C(A.ny, A.nx);
    for (int i = 0; i < C.ny; i++){
        for (int j = 0; j < C.nx; j++){
            C[i][j] = alpha*A[i][j];
        }
    }
    
    return C;
}

Matrix<double> operator+(const Matrix<double>& A, const Matrix<double>& B){
    Matrix<double> C(A.ny, A.nx);
    for (int i = 0; i < C.ny; i++){
        for (int j = 0; j < C.nx; j++){
            C[i][j] = A[i][j]+B[i][j];
        }
    }
    return C;
}

Matrix<double> operator-(const Matrix<double>& A, const Matrix<double>& B){
    Matrix<double> C(A.ny, A.nx);
    for (int i = 0; i < C.ny; i++){
        for (int j = 0; j < C.nx; j++){
            C[i][j] = A[i][j]-B[i][j];
        }
    }
    return C;
}

Var operator*(const Matrix<double>& A, const Var& B){
    Var var_dot;

    var_dot.v1 = A[0][0]*B.v1 + A[0][1]*B.v2 + A[0][2]*B.v3 + A[0][3]*B.v4;
    var_dot.v2 = A[1][0]*B.v1 + A[1][1]*B.v2 + A[1][2]*B.v3 + A[1][3]*B.v4;
    var_dot.v3 = A[2][0]*B.v1 + A[2][1]*B.v2 + A[2][2]*B.v3 + A[2][3]*B.v4;
    var_dot.v4 = A[3][0]*B.v1 + A[3][1]*B.v2 + A[3][2]*B.v3 + A[3][3]*B.v4;

    return var_dot;
}

double get_p(const Var& u){
    double gamma = 1.4;
    double p = (gamma - 1.0) * (u.v4 - 0.5*(u.v2*u.v2+u.v3*u.v3)/u.v1);
    return p;
}

double get_H(const Var& u){
    double gamma = 1.4;
    double p = get_p(u);
    double a = sqrt(gamma*p/u.v1);
    double H = a*a/(gamma-1.0) + 0.5*(u.v2*u.v2+u.v3*u.v3)/(u.v1*u.v1);
    return H;
}

Var get_F(const Var& u){
    Var F;

    double p = get_p(u);
    F.v1 = u.v2;
    F.v2 = u.v2*u.v2/u.v1 + p;
    F.v3 = u.v2*u.v3/u.v1;
    F.v4 = u.v2/u.v1*(u.v4+p);

    return F;
}

Var get_G(const Var& u){
    Var G;

    double p = get_p(u);
    G.v1 = u.v3;
    G.v2 = u.v2*u.v3/u.v1;
    G.v3 = u.v3*u.v3/u.v1 + p;
    G.v4 = u.v3/u.v1*(u.v4+p);

    return G;
}

tuple<double, double, double, double> RoeAvg(const Var& A, const Var& B){
    double u_roe, v_roe, H_roe, a_roe;
    double Ha, Hb;
    const double gamma = 1.4;

    Ha = get_H(A);
    Hb = get_H(B);
    u_roe = (A.v2/sqrt(A.v1)+B.v2/sqrt(B.v1)) / (sqrt(A.v1)+sqrt(B.v1)); 
    v_roe = (A.v3/sqrt(A.v1)+B.v3/sqrt(B.v1)) / (sqrt(A.v1)+sqrt(B.v1));
    H_roe = (sqrt(A.v1)*Ha+sqrt(B.v1)*Hb) / (sqrt(A.v1)+sqrt(B.v1));
    a_roe = sqrt((gamma-1.0)*(H_roe-0.5*(u_roe*u_roe+v_roe*v_roe)));
    return make_tuple(u_roe, v_roe, H_roe, a_roe);
}

// Left EigenMatrix used in Roe method
Matrix<double> get_RoeLEM(double u_roe, double v_roe, double H_roe, double a_roe, double nx, double ny){
    Matrix<double> LeftEigenMat(4, 4);
    const double gamma = 1.4;
    double qn = u_roe*nx + v_roe*ny;
    double lx = -ny;
    double ly = nx;
    double ql = u_roe*lx + v_roe*ly;

    // -----Formation of Left Eigenvectors
    LeftEigenMat[0][0] = 0.5*((gamma-1.0)/(2*a_roe*a_roe)*(u_roe*u_roe+v_roe*v_roe)+qn/a_roe);
    LeftEigenMat[0][1] = -0.5*((gamma-1.0)/(a_roe*a_roe)*u_roe+nx/a_roe);
    LeftEigenMat[0][2] = -0.5*((gamma-1.0)/(a_roe*a_roe)*v_roe+ny/a_roe);
    LeftEigenMat[0][3] = (gamma-1.0)/(2*a_roe*a_roe);
    LeftEigenMat[1][0] = 1.0 - (gamma-1.0)/(2*a_roe*a_roe)*(u_roe*u_roe+v_roe*v_roe);
    LeftEigenMat[1][1] = (gamma-1.0)/(a_roe*a_roe)*u_roe;
    LeftEigenMat[1][2] = (gamma-1.0)/(a_roe*a_roe)*v_roe;
    LeftEigenMat[1][3] = -(gamma-1.0)/(a_roe*a_roe);
    LeftEigenMat[2][0] = 0.5*((gamma-1.0)/(2*a_roe*a_roe)*(u_roe*u_roe+v_roe*v_roe)-qn/a_roe);
    LeftEigenMat[2][1] = -0.5*((gamma-1.0)/(a_roe*a_roe)*u_roe-nx/a_roe);
    LeftEigenMat[2][2] = -0.5*((gamma-1.0)/(a_roe*a_roe)*v_roe-ny/a_roe);
    LeftEigenMat[2][3] = (gamma-1.0)/(2*a_roe*a_roe);
    LeftEigenMat[3][0] = -ql;
    LeftEigenMat[3][1] = lx;
    LeftEigenMat[3][2] = ly;
    LeftEigenMat[3][3] = 0.0;
    
    return LeftEigenMat;
}

// Right EigenMatrix used in Roe method
Matrix<double> get_RoeREM(double u_roe, double v_roe, double H_roe, double a_roe, double nx, double ny){
    Matrix<double> RightEigenMat(4, 4);
    const double gamma = 1.4;
    double qn = u_roe*nx + v_roe*ny;
    double lx = -ny;
    double ly = nx;
    double ql = u_roe*lx + v_roe*ly;

    // -----Formation of Right Eigenvectors
    RightEigenMat[0][0] = 1.0;
    RightEigenMat[0][1] = 1.0;
    RightEigenMat[0][2] = 1.0;
    RightEigenMat[0][3] = 0.0;
    RightEigenMat[1][0] = u_roe - a_roe*nx;
    RightEigenMat[1][1] = u_roe;
    RightEigenMat[1][2] = u_roe + a_roe*nx;
    RightEigenMat[1][3] = lx;
    RightEigenMat[2][0] = v_roe - a_roe*ny;
    RightEigenMat[2][1] = v_roe;
    RightEigenMat[2][2] = v_roe + a_roe*ny;
    RightEigenMat[2][3] = ly;
    RightEigenMat[3][0] = H_roe - a_roe*qn;
    RightEigenMat[3][1] = 0.5*(u_roe*u_roe+v_roe*v_roe);
    RightEigenMat[3][2] = H_roe + a_roe*qn;
    RightEigenMat[3][3] = ql;

    return RightEigenMat;
}

Matrix<double> get_SignEigen(double u_roe, double v_roe, double a_roe, double nx, double ny){
    Matrix<double> SignEigenMat(4, 4);
    double qn = u_roe*nx + v_roe*ny;

    SignEigenMat[0][0] = sign(qn-a_roe);
    SignEigenMat[0][1] = 0.0;
    SignEigenMat[0][2] = 0.0;
    SignEigenMat[0][3] = 0.0;
    SignEigenMat[1][0] = 0.0;
    SignEigenMat[1][1] = sign(qn);
    SignEigenMat[1][2] = 0.0;
    SignEigenMat[1][3] = 0.0;
    SignEigenMat[2][0] = 0.0;
    SignEigenMat[2][1] = 0.0;
    SignEigenMat[2][2] = sign(qn+a_roe);
    SignEigenMat[2][3] = 0.0;
    SignEigenMat[3][0] = 0.0;
    SignEigenMat[3][1] = 0.0;
    SignEigenMat[3][2] = 0.0;
    SignEigenMat[3][3] = sign(qn);

    return SignEigenMat;
}


void gradCalc_2rdJameson(Matrix<Var>& sol, Matrix<Var>& FluxX, Matrix<Var>& FluxY, double dx, double dy, int nx, int ny, const string scheme){
    // x,y for stencil points
    int im2, im1, i0, ip1, ip2, ip3;
    int jm2, jm1, j0, jp1, jp2, jp3;
    // pressure
    double p_im2, p_im1, p_i0, p_ip1, p_ip2, p_ip3;
    double p_jm2, p_jm1, p_j0, p_jp1, p_jp2, p_jp3;
    // viscosity coefficient
    double v_im1, v_i0, v_ip1, v_ip2;
    double v_jm1, v_j0, v_jp1, v_jp2;
    double epsilon_2, epsilon_4;
    // wave speed
    double a_i0, a_ip1;
    double a_j0, a_jp1;
    double lambda_max;
    // flux
    Var g_i0, g_ip1;
    Var f_j0, f_jp1;
    // constant
    const double gamma = 1.4;
    for (int i = 0; i < ny; i++){
        // y coordinate of stencil points 
        im2 = id_map(i-2, 0, ny-1, ny);
        im1 = id_map(i-1, 0, ny-1, ny);
        i0 = id_map(i, 0, ny-1, ny);
        ip1 = id_map(i+1, 0, ny-1, ny);
        ip2 = id_map(i+2, 0, ny-1, ny);
        ip3 = id_map(i+3, 0, ny-1, ny);
        for (int j = 0; j < nx; j++){
            // x coordinate of stencil points 
            jm2 = id_map(j-2, 0, nx-1, nx);
            jm1 = id_map(j-1, 0, nx-1, nx);
            j0 = id_map(j, 0, nx-1, nx);
            jp1 = id_map(j+1, 0, nx-1, nx);
            jp2 = id_map(j+2, 0, nx-1, nx);
            jp3 = id_map(j+3, 0, nx-1, nx);
            // Determine the numerical flux along the x axis
            p_jm2 = get_p(sol[i][jm2]);
            p_jm1 = get_p(sol[i][jm1]);
            p_j0 = get_p(sol[i][j0]);
            p_jp1 = get_p(sol[i][jp1]);
            p_jp2 = get_p(sol[i][jp2]);
            p_jp3 = get_p(sol[i][jp3]);
            v_jm1 = fabs(p_j0-2*p_jm1+p_jm2)/fabs(p_j0+2*p_jm1+p_jm2);
            v_j0 = fabs(p_jp1-2*p_j0+p_jm1)/fabs(p_jp1+2*p_j0+p_jm1);
            v_jp1 = fabs(p_jp2-2*p_jp1+p_j0)/fabs(p_jp2+2*p_jp1+p_j0);
            v_jp2 = fabs(p_jp3-2*p_jp2+p_jp1)/fabs(p_jp3+2*p_jp2+p_jp1);
            epsilon_2 = 0.6*max({v_jm1, v_j0, v_jp1, v_jp2});
            epsilon_4 = max({0.0, 1.0/64.0-epsilon_2});
            a_j0 = sqrt(gamma*p_j0/sol[i][j].v1);
            a_jp1 = sqrt(gamma*p_jp1/sol[i][jp1].v1);
            lambda_max =  max(fabs(sol[i][j].v2/sol[i][j].v1)+a_j0, fabs(sol[i][jp1].v2/sol[i][jp1].v1)+a_jp1);
            f_j0 = get_F(sol[i][j0]);
            f_jp1 = get_F(sol[i][jp1]);
            FluxX[i][j] = 0.5*(f_j0+f_jp1) - lambda_max*epsilon_2*(sol[i][jp1]-sol[i][j]) 
                        + lambda_max*epsilon_4*(sol[i][jp2] - 3.0*sol[i][jp1] + 3.0*sol[i][j] - sol[i][jm1]);
            // Determine the numerical flux along the y axis
            p_im2 = get_p(sol[im2][j]);
            p_im1 = get_p(sol[im1][j]);
            p_i0 = get_p(sol[i0][j]);
            p_ip1 = get_p(sol[ip1][j]);
            p_ip2 = get_p(sol[ip2][j]);
            p_ip3 = get_p(sol[ip3][j]);
            v_im1 = fabs(p_i0-2*p_im1+p_im2)/fabs(p_i0+2*p_im1+p_im2);
            v_i0 = fabs(p_ip1-2*p_i0+p_im1)/fabs(p_ip1+2*p_i0+p_im1);
            v_ip1 = fabs(p_ip2-2*p_ip1+p_i0)/fabs(p_ip2+2*p_ip1+p_i0);
            v_ip2 = fabs(p_ip3-2*p_ip2+p_ip1)/fabs(p_ip3+2*p_ip2+p_ip1);
            epsilon_2 = 0.6*max({v_im1, v_i0, v_ip1, v_ip2});
            epsilon_4 = max({0.0, 1.0/64.0-epsilon_2});
            a_i0 = sqrt(gamma*p_i0/sol[i][j].v1);
            a_ip1 = sqrt(gamma*p_ip1/sol[ip1][j].v1);
            lambda_max =  max(fabs(sol[i][j].v3/sol[i][j].v1)+a_i0, fabs(sol[ip1][j].v3/sol[ip1][j].v1)+a_ip1);
            g_i0 = get_G(sol[i0][j]);
            g_ip1 = get_G(sol[ip1][j]);
            FluxY[i][j] = 0.5*(g_i0+g_ip1) - lambda_max*epsilon_2*(sol[ip1][j]-sol[i][j]) 
                        + lambda_max*epsilon_4*(sol[ip2][j] - 3.0*sol[ip1][j] + 3.0*sol[i][j] - sol[im1][j]);
        }
    }
}


void gradCalc_2rdRoe(Matrix<Var>& sol, Matrix<Var>& FluxX, Matrix<Var>& FluxY, double dx, double dy, int nx, int ny, const string scheme){
    // x,y for stencil points
    int im1, i0, ip1, ip2;
    int jm1, j0, jp1, jp2;
    // flux
    Var g_im1, g_i0, g_ip1, g_ip2;
    Var f_jm1, f_j0, f_jp1, f_jp2;
    // Roe average variable
    double u_roe, v_roe, H_roe, a_roe;
    // Left and right eigenvectors
    Matrix<double> RightEigenMat(4, 4);
    Matrix<double> LeftEigenMat(4, 4);
    Matrix<double> SignEigenMat(4, 4);
    /*Matrix<double> IMat(4, 4);
    for (int i = 0; i < 4; i++){
        for (int j = 0; j < 4; j++){
            if (i==j){
                IMat[i][j] = 1.0;
            }
            else{
                IMat[i][j] = 0.0;
            }
        }
    }*/
    // constant
    const double gamma = 1.4;

    for (int i = 0; i < ny; i++){
        // y coordinate of stencil points 
        im1 = id_map(i-1, 0, ny-1, ny);
        i0 = id_map(i, 0, ny-1, ny);
        ip1 = id_map(i+1, 0, ny-1, ny);
        ip2 = id_map(i+2, 0, ny-1, ny);
        for (int j = 0; j < nx; j++){
            // x coordinate of stencil points 
            jm1 = id_map(j-1, 0, nx-1, nx);
            j0 = id_map(j, 0, nx-1, nx);
            jp1 = id_map(j+1, 0, nx-1, nx);
            jp2 = id_map(j+2, 0, nx-1, nx);
            // Determine the numerical flux along the x axis
            f_jm1 = get_F(sol[i][jm1]);
            f_j0 = get_F(sol[i][j0]);
            f_jp1 = get_F(sol[i][jp1]);
            f_jp2 = get_F(sol[i][jp2]);
            // -----Roe average
            tie(u_roe, v_roe, H_roe, a_roe) = RoeAvg(sol[i][j0], sol[i][jp1]);
            // -----Formation of Eigenvectors
            RightEigenMat = get_RoeREM(u_roe, v_roe, H_roe, a_roe, 1.0, 0.0);
            LeftEigenMat = get_RoeLEM(u_roe, v_roe, H_roe, a_roe, 1.0, 0.0);
            SignEigenMat = get_SignEigen(u_roe, v_roe, a_roe, 1.0, 0.0);
            /*FluxX[i][j] = 0.5*(RightEigenMat*(IMat+SignEigenMat)*LeftEigenMat)*(1.5*f_j0-0.5*f_jm1) 
                        + 0.5*(RightEigenMat*(IMat-SignEigenMat)*LeftEigenMat)*(1.5*f_jp1-0.5*f_jp2);*/
            FluxX[i][j] = 0.25*(-1.0*f_jm1+3.0*f_j0+3.0*f_jp1-f_jp2) 
                        - 0.25*RightEigenMat*SignEigenMat*LeftEigenMat*(f_jm1-3.0*f_j0+3.0*f_jp1-f_jp2);  
            // Determine the numerical flux along the y axis
            g_im1 = get_G(sol[im1][j]);
            g_i0 = get_G(sol[i0][j]);
            g_ip1 = get_G(sol[ip1][j]);
            g_ip2 = get_G(sol[ip2][j]);
            // -----Roe average
            tie(u_roe, v_roe, H_roe, a_roe) = RoeAvg(sol[i0][j], sol[ip1][j]);
            // -----Formation of Eigenvectors
            RightEigenMat = get_RoeREM(u_roe, v_roe, H_roe, a_roe, 0.0, 1.0);
            LeftEigenMat = get_RoeLEM(u_roe, v_roe, H_roe, a_roe, 0.0, 1.0);
            SignEigenMat = get_SignEigen(u_roe, v_roe, a_roe, 0.0, 1.0);
            /*FluxY[i][j] = 0.5*(RightEigenMat*(IMat+SignEigenMat)*LeftEigenMat)*(1.5*g_i0-0.5*g_im1) 
                        + 0.5*(RightEigenMat*(IMat-SignEigenMat)*LeftEigenMat)*(1.5*g_ip1-0.5*g_ip2);*/
            FluxY[i][j] = 0.25*(-1.0*g_im1+3.0*g_i0+3.0*g_ip1-g_ip2) 
                        - 0.25*RightEigenMat*SignEigenMat*LeftEigenMat*(g_im1-3.0*g_i0+3.0*g_ip1-g_ip2);  
        }
    }
}

void step(Matrix<Var>& sol, Matrix<Var>& FluxX, Matrix<Var>& FluxY, double dx, double dy, double dt, int nx, int ny, const string scheme){
    int im1, jm1;
    Matrix<Var> sol_tmp(nx, ny);

    if (scheme=="2rdJameson"){
        gradCalc_2rdJameson(sol, FluxX, FluxY, dx, dy, nx, ny, scheme);
    }
    else if (scheme=="2rdRoe"){
        gradCalc_2rdRoe(sol, FluxX, FluxY, dx, dy, nx, ny, scheme);
    }
    else {
        cout << "Unknown scheme" << endl;
    }

    for (int i = 0; i < ny; i++){
        for (int j = 0; j < nx; j++){
            im1 = id_map(i-1, 0, ny-1, ny);
            jm1 = id_map(j-1, 0, nx-1, nx);
            sol_tmp[i][j] = sol[i][j] - dt/dx*(FluxX[i][j] - FluxX[i][jm1]) - dt/dy*(FluxY[i][j] - FluxY[im1][j]);
        }
    }

    if (scheme=="2rdJameson"){
        gradCalc_2rdJameson(sol_tmp, FluxX, FluxY, dx, dy, nx, ny, scheme);
    }
    else if (scheme=="2rdRoe"){
        gradCalc_2rdRoe(sol_tmp, FluxX, FluxY, dx, dy, nx, ny, scheme);
    }
    else {
        cout << "Unknown scheme" << endl;
    }

    for (int i = 0; i < ny; i++){
        for (int j = 0; j < nx; j++){
            im1 = id_map(i-1, 0, ny-1, ny);
            jm1 = id_map(j-1, 0, nx-1, nx);
            sol_tmp[i][j] = 0.75*sol[i][j] + 0.25*sol_tmp[i][j] - 0.25*dt/dx*(FluxX[i][j]-FluxX[i][jm1]) - 0.25*dt/dy*(FluxY[i][j]-FluxY[im1][j]);
        }
    }

    if (scheme=="2rdJameson"){
        gradCalc_2rdJameson(sol_tmp, FluxX, FluxY, dx, dy, nx, ny, scheme);
    }
    else if (scheme=="2rdRoe"){
        gradCalc_2rdRoe(sol_tmp, FluxX, FluxY, dx, dy, nx, ny, scheme);
    }
    else {
        cout << "Unknown scheme" << endl;
    }

    for (int i = 0; i < ny; i++){
        for (int j = 0; j < nx; j++){
            im1 = id_map(i-1, 0, ny-1, ny);
            jm1 = id_map(j-1, 0, nx-1, nx);
            sol[i][j] = 1.0/3.0*sol[i][j] + 2.0/3.0*sol_tmp[i][j] - 2.0/3.0*dt/dx*(FluxX[i][j]-FluxX[i][jm1]) - 2.0/3.0*dt/dy*(FluxY[i][j]-FluxY[im1][j]);
        }
    }
}


/* ./HW5_isentropic 10.0 100 0.01 0.0 10.0 0.0 10.0 0.1 0.1 2rdJameson */
/* ./HW5_isentropic 10.0 100 0.01 0.0 10.0 0.0 10.0 0.1 0.1 2rdRoe */
int main(int argc, char *argv[]) {
    double elapsedTime; 
    int stepCounter, reportStep;
    double Xstart, Xend;
    double Ystart, Yend;
    double dx, dy, dt, x_i;
    double px, py;
    double temp_e;
    int nx, ny;
    ofstream outFile;

    const double vs = 5.0;
    const double gamma = 1.4;
    const double x0 = 5.0;
    const double y0 = 5.0;
    
    elapsedTime =  atof(argv[1]);
    reportStep = atoi(argv[2]);
    dt = atof(argv[3]);
    Xstart = atof(argv[4]);
    Xend = atof(argv[5]);
    Ystart = atof(argv[6]);
    Yend = atof(argv[7]);
    dx = atof(argv[8]);
    dy = atof(argv[9]);
    const string SPATIAL_DIS = argv[10]; 

    cout << "dx: " << dx << "\tdy: " << dy << "\tdt: " << dt << endl;
    cout << "ReportTime Interval: " << reportStep*dt << endl;
    cout << "Solution Domain: " << "[" << Xstart << ", " << Xend << "]*[" << Ystart << ", " << Yend << "]" << endl;

    nx = int((Xend - Xstart) / dx);
    ny = int((Yend - Ystart) / dy);

    Matrix<Var> Xi(nx, ny);
    Matrix<Var> FluxXiX(nx, ny);
    Matrix<Var> FluxXiY(nx, ny);

    for (int i = 0; i < ny; i++){
        py = Ystart + (i+0.5) * dy - y0;
        for (int j = 0; j < nx; j++){
            px = Xstart + (j+0.5) * dx - x0;
            temp_e = exp(1.0-px*px-py*py);
            Xi[i][j].v1 = pow(1.0-((gamma-1.0)*vs*vs)/(8.0*gamma*M_PI*M_PI)*temp_e, 1.0/(gamma-1));
            Xi[i][j].v2 = Xi[i][j].v1 * (1.0 - vs/2.0/M_PI*sqrt(temp_e)*py);
            Xi[i][j].v3 = Xi[i][j].v1 * (1.0 + vs/2.0/M_PI*sqrt(temp_e)*px);
            Xi[i][j].v4 = pow(Xi[i][j].v1, gamma) / (gamma - 1.0) 
                        + 0.5*(Xi[i][j].v2*Xi[i][j].v2 + Xi[i][j].v3*Xi[i][j].v3)/Xi[i][j].v1;
        }
    }

    stepCounter = 0;
    for (double t = 0; t <= elapsedTime; t+=dt){
        // Spatial Discretization & Time marching
        step(Xi, FluxXiX, FluxXiY, dx, dy, dt, nx, ny, SPATIAL_DIS);
        stepCounter += 1;
        // Report result
        if (stepCounter % reportStep==0){
            string name = "Case1Res_" + SPATIAL_DIS + "_" + to_string(double(stepCounter*dt)) + ".txt";
            cout << "Time: " << t+dt << ", " << name << " recorded" << endl;
            outFile.open(name, ios::out);
            for (int i = 0; i < ny; i++){
                for (int j = 0; j < nx; j++){
                    outFile << Xi[i][j].v1 << " ";
                }
                outFile << endl;
            }
            outFile.close();
        }
    }

    return 0;
}

