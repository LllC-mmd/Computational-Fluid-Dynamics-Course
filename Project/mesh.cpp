#include "mesh.hpp"
#include "global.hpp"
#include <math.h>
#include <set>
#include <iostream>
#include <algorithm>
#include <iomanip>

using namespace std;

Node::Node(int i, double x, double y, double ele){
    nId = i;
    px = x;
    py = y;
    zb = ele;
};


Edge::Edge(int eIdx, int fId, int tId, int flag, char* source, U EdgeU, FluxF EdgeF){
    double theta_t;
    double xf, yf, xt, yt; 

	eId = eIdx;
	nodeId_Edge[0] = fId;
	nodeId_Edge[1] = tId;
	cellId_Edge[0] = -1;
	cellId_Edge[1] = -1;
    // boundaryFlag = 0: inner edge, i.e., no BC
    // boundaryFlag = 1: solid wall boundary
    // boundaryFlag = 2: open boundary where h is given
    // boundaryFlag = 3: open boundary where uN is given
    // boundaryFlag = 4: open boundary where h*uN is given
    boundaryFlag = flag; 
    BCsource = source;

	int nGP;
	// For MUSCL Reconstruction, we need 1 Gaussian Points
	if (recMethod == "MUSCL"){
		nGP = 1;
	}
	// For Varitional Reconstruction, we need 4 Gaussian Points
	else if (recMethod == "VariationalRec"){
		nGP = 3;
	}
    else{
        cout << "Not Implemented Reconstruction method!" << endl;
    }
    // Boundary edge with double face
    if (flag > 5){
        nGP *= 2;
    }
	UL_GP.resize(nGP);
	UR_GP.resize(nGP);
	F_GP.resize(nGP);
	fill(UL_GP.begin(), UL_GP.end(), EdgeU);
	fill(UR_GP.begin(), UR_GP.end(), EdgeU);
	fill(F_GP.begin(), F_GP.end(), EdgeF);

    xf = globalNode[fId-1].px;
    yf = globalNode[fId-1].py;
    xt = globalNode[tId-1].px;
    yt = globalNode[tId-1].py;

    cx = 0.5 * (xf + xt);
    cy = 0.5 * (yf + yt);
    length = sqrt(pow(xf-xt, 2)+pow(yf-yt, 2));
    z_Edge[0] = min(globalNode[fId-1].zb, globalNode[tId-1].zb);
    z_Edge[1] = max(globalNode[fId-1].zb, globalNode[tId-1].zb);
    zm = 0.5 * (z_Edge[0] + z_Edge[1]);
    if (xf == xt){
        theta_t = pi/2;
    }
    else{
        theta_t = atan((yf - yt)/(xf - xt));
    }
    theta_n = pi/2 + theta_t;
};


// constructor for Cell
Cell::Cell(int cIdx, int e1, int e2, int e3, double n_r, double rc, U u){
    int N1, N2, N3;
    int n1, n2, n3;
    double x1, y1, z1, x2, y2, z2, x3, y3, z3;

    cId = cIdx;
    n = n_r;
    runoff_coef = rc;
    _Uc = u;
    vfv_order = 3;
    if (steppingMethod == "RK-3"){
        _U_alpha.resize(4);
    }
    else{
        _U_alpha.resize(3);
    }
    for (int i = 0; i < _U_alpha.size(); i++){
        _U_alpha[i].h = 0.0;
        _U_alpha[i].hu = 0.0;
        _U_alpha[i].hv = 0.0;
    }
    if (recMethod == "VariationalRec"){
        _Res_alpha.resize(2);
        for (int i = 0; i < _Res_alpha.size(); i++){
            _Res_alpha[i].v1 = 0.0;
            _Res_alpha[i].v2 = 0.0;
            _Res_alpha[i].v3 = 0.0;
        }
    }
    
    // Set the nodes of the cell
    N1 = globalEdge[e1-1].getFrom();
    N2 = globalEdge[e1-1].getTo();
    if ((globalEdge[e2-1].getFrom() == N1) || (globalEdge[e2-1].getFrom() == N2)){
        N3 = globalEdge[e2-1].getTo();
    }
    else{
        N3 = globalEdge[e2-1].getFrom();
    }

    x1 = globalNode[N1-1].px;
    y1 = globalNode[N1-1].py;
    x2 = globalNode[N2-1].px;
    y2 = globalNode[N2-1].py;
    x3 = globalNode[N3-1].px;
    y3 = globalNode[N3-1].py;

    if ((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1) > 0){
        nodeId_Cell[0] = N1;
        nodeId_Cell[1] = N2;
        nodeId_Cell[2] = N3;
    }
    else{
        nodeId_Cell[0] = N1;
        nodeId_Cell[1] = N3;
        nodeId_Cell[2] = N2;
    }

    // Set the edges of the cell
    if (((globalEdge[e1-1].getFrom() == nodeId_Cell[0]) && (globalEdge[e1-1].getTo() == nodeId_Cell[1])) || ((globalEdge[e1-1].getFrom() == nodeId_Cell[1]) && (globalEdge[e1-1].getTo() == nodeId_Cell[0]))){
        edgeId_Cell[2] = e1;
        edgeBCTag_Cell[2] = globalEdge[e1-1].boundaryFlag;
    }
    else if (((globalEdge[e1-1].getFrom() == nodeId_Cell[1]) && (globalEdge[e1-1].getTo() == nodeId_Cell[2])) || ((globalEdge[e1-1].getFrom() == nodeId_Cell[2]) && (globalEdge[e1-1].getTo() == nodeId_Cell[1]))){
        edgeId_Cell[0] = e1;
        edgeBCTag_Cell[0] = globalEdge[e1-1].boundaryFlag;
    }
    else{
        edgeId_Cell[1] = e1;
        edgeBCTag_Cell[1] = globalEdge[e1-1].boundaryFlag;
    }

    if (((globalEdge[e2-1].getFrom() == nodeId_Cell[0]) && (globalEdge[e2-1].getTo() == nodeId_Cell[1])) || ((globalEdge[e2-1].getFrom() == nodeId_Cell[1]) && (globalEdge[e2-1].getTo() == nodeId_Cell[0]))){
        edgeId_Cell[2] = e2;
        edgeBCTag_Cell[2] = globalEdge[e2-1].boundaryFlag;
    }
    else if (((globalEdge[e2-1].getFrom() == nodeId_Cell[1]) && (globalEdge[e2-1].getTo() == nodeId_Cell[2])) || ((globalEdge[e2-1].getFrom() == nodeId_Cell[2]) && (globalEdge[e2-1].getTo() == nodeId_Cell[1]))){
        edgeId_Cell[0] = e2;
        edgeBCTag_Cell[0] = globalEdge[e2-1].boundaryFlag;
    }
    else{
        edgeId_Cell[1] = e2;
        edgeBCTag_Cell[1] = globalEdge[e2-1].boundaryFlag;
    }

    if (((globalEdge[e3-1].getFrom() == nodeId_Cell[0]) && (globalEdge[e3-1].getTo() == nodeId_Cell[1])) || ((globalEdge[e3-1].getFrom() == nodeId_Cell[1]) && (globalEdge[e3-1].getTo() == nodeId_Cell[0]))){
        edgeId_Cell[2] = e3;
        edgeBCTag_Cell[2] = globalEdge[e3-1].boundaryFlag;
    }
    else if (((globalEdge[e3-1].getFrom() == nodeId_Cell[1]) && (globalEdge[e3-1].getTo() == nodeId_Cell[2])) || ((globalEdge[e3-1].getFrom() == nodeId_Cell[2]) && (globalEdge[e3-1].getTo() == nodeId_Cell[1]))){
        edgeId_Cell[0] = e3;
        edgeBCTag_Cell[0] = globalEdge[e3-1].boundaryFlag;
    }
    else{
        edgeId_Cell[1] = e3;
        edgeBCTag_Cell[1] = globalEdge[e3-1].boundaryFlag;
    }

    // Initialize the neighbour cell of the cell
    neighborCellId_Cell.fill(-255);

    n1 = nodeId_Cell[0];
    n2 = nodeId_Cell[1];
    n3 = nodeId_Cell[2];

    x1 = globalNode[n1-1].px;
    y1 = globalNode[n1-1].py;
    z1 = globalNode[n1-1].zb;
    x2 = globalNode[n2-1].px;
    y2 = globalNode[n2-1].py;
    z2 = globalNode[n2-1].zb;
    x3 = globalNode[n3-1].px;
    y3 = globalNode[n3-1].py;
    z3 = globalNode[n3-1].zb;

    double x_max = max(x1, max(x2, x3));
    double x_min = min(x1, min(x2, x3));
    double y_max = max(y1, max(y2, y3));
    double y_min = min(y1, min(y2, y3));
    dx_max = 0.5*(x_max-x_min);
    dy_max = 0.5*(y_max-y_min);

    // Set the center elevation and coordinates
    cpx = (x1+x2+x3)/3;
    cpy = (y1+y2+y3)/3;
    cpz = (z1+z2+z3)/3;

    z_node[0] = z1;
    z_node[1] = z2;
    z_node[2] = z3;
    sort(z_node.begin(), z_node.end());

    _Uc = VFR_set_eta(_Uc);
    // Set the area and slope
    area = fabs((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1))/2;
    ix = 1.0/(2.0*area)*((y2-y3)*z1 + (y3-y1)*z2 + (y1-y2)*z3);
    iy = 1.0/(2.0*area)*((x3-x2)*z1 + (x1-x3)*z2 + (x2-x1)*z3);

    globalNode[n1-1].CellId_Node.push_back(cIdx);
    globalNode[n2-1].CellId_Node.push_back(cIdx);
    globalNode[n3-1].CellId_Node.push_back(cIdx);
}


// Updata free surface elevation according to new h
U Cell::VFR_set_eta(const U& u){
    U u_new;
    u_new.h = u.h;
    u_new.hu = u.hu;
    u_new.hv = u.hv;
    double z1 = z_node[0];
    double z2 = z_node[1];
    double z3 = z_node[2];
    double r1 = z3 - 3.0*z1;
    double r2 = 3.0*u_new.h*z1 - 3.0*u_new.h*z3 - z3*z2 + z1*z2 + z1*z1;
    
    double eta_1 = z1 + pow(3.0*u_new.h*(z2-z1)*(z3-z1), 1.0/3.0);
    double eta_2 = 0.5*(-r1 + sqrt(r1*r1 - 4.0*r2));
    double eta_3 = u_new.h + (z1+z2+z3)/3.0;

    if (u_new.h <= machine_epsilon){
        u_new.eta = z1;
    }
    else if ((z1 < eta_1) && (eta_1 <= z2)){
        u_new.eta = eta_1;
    }
    else if ((z2 < eta_2) && (eta_2 <= z3)){
        u_new.eta = eta_2;
    }
    else{
        u_new.eta = eta_3;
    }
    return u_new;
}


double Cell::VFR_set_h(double eta){
    double z1 = z_node[0];
    double z2 = z_node[1];
    double z3 = z_node[2];

    if (eta > z3){
        return eta - (z1+z2+z3)/3.0;
    }
    else if (eta > z2){
        return (eta*eta + eta*z3 - 3.0*eta*z1 - z3*z2 + z1*z2 + z1*z1) / 3.0 / (z3-z1);
    }
    else if (eta > z1){
        return pow(eta-z1, 3.0) / 3.0 / (z2-z1) / (z3-z1);
    }
    else{
        return 0.0;
    }
}

// Set Variational Reconstruction coefficients
void Cell::VarRecInit(int Nc){
    // mean value for basis functions over the element
    basis_mean.resize(Nc);
    // A for h, hu, hv
	A_eta.resize(Nc, Nc);
    A_eta.fill(0.0);
    A_hu.resize(Nc, Nc);
    A_hu.fill(0.0);
    A_hv.resize(Nc, Nc);
    A_hv.fill(0.0);
    // B for h, hu, hv with each neighbour cell
    B_eta.resize(getNeighborCellCount());
    for (int i = 0; i < B_eta.size(); i++){
        B_eta[i].resize(Nc, Nc);
        B_eta[i].fill(0.0);
    }
    B_hu.resize(getNeighborCellCount());
    for (int i = 0; i < B_hu.size(); i++){
        B_hu[i].resize(Nc, Nc);
        B_hu[i].fill(0.0);
    }
    B_hv.resize(getNeighborCellCount());
    for (int i = 0; i < B_hv.size(); i++){
        B_hv[i].resize(Nc, Nc);
        B_hv[i].fill(0.0);
    }
    // reconstruction coefficients for h, hu, hv
    coef_eta.resize(Nc);
    coef_eta.fill(0.0);
    coef_hu.resize(Nc);
    coef_hu.fill(0.0);
    coef_hv.resize(Nc);
    coef_hv.fill(0.0);
    // b for h, hu, hv
    b_eta_unit.resize(Nc, 3);
    b_hu_unit.resize(Nc, 3);
    b_hv_unit.resize(Nc, 3);
}

// calculate the mean value of basis functions used in the Variational Reconstruction
void Cell::VarRecPrepare_basisMean(){
    int m_m, n_m, m_l, n_l;
    int xf, yf, xt, yt;
    // 4th order by default
    int k = vfv_order;
    int Nc = (k+1)*(k+2)/2-1;
    double x1 = globalNode[nodeId_Cell[0]-1].px;
    double y1 = globalNode[nodeId_Cell[0]-1].py;
    double x2 = globalNode[nodeId_Cell[1]-1].px;
    double y2 = globalNode[nodeId_Cell[1]-1].py;
    double x3 = globalNode[nodeId_Cell[2]-1].px;
    double y3 = globalNode[nodeId_Cell[2]-1].py;
    // Initialization
    VarRecInit(Nc);
    // Determine the zero-basis functions
    int c = 0;
    for (int p = 1; p <= k; p++){
        for (int m = 0; m <= p; m++){
            basis_mean(c) = 1.0/area * integral_triGP(delta_xy, x1, y1, x2, y2, x3, y3, {1.0/dx_max, cpx, 1.0/dy_max, cpy, double(m), double(p-m)}, 3);
            c += 1;
        }
    }
}

// calculate A, B, b (for unit u) used in the Variational Reconstruction
void Cell::VarRecPrepare_coefMat(){
    int eId, ncId;
    double d_ij, theta;
    int m_m, n_m, m_l, n_l;
    int xf, yf, xt, yt;
    double basis_mean_j, dx_max_j, dy_max_j, cpx_j, cpy_j; 
    double val_tmp_A = 0.0;
    double val_tmp_B = 0.0;
    double val_tmp_b = 0.0;
    // 4th order by default
    int k = vfv_order;
    int Nc = (k+1)*(k+2)/2-1;

    // Determine the matrix A using line integral
    for (int l = 0; l < Nc; l++){
        m_l = taylor_table.at(k)[l][0];  // order of (x-x_i)/dx in phi_m
        n_l = taylor_table.at(k)[l][1];  // order of (y-y_i)/dy in phi_m
        for (int m = 0; m < Nc; m++){
            m_m = taylor_table.at(k)[m][0];
            n_m = taylor_table.at(k)[m][1];
            // ---sum over edges
            int num_nc = 0;
            for (int j = 0; j < 3; j++){
                eId = edgeId_Cell[j];
                theta = globalEdge[eId-1].theta_n;
                xf = globalNode[globalEdge[eId-1].getFrom()-1].px;
                yf = globalNode[globalEdge[eId-1].getFrom()-1].py;
                xt = globalNode[globalEdge[eId-1].getTo()-1].px;
                yt = globalNode[globalEdge[eId-1].getTo()-1].py;
                // ------Inner edge
                if (edgeBCTag_Cell[j] == 0){
                    ncId = neighborCellId_Cell[j];
                    d_ij = sqrt(pow(cpx-globalCell[ncId-1].cpx, 2)+ pow(cpy-globalCell[ncId-1].cpy, 2));
                    dx_max_j = globalCell[ncId-1].dx_max;
                    dy_max_j = globalCell[ncId-1].dy_max;
                    cpx_j = globalCell[ncId-1].cpx;
                    cpy_j = globalCell[ncId-1].cpy;
                    // sum over all of the derivatives whose order <= k
                    for (int p = 0; p <= k; p++){
                        if (p==0){
                            basis_mean_j = globalCell[ncId-1].basis_mean(l);
                            val_tmp_A = 1.0/d_ij * integral_lineGP(grad_phi_zero_ml, xf, yf, xt, yt, {basis_mean(m), double(m_m), double(n_m), 1.0/dx_max, cpx, 1.0/dy_max, cpy,
                                                                                                        basis_mean(l), double(m_l), double(n_l), 1.0/dx_max, cpx, 1.0/dy_max, cpy}, 3);
                            val_tmp_B = 1.0/d_ij * integral_lineGP(grad_phi_zero_ml, xf, yf, xt, yt, {basis_mean(m), double(m_m), double(n_m), 1.0/dx_max, cpx, 1.0/dy_max, cpy, 
                                                                                                        basis_mean_j, double(m_l), double(n_l), 1.0/dx_max_j, cpx_j, 1.0/dy_max_j, cpy_j}, 3);                                               
                        }
                        else{
                            val_tmp_A = pow(d_ij, 2*p-1)/pow(ReFac(p, p), 2) * integral_lineGP(grad_phi_ml, xf, yf, xt, yt, {double(m_m), double(n_m), 1.0/dx_max, cpx, 1.0/dy_max, cpy, 
                                                                                                                                double(m_l), double(n_l), 1.0/dx_max, cpx, 1.0/dy_max, cpy, theta, double(p)}, 3);
                            val_tmp_B = pow(d_ij, 2*p-1)/pow(ReFac(p, p), 2) * integral_lineGP(grad_phi_ml, xf, yf, xt, yt, {double(m_m), double(n_m), 1.0/dx_max, cpx, 1.0/dy_max, cpy, 
                                                                                                                                double(m_l), double(n_l), 1.0/dx_max_j, cpx_j, 1.0/dy_max_j, cpy_j, theta, double(p)}, 3);
                        }
                        A_eta(m, l) += val_tmp_A;
                        A_hu(m, l) += val_tmp_A;
                        A_hv(m, l) += val_tmp_A;
                        B_eta[num_nc](m, l) += val_tmp_B;
                        B_hu[num_nc](m, l) += val_tmp_B;
                        B_hv[num_nc](m, l) += val_tmp_B;
                        // unit b
                        val_tmp_b = 1.0/d_ij * integral_lineGP(phi_single, xf, yf, xt, yt, {basis_mean(l), double(m_l), double(n_l), 1.0/dx_max, cpx, 1.0/dy_max, cpy}, 3);
                        b_eta_unit(l, j) += val_tmp_b;
                        b_hu_unit(l, j) += val_tmp_b;
                        b_hv_unit(l, j) += val_tmp_b;
                    }
                    num_nc += 1;
                }
                // ------Boundary edge with no-slip solid BC
                else if (edgeBCTag_Cell[j] == 2){
                    d_ij = sqrt(pow(cpx-0.5*(xf+xt), 2)+ pow(cpy-0.5*(yf+yt), 2));
                    // A
                    val_tmp_A = 4.0/d_ij * integral_lineGP(grad_phi_zero_ml, xf, yf, xt, yt, {basis_mean(m), double(m_m), double(n_m), 1.0/dx_max, cpx, 1.0/dy_max, cpy, 
                                                                                                basis_mean(l), double(m_l), double(n_l), 1.0/dx_max, cpx, 1.0/dy_max, cpy}, 3);
                    A_hu(m, l) += val_tmp_A;
                    A_hv(m, l) += val_tmp_A;
                    // unit b
                    val_tmp_b = -4.0/d_ij * integral_lineGP(phi_single, xf, yf, xt, yt, {basis_mean(l), double(m_l), double(n_l), 1.0/dx_max, cpx, 1.0/dy_max, cpy}, 3);
                    b_hu_unit(l, j) += val_tmp_b;
                    b_hv_unit(l, j) += val_tmp_b;
                }
                // ------Boundary edge with free outflow BC, IJF==0
            }
        }
    }
}


// Given Cell-Edge connection, construct topological relationship for Edge
void topSearch(){
    int o, n, nId;
    float cmx, cmy;
    int eId;

    // determine the topological relationship of the edge
    for (int i=0; i<nCell; i++){
        // set the owner cell, neightbour cell and the opposite points of the edge
        for (int j = 0; j < globalCell[i].edgeId_Cell.size(); j++){
            eId = globalCell[i].edgeId_Cell[j];
            cmx = globalCell[i].cpx - globalEdge[eId-1].cx;
            cmy = globalCell[i].cpy - globalEdge[eId-1].cy;
            // If the owner cell of this edge has not been set
            // then set this cell as its owner cell and also set its left opposite point
            if (globalEdge[eId-1].getOCell() == -1){
                globalEdge[eId-1].setOCell(globalCell[i].cId);
                // Make sure the normal vector points from the owner to the neighbour 
                if (cos(globalEdge[eId-1].theta_n)*cmx+sin(globalEdge[eId-1].theta_n)*cmy > 0){
                    globalEdge[eId-1].theta_n += pi;
                }
            }
            // else the owner cell of this edge has been set
            // then set this cell as its neighbour cell and also set its right opposite point
            else{
                globalEdge[eId-1].setNCell(globalCell[i].cId);
            }
        }
    }
    // Determine the neighbour cell of a cell
    for (int i=0; i<nEdge; i++){
        o = globalEdge[i].getOCell();
        n = globalEdge[i].getNCell();
        // For boundary egdes, n = -1
        if (n >= 0){
            globalCell[n-1].setNeighbourCell(o, globalEdge[i].eId);
            globalCell[o-1].setNeighbourCell(n, globalEdge[i].eId);
        }
    }
}
