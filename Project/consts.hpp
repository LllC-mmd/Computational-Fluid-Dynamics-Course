// Constants value for IUHM
#include <map>
#include <array>
#include <vector>
#include <math.h>

using namespace std;

//---------------- Flow2D model part ----------------//
const double Newton_Raphson_tol = 1e-9;

// slope limiter 
const double epsilon_limiter = 1e-12;

// Wetting-Drying processing
const double epsilon_wd = 1e-6;

// h_FLOOD: 1e-2; h_DRY: 5e-4
// h_FLOOD: 5e-3; h_DRY: 1e-4

// global index initialization for mesh 
const double g = 9.81;
const double pi = M_PI;

// Taylor expansion reference
const map<int, vector<vector<int> > > taylor_table = {
    {3, {
        {0, 1}, {1, 0},
        {0, 2}, {1, 1}, {2, 0},
        {0, 3}, {1, 2}, {2, 1}, {3, 0},},
    },
    {4, {{0, 1}, {1, 0},
        {0, 2}, {1, 1}, {2, 0},
        {0, 3}, {1, 2}, {2, 1}, {3, 0},
        {0, 4}, {1, 3}, {2, 2}, {3, 1}, {4, 0},}
    }
};

const map<int, vector<vector<int> > > binomial_table = {
    {1, {
        {1, 1, 0}, {1, 0, 1}}
    },
    {2, {
        {1, 2, 0}, {2, 1, 1}, {1, 0, 2}}
    },
    {3, {
        {1, 3, 0}, {3, 2, 1}, {3, 1, 2}, {1, 0, 3}}
    },
};


//---------------- Gaussian-Quadrature Integral part ----------------//
// Line Integral using G-Q
const map<int, vector<vector<double> > > LineGQ_table = {
    {3, {
        // Weight, alpha: [-1, 1]
        {0.888888888888889, 0.000000000000000},
        {0.555555555555556, 0.774596669241483},
        {0.555555555555556, -0.774596669241483},}
    },
    {4, {
        {0.652145154862546, 0.339981043584856},
        {0.652145154862546, -0.339981043584856},
        {0.347854845137454, 0.861136311594053},
        {0.347854845137454, -0.861136311594053},}
    }
};

// Triangle Integral using G-Q
const map<int, vector<vector<double> > > TriGQ_table = {
    {3, {
        // Weight, zeta_1, zeta_2, zeta_3
        {-0.562500000000000, 0.333333333333333, 0.333333333333333, 0.333333333333333},
        {0.520833333333333, 0.600000000000000, 0.200000000000000, 0.200000000000000},
        {0.520833333333333, 0.200000000000000, 0.600000000000000, 0.200000000000000},
        {0.520833333333333, 0.200000000000000, 0.200000000000000, 0.600000000000000},}
    },
    {4, {
        {0.109951743655322, 0.816847572980459, 0.091576213509771, 0.091576213509771},
        {0.109951743655322, 0.091576213509771, 0.816847572980459, 0.091576213509771},
        {0.109951743655322, 0.091576213509771, 0.091576213509771, 0.816847572980459},
        {0.223381589678011, 0.108103018168070, 0.445948490915965, 0.445948490915965},
        {0.223381589678011, 0.445948490915965, 0.108103018168070, 0.445948490915965},
        {0.223381589678011, 0.445948490915965, 0.445948490915965, 0.108103018168070},}
    }
};

//---------------- Equation solution part ----------------//
// CFL constraints
const double cfl = 0.8;

// Explicit Runge-Kutta: RK3
const double dt_default = 0.1;
const vector<vector<double> > RK3_table = {
    {1.0, 0.0, 1.0},
    {3.0/4.0, 1.0/4.0, 1.0/4.0},
    {1.0/3.0, 2.0/3.0, 2.0/3.0}
};

// Implicit Runge-Kutta: SDIRK4
const double dtau_default = 0.1;
const double zeta = 0.128886400515;
const double omega = 1.3;
const double epsilon_LUSGS_default = 1e-6;
const int maxiter_LUSGS_default = 100;

const array<double, 3> SDIRK4_weight {1.0/(6*pow(2*zeta-1.0, 2)), (4.0*zeta*zeta - 4*zeta + 2.0/3.0)/pow(2*zeta-1.0, 2), 1.0/(6*pow(2*zeta-1.0, 2))};

const map<int, vector<double> > SDIRK4_table = {
    {0, {zeta}},
    {1, {0.5-zeta, zeta}},
    {2, {2*zeta, 1.0-4.0*zeta, zeta}},
};


//---------------- Boundary Condition Reference ----------------//
// 0: inner edge
// 1: boundary edge: free outflow
// 2: boundary edge: no-slip solid wall
// 3: boundary edge with a given h
// 4: boundary edge with a given u
// 5: boundary edge with a given hu
// 6: boundary edge: no-slip solid wall (double face)