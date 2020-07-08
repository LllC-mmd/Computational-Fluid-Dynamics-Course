#include "consts.hpp"
#include "flowMath.hpp"
#include "math.h"
#include <string>


double integral_triGP(double (*f)(double, double, const vector<double>&), 
                        double x1, double y1, double x2, double y2, double x3, double y3, 
                        const vector<double>& argv, int order);

double integral_lineGP(double (*f)(double, double, const vector<double>&),
                        double x1, double y1, double x2, double y2, const vector<double>& argv, int order);


double delta_xy(double x, double y, const vector<double>& argv);
double zero_basis_check(double x, double y, const vector<double>& argv);
double phi_single(double x, double y, const vector<double>& argv);
double grad_phi_zero_ml(double x, double y, const vector<double>& argv);
double grad_phi_ml(double x, double y, const vector<double>& argv);