#include "integration.hpp"

#include <iostream>
#include <iomanip>


double integral_triGP(double (*f)(double, double, const vector<double>&), 
                        double x1, double y1, double x2, double y2, double x3, double y3, 
                        const vector<double>& argv, int order){
    double val = 0.0;
    double w, zeta_1, zeta_2;
    double x, y, s, t;

    double area = 0.5*fabs((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1));
    for (auto i = TriGQ_table.at(order).begin(); i != TriGQ_table.at(order).end(); i++){
        w = (*i)[0];
        zeta_1 = (*i)[1];
        zeta_2 = (*i)[2];
        s = zeta_2;
        t = zeta_1;
        x = x1 + (x2-x1)*s + (x3-x1)*t;
        y = y1 + (y2-y1)*s + (y3-y1)*t;
        val += w*f(x, y, argv);
    }
    val *= area;
    return val;
}


double integral_lineGP(double (*f)(double, double, const vector<double>&),
                        double x1, double y1, double x2, double y2, const vector<double>& argv, int order){
    double w, alpha;
    double x_GP, y_GP;
    double val = 0.0;
    double xm = 0.5*(x1+x2);
    double ym = 0.5*(y1+y2);

    double length = sqrt(pow(x1-x2, 2) + pow(y1-y2, 2));
    for (auto i = LineGQ_table.at(order).begin(); i != LineGQ_table.at(order).end(); i++){
        w = (*i)[0];
        alpha = (*i)[1];
        x_GP = xm + 0.5*alpha*(x2-x1);
        y_GP = ym + 0.5*alpha*(y2-y1);
        val += 0.5 * w * f(x_GP, y_GP, argv);
    }
    val *= length;
    return val;
}

// ********************************* Zero-basis function *********************************
double delta_xy(double x, double y, const vector<double>& argv){
    // we can check if argv.size() == 6
    double k_x = argv[0];
    double b_x = argv[1];
    double k_y = argv[2];
    double b_y = argv[3];
    double m = argv[4];
    double n = argv[5];
    return pow(k_x*(x-b_x), m)*pow(k_y*(y-b_y), n);
}


double zero_basis_check(double x, double y, const vector<double>& argv){
    double k_x = argv[0];
    double b_x = argv[1];
    double k_y = argv[2];
    double b_y = argv[3];
    double m = argv[4];
    double n = argv[5];
    double mean = argv[6];
    double val = pow(k_x*(x-b_x), m)*pow(k_y*(y-b_y), n) - mean;
    return val;
}


double phi_single(double x, double y, const vector<double>& argv){
    // parameters for phi_m
    double mean_phi = argv[0]; 
    double m_phi = argv[1];
    double n_phi = argv[2];
    // other parameters
    double k_x = argv[3];
    double b_x = argv[4];
    double k_y = argv[5];
    double b_y = argv[6]; 

    double val = pow(k_x*(x-b_x), m_phi)*pow(k_y*(y-b_y), n_phi) - mean_phi;

    return val;
}


double grad_phi_zero_ml(double x, double y, const vector<double>& argv){
    // parameters for phi_m
    double mean_phim = argv[0]; 
    double m_phim = argv[1];
    double n_phim = argv[2];
    double k_x_m = argv[3];
    double b_x_m = argv[4];
    double k_y_m = argv[5];
    double b_y_m = argv[6];
    // parameters for phi_l
    double mean_phil = argv[7]; 
    double m_phil = argv[8];
    double n_phil = argv[9];
    double k_x_l = argv[10];
    double b_x_l = argv[11];
    double k_y_l = argv[12];
    double b_y_l = argv[13]; 

    double val = (pow(k_x_m*(x-b_x_m), m_phim)*pow(k_y_m*(y-b_y_m), n_phim) - mean_phim) 
                    * (pow(k_x_l*(x-b_x_l), m_phil)*pow(k_y_l*(y-b_y_l), n_phil) - mean_phil);

    return val;
}


double grad_phi_ml(double x, double y, const vector<double>& argv){
     // parameters for phi_m
    double m_phim = argv[0];
    double n_phim = argv[1];
    double k_x_m = argv[2];
    double b_x_m = argv[3];
    double k_y_m = argv[4];
    double b_y_m = argv[5];
    // parameters for phi_l
    double m_phil = argv[6];
    double n_phil = argv[7];
    double k_x_l = argv[8];
    double b_x_l = argv[9];
    double k_y_l = argv[10];
    double b_y_l = argv[11]; 
    // other parameters
    double theta = argv[12];
    int order = argv[13];

    double val_m = 0.0;
    double val_l = 0.0;
    double cos_n = cos(theta);
    double sin_n = sin(theta);

    int b_coef, s, r;
    for (auto iter = binomial_table.at(order).begin(); iter != binomial_table.at(order).end(); iter++){
        b_coef = (*iter)[0];
        s = (*iter)[1];
        r = (*iter)[2];

        val_m += b_coef * pow(cos_n, s) * pow(sin_n, r) * dReLU(m_phim-s) * dReLU(n_phim-r) 
                    * ReFac(m_phim, s)*pow(k_x_m, s) * ReFac(n_phim, r)*pow(k_y_m, r)
                    * pow(k_x_m*(x-b_x_m), m_phim-s) * pow(k_y_m*(y-b_y_m), n_phim-r);
        val_l += b_coef * pow(cos_n, s) * pow(sin_n, r) * dReLU(m_phil-s) * dReLU(n_phil-r) 
                    * ReFac(m_phil, s)*pow(k_x_l, s) * ReFac(n_phil, r)*pow(k_y_l, r)
                    * pow(k_x_l*(x-b_x_l), m_phil-s) * pow(k_y_l*(y-b_y_l), n_phil-r);
    }

    return val_m * val_l;
}



/*
// Functions to test the accuracy of integral_triGP and integral_lineGP
double f_testTri(double x, double y, const vector<double>& argv){
    return sin(x)*sin(y);
}

double f_testLine(double x, double y, const vector<double>& argv){
    return x+y;
}

int main(int argc, char *argv[]) {
    double val_tri = integral_triGP(f_testTri, 0.0, 0.0, 0.7, 0.5, 1.0, -1.0, {0.0, 0.0}, 4);
    std::cout << "C++ Integral over Triangle: " << std::setprecision(10) << val_tri << std::endl;
    double val_line = integral_lineGP(f_testLine, 1.0, 1.0, 0.0, 0.0, {0.0, 0.0}, 4);
    std::cout << "C++ Integral along Line: " << val_line << std::endl;
    return 0;
}*/
