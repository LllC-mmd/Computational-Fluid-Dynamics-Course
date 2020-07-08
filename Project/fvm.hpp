#include "mesh.hpp"

U rotate(U u, double theta);
U rotate_back(U u, double theta);
FluxF rotate_back(FluxF f, double theta);

tuple<U, U> grad_Calc(const U& u1, double x1, double y1, const U& u2, double x2, double y2, const U& u3, double x3, double y3);
tuple<double, double> RoeAvg(const U& A, const U& B);