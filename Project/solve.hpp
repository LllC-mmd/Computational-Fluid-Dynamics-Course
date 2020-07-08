#include "mesh.hpp"

void solve_RK3(double dt);
void RK3_stepping(double dt);

void solve_SDIRK4(double dt, double dtau, double epsilon, int max_inner_iter);
