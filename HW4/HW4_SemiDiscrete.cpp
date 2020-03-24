#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <string>

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


void gradCalc(double* sol, double* grad, double dx, int num_grid, const string scheme){
    if (scheme=="2rdUpWind"){
        for (int i=0; i<num_grid; i++){
            grad[i] = -(3*sol[id_map(i, 0, num_grid-1, num_grid)] - 4*sol[id_map(i-1, 0, num_grid-1, num_grid)]+
                        sol[id_map(i-2, 0, num_grid-1, num_grid)]) / (2*dx);
        }
    }
    else if (scheme=="2rdCenter"){
        for (int i=0; i<num_grid; i++){
            grad[i] = -(sol[id_map(i+1, 0, num_grid-1, num_grid)] - sol[id_map(i-1, 0, num_grid-1, num_grid)]) / (2*dx);
        }
    }
    else {
        cout << "Unknown scheme" << endl;
    }

}

void step(double* sol, double* grad, double dx, double dt, int num_grid, const string scheme){
    double* sol_tmp = new double[num_grid];

    gradCalc(sol, grad, dx, num_grid, scheme);
    for (int i=0; i<num_grid; i++){
        sol_tmp[i] = sol[i] + dt * grad[i];
    }

    gradCalc(sol_tmp, grad, dx, num_grid, scheme);
    for (int i=0; i<num_grid; i++){
        sol_tmp[i] = 3.0/4.0*sol[i] + 1.0/4.0*(sol_tmp[i]+dt*grad[i]);
    }

    gradCalc(sol_tmp, grad, dx, num_grid, scheme);
    for (int i=0; i<num_grid; i++){
        sol[i] = 1.0/3.0*sol[i] + 2.0/3.0*(sol_tmp[i]+dt*grad[i]);
    }

    delete[] sol_tmp;
}

int main(int argc, char *argv[]) {
    double elapsedTime; 
    int stepCounter, reportStep;
    double start, end;
    double dx, dt, x_i;
    int num_grid;
    double *sol, *grad_X;
    ofstream outFile;
    
    elapsedTime =  atof(argv[1]);
    reportStep = atoi(argv[2]);
    dt = atof(argv[3]);
    start = atof(argv[4]);
    end = atof(argv[5]);
    dx = atof(argv[6]);
    const string SPATIAL_DIS = argv[7]; 
    
    cout << "dx: " << dx << "\t" << "dt: " << dt << endl;
    cout << "ReportTime Interval: " << reportStep*dt << endl;
    cout << "Solution Domain: " << "[" << start << ", " << end << "]" << endl;

    num_grid = int((end-start)/dx) + 1;
    sol = new double[num_grid];
    grad_X = new double[num_grid];

    /* ./HW4_SemiDiscrete 5.0 100 0.01 -0.5 0.5 0.05 2rdUpWind */
    /********Case 1********/
    for (int i=0; i<num_grid; i++){
        x_i = start + i*dx;
        if (x_i < -0.25){
            sol[i] = 0.0;
        }
        else if (x_i <= 0.25){
            sol[i] = 1.0;
        }
        else {
            sol[i] = 0.0;
        }
    }

    stepCounter = 0;
    for (double t = 0; t <= elapsedTime; t+=dt){
        // Spatial Discretization & Time marching
        step(sol, grad_X, dx, dt, num_grid, SPATIAL_DIS);
        stepCounter += 1;
        // Report result
        if (stepCounter % reportStep==0){
            string name = "Case1Res_" + SPATIAL_DIS + "_" + to_string(int(stepCounter*dt)) + ".txt";
            cout << "Time: " << t+dt << ", " << name << " recorded" << endl;
            outFile.open(name, ios::out);
            for (int i = 0; i < num_grid; i++){
                outFile << sol[i] << endl;
            }
            outFile.close();
        }
    }

    /********Case 2********/
    for (int i=0; i<num_grid; i++){
        x_i = start + i*dx;
        sol[i] = sin(4*M_PI*x_i);
    }

    stepCounter = 0;
    for (double t = 0; t <= elapsedTime; t+=dt){
        // Spatial Discretization & Time marching
        step(sol, grad_X, dx, dt, num_grid, SPATIAL_DIS);
        stepCounter += 1;
        // Report result
        if (stepCounter % reportStep==0){
            string name = "Case2Res_" + SPATIAL_DIS + "_" + to_string(int(stepCounter*dt)) + ".txt";
            cout << "Time: " << t+dt << ", " << name << " recorded" << endl;
            outFile.open(name, ios::out);
            for (int i = 0; i < num_grid; i++){
                outFile << sol[i] << endl;
            }
            outFile.close();
        }
    }

    return 0;
}

