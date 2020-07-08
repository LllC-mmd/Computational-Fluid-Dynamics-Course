#include "global.hpp"
#include "solve.hpp"
#include <armadillo>


void solve_RK3(double dt){
    double dt_RK = dt;
    double dt_left = dt;
    // ---Determine the true time step under CFL constraints
    while (dt_left > 0.0){
        dt_RK = dt;
        for (int i = 0; i < nCell; i++){
            dt_RK = min(dt_RK, globalCell[i].getCellLocalTime(dt_RK));
        }
        cout << dt_left << "\t" << dt_RK << endl;
        dt_RK = max(dt_min, dt_RK);
        dt_RK = min(dt_RK, dt_left);
        RK3_stepping(dt_RK);
        dt_left -= dt_RK;
    }
}

void RK3_stepping(double dt){
    int oId, nId;
    double eta;

    for (int i = 0; i < nCell; i++){
        globalCell[i].update_RK3_Init();
    }
    
    for (int a = 0; a < 3; a++){
        // ---Variable reconstruction
        for (int i = 0; i < nCell; i++){
            /*if (i==2322){
                cout << globalCell[i].cId << endl;
            }*/
            eta = globalCell[i].getValU_alpha(a).eta;
            if (eta < globalCell[i].z_node[0]+epsilon_wd){
                globalCell[i].CellRec_dry(globalCell[i].getValU_alpha(a));
            }
            else if (eta < globalCell[i].z_node[2]){
                globalCell[i].CellRec_partial(globalCell[i].getValU_alpha(a));
            }
            else{
                globalCell[i].CellRec_MUSCL(a);
            }
        }
        for (int i = 0; i < nEdge; i++){
            // ------------Adjust for robustness
            if (globalEdge[i].boundaryFlag==0){
                oId = globalEdge[i].getOCell();
                nId = globalEdge[i].getNCell();
                globalEdge[i].reconstruct_adjust(globalCell[oId-1].getValU_alpha(a), globalCell[nId-1].getValU_alpha(a));
            }
            // ------------Boundary reconstruction
            else{
                globalEdge[i].reconstruct_BC(a, Newton_Raphson_tol);
            }
        }
        // ---Determine the type of cell according to the reconstructed state 
        for (int i = 0; i < nCell; i++){
            globalCell[i].WetDryCheck(globalCell[i].getValU_alpha(a));
        }
        // ---Wet-dry fronts sharpen
        for (int i = 0; i < nCell; i++){
            globalCell[i].update_RK3(dt, a);
        }
    }
    
    for (int i = 0; i < nCell; i++){
        globalCell[i].update_RK3_Step();
    }
}


void solve_SDIRK4(double dt, double dtau, double epsilon, int max_inner_iter){
    int oId, nId;
    double eta;

    double res_h_ini = 0.0;
    double res_h = 0.0;
    double res_tmp;
    
    // For SDIRK4, we need to solve the implicit equation 3 times using Dual Time Stepping
    for (int a = 0; a < 3; a++){
        for (int i = 0; i < nCell; i++){
            globalCell[i].update_SDIRK4_InnerInit();
        }
        // ---Solve the inner problem via LU-SGS where s is the inner time step
        for (int s = 0; s < max_inner_iter; s++){
            // ------Update variational reconstruction coefficients for each cell
            for (int i = 0; i < nCell; i++){
                eta = globalCell[i].getValU_inner().eta;
                if (eta >= globalCell[i].z_node[2]){
                    globalCell[i].VarRecUpdate();
                }
            }
            // ------LU-SGS sweeping
            // ---------call reconstruction procedure for Wet-Dry check and Residual calculation for inner loop
            for (int i = 0; i < nCell; i++){
                eta = globalCell[i].getValU_inner().eta;
                if (eta < globalCell[i].z_node[0]+epsilon_wd){
                    globalCell[i].CellRec_dry(globalCell[i].getValU_inner());
                }
                else if (eta < globalCell[i].z_node[2]){
                    globalCell[i].CellRec_partial(globalCell[i].getValU_inner());
                }
                else{
                    globalCell[i].CellRec_VarRec(globalCell[i].getValU_inner());
                }
            } 
            for (int i = 0; i < nEdge; i++){
                // ------------Adjust for robustness
                if (globalEdge[i].boundaryFlag==0){
                    oId = globalEdge[i].getOCell();
                    nId = globalEdge[i].getNCell();
                    globalEdge[i].reconstruct_adjust(globalCell[oId-1].getValU_inner(), globalCell[nId-1].getValU_inner());
                }
                // ------------Boundary reconstruction
                else{
                    globalEdge[i].reconstruct_BC(a, Newton_Raphson_tol);
                }
            }
            for (int i = 0; i < nCell; i++){
                globalCell[i].WetDryCheck(globalCell[i].getValU_inner());
            }
            // ---------LU-SGS forward sweeping
            for (int i = 0; i < nCell; i++){
                globalCell[i].update_LUSGS_forward(dt, dtau, a);
            }
            // ---------LU-SGS backward sweeping & update
            for (int i = 0; i < nCell; i++){
                res_tmp = fabs(globalCell[i].update_LUSGS_backward(dt, dtau, a));
                if (s == 0){
                    res_h_ini += res_tmp;
                }
                else{
                    res_h += res_tmp;
                }
                globalCell[i].update_SDIRK4_InnerUpdate();
            }
            // ------Stop if the residual is converged
            if (s > 0){
                cout << "Inner iteration " << s << "\tRelative Residual: " << res_h/res_h_ini << "\tAbsolute Residual: " << res_h << endl;
                if (res_h/res_h_ini < epsilon){
                    break;
                }
            }
            if (s==max_inner_iter-1){
                cout << "Maximum Inner Iteration exceeds" << endl;
            }
        }
        // ---Store the middle solution of the inner problem for a-th solution
        for (int i = 0; i < nCell; i++){
            globalCell[i].update_SDIRK4_StoreU(a);
        }
        // ---Store the residual of the middle solution
        // ------Reconstruct the Left and Right side Gaussian points' value using the middle solution
        for (int i = 0; i < nCell; i++){
            globalCell[i].CellRec_VarRec(globalCell[i].getValU_alpha(a));
        }
        // ------Calculate the residual using the reconstructed value
        for (int i = 0; i < nCell; i++){
            globalCell[i].setRes_alpha(a);
        }
    }
    // Update the solution for the outer problem
    for (int i = 0; i < nCell; i++){
        globalCell[i].update_SDIRK4_Step(dt);
        globalCell[i].update_SDIRK4_ClearU();
    }
};
