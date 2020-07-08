#include "global.hpp"
#include "preproc.hpp"
#include "solve.hpp"
#include "postprec.hpp"

#include <stdlib.h>
#include <string>
#include <fstream>
#include <iostream>
#include <regex>
#include <armadillo>

using namespace std;

int nEdge;
int nCell;
int nNode;

vector<Node> globalNode;
vector<Edge> globalEdge;
vector<Cell> globalCell;

double dt_min;
string recMethod;
string fluxMethod;
string steppingMethod;
double h_FLOOD;
double h_DRY;
double u_max;

int main(int argc, char *argv[]) {
    // Main part of IUHM

    //---------- Initialize Simulation Settings ----------
    string NodeGeometryFile;
    string EdgeGeometryFile;
    string CellGeometryFile;
    
    double dt_sim;
    double dtau;
    double dt_report;
    double elapsedTime;
    double epsilon_LUSGS;
    double maxiter_LUSGS;
    string rainfallFile;
    double dt_rainfall;
    string resDir;
    ifstream prm_file;
    string line;
    smatch m;
    regex r_nNode("Number of Node:\\s+(.*)");
    regex r_NodeGeometryFile("Node Geometry File:\\s+(.*)");
    regex r_nEdge("Number of Edge:\\s+(.*)");
    regex r_EdgeGeometryFile("Edge Geometry File:\\s+(.*)");
    regex r_nCell("Number of Cell:\\s+(.*)");
    regex r_CellGeometryFile("Cell Geometry File:\\s+(.*)");
    regex r_recMethod("Reconstruction Method:\\s+(.*)");
    regex r_fluxMethod("Flux Calculator:\\s+(.*)");
    regex r_hFlood("h_FLOOD:\\s+(.*)");
    regex r_hDry("h_DRY:\\s+(.*)");
    regex r_uMax("U_max:\\s+(.*)");
    regex r_steppingMethod("Time Stepping Method:\\s+(.*)");
    regex r_dt_sim("Maximum TimeStep:\\s+(.*)");
    regex r_dt_min("Minimum TimeStep:\\s+(.*)");
    regex r_dtau("Inner TimeStep:\\s+(.*)");
    regex r_dt_report("Report TimeStep:\\s+(.*)");
    regex r_elapsedTime("Elapsed Time:\\s+(.*)");
    regex r_epsilon_LUSGS("LU-SGS residual:\\s+(.*)");
    regex r_maxiter_LUSGS("LU-SGS max iteration:\\s+(.*)");
    regex r_rainfallFile("Rainfall Record File:\\s+(.*)");
    regex r_dt_rainfall("Rainfall TimeStep:\\s+(.*)");
    regex r_result("Results saving Directory:\\s+(.*)");

    prm_file.open(argv[1], ios::out);
    while(getline(prm_file, line)){
        if (line[0] != '#'){
            if (regex_match(line, m, r_nNode)){
                nNode = atoi(m.str(m.size()-1).c_str());
            }
            if (regex_match(line, m, r_NodeGeometryFile)){
                NodeGeometryFile = m.str(m.size()-1);
            }
            if (regex_match(line, m, r_nEdge)){
                nEdge = atoi(m.str(m.size()-1).c_str());
            }
            if (regex_match(line, m, r_EdgeGeometryFile)){
                EdgeGeometryFile = m.str(m.size()-1);
            }
            if (regex_match(line, m, r_nCell)){
                nCell = atoi(m.str(m.size()-1).c_str());
            }
            if (regex_match(line, m, r_CellGeometryFile)){
                CellGeometryFile = m.str(m.size()-1);
            }
            if (regex_match(line, m, r_recMethod)){
                recMethod = m.str(m.size()-1);
            }
            if (regex_match(line, m, r_fluxMethod)){
                fluxMethod = m.str(m.size()-1);
            }
            if (regex_match(line, m, r_hFlood)){
                h_FLOOD = atof(m.str(m.size()-1).c_str());
            }
            if (regex_match(line, m, r_hDry)){
                h_DRY = atof(m.str(m.size()-1).c_str());
            }
            if (regex_match(line, m, r_uMax)){
                u_max = atof(m.str(m.size()-1).c_str());
            }
            if (regex_match(line, m, r_steppingMethod)){
                steppingMethod = m.str(m.size()-1);
            }
            if (regex_match(line, m, r_dt_sim)){
                dt_sim = atof(m.str(m.size()-1).c_str());
            }
            if (regex_match(line, m, r_dt_min)){
                dt_min = atof(m.str(m.size()-1).c_str());
            }
            if (regex_match(line, m, r_dtau)){
                dtau = atof(m.str(m.size()-1).c_str());
            }
            if (regex_match(line, m, r_dt_report)){
                dt_report = atof(m.str(m.size()-1).c_str());
            }
            if (regex_match(line, m, r_elapsedTime)){
                elapsedTime = atof(m.str(m.size()-1).c_str());
            }
            if (regex_match(line, m, r_epsilon_LUSGS)){
                epsilon_LUSGS = atof(m.str(m.size()-1).c_str());
            }
            if (regex_match(line, m, r_maxiter_LUSGS)){
                maxiter_LUSGS = atoi(m.str(m.size()-1).c_str());
            }
            if (regex_match(line, m, r_rainfallFile)){
                rainfallFile = m.str(m.size()-1);
            }
            if (regex_match(line, m, r_dt_rainfall)){
                dt_rainfall = atof(m.str(m.size()-1).c_str());
            }
            if (regex_match(line, m, r_result)){
                resDir = m.str(m.size()-1);
            }
        }
    }

    cout << "********************* Input Files *********************" << endl;
    cout << "Node Geometry File: " << NodeGeometryFile << "\t" << "Number of Node: " << nNode << endl;
    cout << "Edge Geometry File: " << EdgeGeometryFile << "\t" << "Number of Edge: " << nEdge << endl;
    cout << "Cell Geometry File: " << CellGeometryFile << "\t" << "Number of Cell: " << nCell << endl;
    cout << "Rainfall recorde file: " << rainfallFile << endl;
    cout << "Rainfall TimeStep: " << dt_rainfall << endl;
    cout << "********************* Simulation Settings *********************" << endl;
    cout << "Reconstruction Method: " << recMethod << endl;
    cout << "Flux Calculator: " << fluxMethod << endl;
    cout << "h_Flood: " << h_FLOOD << "\t" << "h_Dry: " << h_DRY << "\t" << "u_max: " << u_max << endl;
    cout << "Time Stepping Method: " << steppingMethod << endl;
    cout << "Maximum TimeStep: " << dt_sim << endl;
    cout << "Minimum TimeStep: " << dt_min << endl;
    cout << "Inner TimeStep: " << dtau << endl;
    cout << "Report TimeStep: " << dt_report << endl;
    cout << "Elapsed Time: " << elapsedTime << endl;
    cout << "LU-SGS residual: " << epsilon_LUSGS << endl;
    cout << "LU-SGS max iteration: " << maxiter_LUSGS << endl;
    cout << "********************* Output Settings *********************" << endl;
    cout << "Simulation Results would be saved under: " << resDir << endl;
    cout << "***********************************************************" << endl;
    cout << endl;
    
    //---------- Initialize 2D flow objects ----------
    int *NodeID;
    double *NodeX, *NodeY, *NodeEle;
    int *EdgeID, *EdgeF, *EdgeT, *BDFlag;
    char **BCsource;
    int *CellID, *CellE1, *CellE2, *CellE3;
    double *Celln, *Cellrc, *Cell_hIni, *Cell_uIni, *Cell_vIni;
    
    //-----Node part
    NodeID = new int[nNode];
    NodeX = new double[nNode];
    NodeY = new double[nNode];
    NodeEle = new double[nNode];
    hid_t NodeH5file = H5Fopen(NodeGeometryFile.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

    readHDF5(NodeH5file, "/NodeID", NodeID, "int");
    readHDF5(NodeH5file, "/NodeX", NodeX, "double");
    readHDF5(NodeH5file, "/NodeY", NodeY, "double");
    readHDF5(NodeH5file, "/NodeEle", NodeEle, "double");

    for (int i = 0; i < nNode; i++){
        globalNode.push_back(Node(NodeID[i], NodeX[i], NodeY[i], NodeEle[i]));
    }

    H5Fclose(NodeH5file);
    delete[] NodeID;
    delete[] NodeX;
    delete[] NodeY;
    delete[] NodeEle;
    cout << nNode << " nodes are created successfully" << endl;

    //-----Edge part
    EdgeID = new int[nEdge];
    EdgeF = new int[nEdge];
    EdgeT = new int[nEdge];
    BDFlag = new int[nEdge];
    BCsource = new char*[nEdge];
    hid_t EdgeH5file = H5Fopen(EdgeGeometryFile.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    H5File EdgeH5Str(EdgeGeometryFile.c_str(), H5F_ACC_RDONLY);  

    readHDF5(EdgeH5file, "/EdgeID", EdgeID, "int");
    readHDF5(EdgeH5file, "/EdgeF", EdgeF, "int");
    readHDF5(EdgeH5file, "/EdgeT", EdgeT, "int");
    readHDF5(EdgeH5file, "/BDFlag", BDFlag, "int");
    readHDF5_str(EdgeH5Str, "/BCsource", BCsource);

    for (int i = 0; i < nEdge; i++){
        globalEdge.push_back(Edge(EdgeID[i], EdgeF[i], EdgeT[i], BDFlag[i], BCsource[i], {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}));
    }

    H5Fclose(EdgeH5file);
    EdgeH5Str.close();
    delete[] EdgeID;
    delete[] EdgeF;
    delete[] EdgeT;
    delete[] BDFlag;
    delete[] BCsource;
    cout << nEdge << " edges are created successfully" << endl;

    //-----Cell part
    CellID = new int[nCell];
    CellE1 = new int[nCell];
    CellE2 = new int[nCell];
    CellE3 = new int[nCell];
    Celln = new double[nCell];
    Cellrc = new double[nCell];
    Cell_hIni = new double[nCell];
    Cell_uIni = new double[nCell];
    Cell_vIni = new double[nCell];
    hid_t CellH5file = H5Fopen(CellGeometryFile.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

    readHDF5(CellH5file, "/CellID", CellID, "int");
    readHDF5(CellH5file, "/CellE1", CellE1, "int");
    readHDF5(CellH5file, "/CellE2", CellE2, "int");
    readHDF5(CellH5file, "/CellE3", CellE3, "int");
    readHDF5(CellH5file, "/Celln", Celln, "double");
    readHDF5(CellH5file, "/Cellrc", Cellrc, "double");
    readHDF5(CellH5file, "/h", Cell_hIni, "double");
    readHDF5(CellH5file, "/U", Cell_uIni, "double");
    readHDF5(CellH5file, "/V", Cell_vIni, "double");

    for (int i = 0; i < nCell; i++){
        globalCell.push_back(Cell(CellID[i], CellE1[i], CellE2[i], CellE3[i], Celln[i], Cellrc[i], {Cell_hIni[i], Cell_hIni[i]*Cell_uIni[i], Cell_hIni[i]*Cell_vIni[i], 0.0}));
    }
    
    H5Fclose(CellH5file);
    delete[] CellID;
    delete[] CellE1;
    delete[] CellE2;
    delete[] CellE3;
    delete[] Celln;
    delete[] Cellrc;
    delete[] Cell_hIni;
    delete[] Cell_uIni;
    delete[] Cell_vIni;
    cout << nCell << " cells are created successfully\n" << endl;

    //---------- Topological search for connectivity ----------
    topSearch();
    cout << "Mesh topology are constructed successfully\n" << endl;

    /*nodeInfoSave(resDir);
    cout << "Node Info is saved\n" << endl;
    edgeInfoSave(resDir);
    cout << "Edge Info is saved\n" << endl;
    cellInfoSave(resDir);
    cout << "Cell Info is saved\n" << endl;*/

    //---------- Initialize coefficients for Variational Reconstruction ----------
    if (recMethod == "VariationalRec"){
        for (int i = 0; i < nCell; i++){
            globalCell[i].VarRecPrepare_basisMean();
        }
        for (int i = 0; i < nCell; i++){
            globalCell[i].VarRecPrepare_coefMat();
        }
        cout << "Coefficients of varitional reconstruction are initialized" << endl;
    }
    
    //---------- Simulation stepping ----------
    double *rain;
    hid_t RainH5file = H5Fopen(rainfallFile.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    rain = new double[nCell];

    double reportLeft = 0.0;
    double rainfallLeft = 0.0;
    int rain_counter = 0;
    for (double t = 0; t < elapsedTime; t+=dt_sim){
        // check if results need to be saved
        if (reportLeft < dt_sim){
            simulationSave(t, resDir);
            reportLeft = dt_report;
        }
        // check if rainfall need to be updated
        if (rainfallLeft < dt_sim){
            cout << "Reading Rainfall Record at " << rain_counter*dt_rainfall << "s" << endl;
            readHDF5(RainH5file, to_string(int(rain_counter*dt_rainfall)).c_str(), rain, "double");
            for (int i = 0; i < nCell; i++){
                // mm/h => m/s
                globalCell[i].S_p = rain[i] * globalCell[i].runoff_coef / 1000.0 / 3600.0;
            }
            rainfallLeft = dt_rainfall;
            rain_counter += 1;
        }
        // SWE solving
        if (steppingMethod == "RK-3"){
            solve_RK3(dt_sim);
        }
        else if (steppingMethod == "SDIRK4"){
            solve_SDIRK4(dt_sim, dtau, epsilon_LUSGS, maxiter_LUSGS);
        }
        else {
            cout << "Unknown Stepping Method" << endl;
        }

        cout << "Time now: " << t+dt_sim << endl;

        reportLeft -= dt_sim;
        rainfallLeft -= dt_sim;
    }

    return 0;
}