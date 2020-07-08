#include "postprec.hpp"
#include "global.hpp"
#include "mesh.hpp"

#include <fstream>
#include <iostream>
#include <iomanip>

using namespace std;


void simulationSave(double t_present, string result_Dir){
    ofstream outFile;
    string name = result_Dir + "/Res_" + to_string(t_present) + ".txt";
    cout << "Time: " << t_present << ", " << name << " recorded" << endl;
    
    outFile.open(name, ios::out);
    
    outFile << "cId" << "\t" << "h" << "\t" << "hu" << "\t" << "hv" << "\teta" << endl;
    for (int i = 0; i < nCell; i++){
        outFile << globalCell[i].cId << "\t" << setprecision(8) << globalCell[i].getValUc().h << "\t" 
        << setprecision(8) << globalCell[i].getValUc().hu << "\t" << setprecision(8) << globalCell[i].getValUc().hv << "\t" << setprecision(8) << globalCell[i].getValUc().eta << endl;
    }
    outFile.close();
};

void cellInfoSave(string result_Dir){
    ofstream outFile;
    string name = result_Dir + "/cellInfo.txt";
    
    outFile.open(name, ios::out);

    outFile << "CellID\tN1\tN2\tN3\tE1\tE2\tE3" << endl; 
    for (int i = 0; i < nCell; i++){
        outFile << globalCell[i].cId << "\t" <<  globalCell[i].nodeId_Cell[0] << "\t" << globalCell[i].nodeId_Cell[1] << "\t" << globalCell[i].nodeId_Cell[2] << "\t" 
        << globalCell[i].edgeId_Cell[0] << "\t" << globalCell[i].edgeId_Cell[1] << "\t" << globalCell[i].edgeId_Cell[2] << endl;
    }
}

void edgeInfoSave(string result_Dir){
    ofstream outFile;
    string name = result_Dir + "/edgeInfo.txt";
    
    outFile.open(name, ios::out);

    outFile << "EdgeID\tOCell\tNCell\ttheta_n\tBoundaryFlag" << endl; 
    for (int i = 0; i < nEdge; i++){
        outFile << globalEdge[i].eId << "\t" <<  globalEdge[i].getOCell() << "\t" << globalEdge[i].getNCell() << "\t" << globalEdge[i].theta_n/3.1415926*180 << "\t" << globalEdge[i].boundaryFlag << endl;
    }
}

void nodeInfoSave(string result_Dir){
    ofstream outFile;
    string name = result_Dir + "/nodeInfo.txt";
    
    outFile.open(name, ios::out);

    outFile << "NodeID\tX\tY" << endl; 
    for (int i = 0; i < nNode; i++){
        outFile << globalNode[i].nId << "\t" <<  globalNode[i].px << "\t" << globalNode[i].py << endl;
    }
}