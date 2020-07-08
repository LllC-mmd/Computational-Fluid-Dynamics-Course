#include <vector>
#include <string>

using namespace std;

class Node;
class Edge;
class Cell;

// global variable for geometry
extern int nNode;
extern int nEdge;
extern int nCell;

// global variable for FVM
extern double dt_min;
extern string recMethod;
extern string fluxMethod;
extern string steppingMethod;

extern double h_FLOOD;
extern double h_DRY;
extern double u_max;

// global variables are declared in main program, i.e., SWE_solver.cpp
extern vector<Node> globalNode;
extern vector<Edge> globalEdge;
extern vector<Cell> globalCell;