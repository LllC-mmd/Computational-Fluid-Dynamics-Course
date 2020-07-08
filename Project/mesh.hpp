#include <string>
#include <queue>
#include <array>
#include <vector>
#include <armadillo>
#include "variable.hpp"
#include "integration.hpp"

using namespace std;


class Node{
	public:
		int nId;
		double px;
		double py;
		double zb;
		// Mesh Topolpgy for Node
		// ---Global Index for Cells owning the Node; sizeof: 5~10
		vector<int> CellId_Node;
		// Node constructor
		Node();
		Node(int i, double x, double y, double ele);
		// Destructor
		~Node(){};
};

class Edge{
	private:
		// boundary condition
		double _hBC;
		double _uBC;
		double _huBC;
	public:
		int eId;
		// Mesh Topolpgy for Edge
		// ---Global Index for Nodes of the Edge: [0] fromNode; [1] toNode
		array<int, 2> nodeId_Edge;
		// ---Global Index for Cells possessing the Edge: [0] ownerCell; [1] neighbourCell;
		// ------Flux goes from the owner cell to the neighbor cell
		array<int, 2> cellId_Edge;
		double theta_n;
		double length;
		double zm;
		array<double, 2> z_Edge;
		// for MUSCL reconstruction e.g. slope computation
		double cx;
		double cy;
		int boundaryFlag;
		string BCsource;
		// Gaussian Points
		vector<U> UL_GP;
		vector<U> UR_GP;
		vector<FluxF> F_GP; 
		// Edge constructor
		Edge();
		Edge(int eIdx, int fId, int tId, int flag, char* source, U EdgeU, FluxF EdgeF);
		// get and set value
		int isFlood(int cId);
		int getFrom(){return nodeId_Edge[0];};
		int getTo(){return nodeId_Edge[1];};
		int getOCell(){return cellId_Edge[0];};
		int getNCell(){return cellId_Edge[1];};
		void setOCell(int p){cellId_Edge[0]=p;};
		void setNCell(int p){cellId_Edge[1]=p;};
		double getEdgeSpectralRadius(const U& UL, const U& UR);
		// Finite Volume Method
		// ---Reconstruction
		void reconstruct_BC(int alpha, double tol); 
		void reconstruct_free(int i); // reconstruction at the free outflow boundary
		void reconstruct_OpenBC_super(int i); // reconstruction at the open boundary for supercritical flow
		void reconstruct_solid(int i);  // reconstruction at the no-slip solid wall boundary
		void reconstruct_OpenBC_h(int i);  // reconstruction at the open boundary where h is given
		void reconstruct_OpenBC_u(int i);  // reconstruction at the open boundary where u is given
		void reconstruct_OpenBC_hu(int i, double tol);  // reconstruction at the open boundary where Q=hu is given
		void reconstruct_double_solid(int i, int alpha);  // reconstruction at the no-slip solid wall boundary (double)
		void reconstruct_adjust(const U& u_owner, const U& u_neighbour);
		// ---Flux computation
		void fluxCalc(double flux_dir);
		void fluxCalc_hllc(double flux_dir);  // determine the flux caused by source terms using HLLCS solver
		void fluxCalc_BC(double flux_dir);
		void fluxCalc_double_solid(double flow_dir);
		// Destructor
		~Edge(){};
};

class Cell{
	private:
		int _FloodedEdgeCount;
		U _Uc;
		Var _dU;
		U _u_inner;
		vector<U> _U_alpha;
		vector<Var> _Res_alpha;
	public:
		int cId;
		double cpx;
		double cpy;
		double cpz;
		// Mesh Topolpgy for Cell
		// ---Global Index for Nodes of the Cell
		array<int, 3> nodeId_Cell;
		// ---Global Index for Edges of the Cell 
		array<int, 3> edgeId_Cell;
		array<int, 3> edgeBCTag_Cell;
		// ---Global Index for Neighboring Cells of the Cell
		array<int, 3> neighborCellId_Cell;
		// Variational Reconstruction coefficients
		int vfv_order;  // vfv_order=3 for 4th reconstruction
		double dx_max;
		double dy_max;
		arma::Col<double> basis_mean;
		arma::Mat<double> A_eta;
		arma::Mat<double> A_hu;
		arma::Mat<double> A_hv;
		arma::Col<double> coef_eta;
		arma::Col<double> coef_hu;
		arma::Col<double> coef_hv;
		vector<arma::Mat<double> > B_eta;
		vector<arma::Mat<double> > B_hu;
		vector<arma::Mat<double> > B_hv;
		arma::Mat<double> b_eta_unit;
		arma::Mat<double> b_hu_unit;
		arma::Mat<double> b_hv_unit;
		// physical property
		array<double, 3> z_node;
		double runoff_coef;
		double S_p; // source term caused by net rainfall
		int label;
		double area;
		double n;
		double ix;   // slope along the x direction
		double iy;   // slope along the y direction
		// Cell constructor
		Cell();
		Cell(int cIdx, int e1, int e2, int e3, double n_r, double rc, U u);
		// get and set value
		int getFloodedEdgeCount() const {return _FloodedEdgeCount;};
		void setNeighbourCell(int cId, int eId){
			int it = 0;
			while (edgeId_Cell[it] != eId && it<edgeId_Cell.size()){
				it++;
			}
			neighborCellId_Cell[it] = cId;
		};
		int getNeighborCellCount(){
			int n = 0;
			for (int i = 0; i < neighborCellId_Cell.size(); i++){
				if (neighborCellId_Cell[i] > 0){
					n +=1;
				}
			}
			return n;
		}
		void FloodEdgeCount();
		U getValUc() const {return _Uc;};
		U getValU_inner() const {return _u_inner;};
		U getValU_alpha(int alpha) const {return _U_alpha[alpha];};
		Var getValdU() const {return _dU;};
		void setValUc(U u){_Uc = u;};
		void setValdU(Var du){_dU = du;}
		void setRes_alpha(int alpha){_Res_alpha[alpha] = getRes_VarRec(_U_alpha[alpha]);}
		void setValSp(double S){S_p=S;};
		void increValh_Uc(double delta_h){_Uc.h += delta_h;};
		void increValh_alpha(double delta_h, int alpha){_U_alpha[alpha].h += delta_h;};
		U VFR_set_eta(const U& u);
		double VFR_set_h(double eta);
		// MUSCL Reconstruction
		void CellRec_MUSCL(int alpha);
		// Variational Reconstruction
		void CellRec_VarRec(const U& u_c);
		double rec_hu(const U& u_c, double x1, double y1, double x2, double y2, double alpha);
		double rec_hv(const U& u_c, double x1, double y1, double x2, double y2, double alpha);
		double rec_eta(const U& u_c, double x1, double y1, double x2, double y2, double alpha);
		void VarRecInit(int Nc);
		void VarRecPrepare_basisMean();
		void VarRecPrepare_coefMat();
		void VarRecUpdate();
		// Reconstruction for dry and partially submerged cell
		void CellRec_dry(const U& u_c);
		void CellRec_partial(const U& u_c);
		// Use fluxes of edges to update the cell-averaged variable
		double getCellLocalTime(double dt);
		double getCellSpectralRadius();
		Var getRes_VarRec(const U& u_tmp);
		Var getResFlux_MUSCL(const U& u_tmp);
		void update_RK3(double dt, int alpha);
		void update_RK3_Init(){_U_alpha[0]=_Uc;};
		void update_RK3_Step();
		void WetDryCheck(const U& u_tmp);
		void update_Dry(const U& u_rhs, int alpha);
		void update_SemiFlooded(const U& u_rhs, int alpha);
		void update_Others(int alpha);
		void update_Others_Uc();
		// Implicit update
		void update_SDIRK4_Step(double dt);
		void update_SDIRK4_StoreU(int alpha){
			_U_alpha[alpha].h = _u_inner.h;
			_U_alpha[alpha].hu = _u_inner.hu;
			_U_alpha[alpha].hv = _u_inner.hv;
			_U_alpha[alpha] = VFR_set_eta(_U_alpha[alpha]);
		};
		void update_SDIRK4_ClearU(){
			for (int i = 0; i < _U_alpha.size(); i++){
        		_U_alpha[i].h = 0.0;
        		_U_alpha[i].hu = 0.0;
        		_U_alpha[i].hv = 0.0;
				_U_alpha[i] = VFR_set_eta(_U_alpha[i]);
    		}
		};
		void update_SDIRK4_InnerInit(){
			_u_inner.h = _Uc.h;
    		_u_inner.hu = _Uc.hu;
    		_u_inner.hv = _Uc.hv;
			_u_inner = VFR_set_eta(_u_inner);
			
			_dU.v1 = 0.0;
    		_dU.v2 = 0.0;
    		_dU.v3 = 0.0;
		};
		void update_SDIRK4_InnerUpdate();
		void update_LUSGS_forward(double dt, double dtau, int alpha);
		double update_LUSGS_backward(double dt, double dtau, int alpha);  
		// int inundation_check();  
		// Destructor
		~Cell(){};
};

void topSearch();