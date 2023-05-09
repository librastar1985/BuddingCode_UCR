

#ifndef SYSTEM_H_
#define SYSTEM_H_

#pragma once

//#include <gsl/gsl_matrix.h>
#include <fstream>
//#include <gsl/gsl_vector.h>
//#include <gsl/gsl_math.h>
//#include <gsl/gsl_multimin.h>
//#include <gsl/gsl_errno.h>

#include <memory>
#include <math.h>
#include <thrust/extrema.h>
#include <thrust/sort.h>
#include <thrust/copy.h>
#include <thrust/for_each.h>
#include <thrust/transform.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/tuple.h>
#include <thrust/device_vector.h>
#include <thrust/device_ptr.h>
#include <thrust/remove.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/copy.h>
#include <thrust/execution_policy.h>
#include <thrust/pair.h>
#include <stdint.h>
#include "my_vector.h"

#if !defined(ONED) && !defined(TWOD) && !defined(THREED)
#   define ONED
#   define TWOD
#   define THREED
#endif /* !defined(ONED) && !defined(TWOD) && !defined(THREED) */

enum {
#if defined(THREED)
	MAXD = 3
#elif defined(TWOD)
	MAXD = 2
#elif defined(ONED)
	MAXD = 1
#endif /* defined(THREED) */
};

#if !defined(RK_STEP)
enum {RK_STEP = 2};             // Should be 3rd or 4th order TVD RK
#endif
#if !defined(MAX_N_COEF)
enum {MAX_N_COEF = 3};          // 3: linear polynomial (2D)
#endif
#if !defined(N_EQN)
enum {N_EQN = 1};               // 4: correspond to  (2D) Euler
#endif
#if !defined(YES)
enum {YES = 1};
#endif
#if !defined(NO)
enum {NO = 0};                  // 1: const. 2D Burger's
#endif
#if !defined(MAX_AD_MV_ELE)
enum {MAX_AD_MV_ELE = 4};
#endif
#if !defined(MAX_WL_AD_ELE)
enum {MAX_WL_AD_ELE = 4};
#endif

// struct loc_EleNumber: public std::binary_function<ELEMENT, int, bool>
// {
//   bool operator () ( const ELEMENT &ele, const int &number ) const {
//     return ele.Id == number;
//     }
// };
//
struct CapsidInfoVecs {
	thrust::device_vector<int> id_bucket;	//bucket id
	// bucket value means global rank of a certain point
	thrust::device_vector<int> id_value;//node id

 	thrust::device_vector<double> nodeLocX;
	thrust::device_vector<double> nodeLocY;
	thrust::device_vector<double> nodeLocZ;


	//used for capside linear springs binding to membrane
	thrust::device_vector<int> tempMembraneId;
	thrust::device_vector<int> tempCapsideId;

	thrust::device_vector<double> tempLengthsPairs;
	thrust::device_vector<double> tempNodeForceX;
	thrust::device_vector<double> tempNodeForceY;
	thrust::device_vector<double> tempNodeForceZ;

	//used for capside repulsion
	int factor = 10;//number of default nodes repulsion can interact with.
	thrust::device_vector<int> tempNodeIdUnreduced;
	thrust::device_vector<double> tempNodeForceXUnreduced;
	thrust::device_vector<double> tempNodeForceYUnreduced;
	thrust::device_vector<double> tempNodeForceZUnreduced;

	thrust::device_vector<int> tempNodeIdReduced;
	thrust::device_vector<double> tempNodeForceXReduced;
	thrust::device_vector<double> tempNodeForceYReduced;
	thrust::device_vector<double> tempNodeForceZReduced;

	double spring_constant = 10.0;
	double length_zero = 0.97;
	double length_cutoff = 1.63;
	int maxNodeCount;
	double viscosity = 1.0;
	double forceX = 0.0;
	double forceY = 0.0;
	double forceZ = 0.0;

	int num_connections=0;

};

//Data Structure for node location. velocity and force
struct CoordInfoVecs {

	double k_0;
	double k_1;
	double k_2;
	double k_3;
	double k_4;
	double k_ss;//10.75;
	double beta;///1.45;
	double gamma;
	double q1;
	double h;

	thrust::device_vector<double> b_per_triangle;

	thrust::device_vector<double> u;

	thrust::device_vector<double> soln_per_triangle;

	thrust::device_vector<double> scaling_per_edge;

	thrust::device_vector<int> nodes2Triangles_1;
	thrust::device_vector<int> nodes2Triangles_2;
	thrust::device_vector<int> nodes2Triangles_3;
	thrust::device_vector<int> nodes2Triangles_4;
	thrust::device_vector<int> nodes2Triangles_5;
	thrust::device_vector<int> nodes2Triangles_6;
	thrust::device_vector<int> nodes2Triangles_7;
	thrust::device_vector<int> nodes2Triangles_8;
	thrust::device_vector<int> nodes2Triangles_9;
//	thrust::device_vector<int> nodes2Triangles_10;
//	thrust::device_vector<int> nodes2Triangles_11;
//	thrust::device_vector<int> nodes2Triangles_12;

	thrust::device_vector<double> SurfaceNormalX;
	thrust::device_vector<double> SurfaceNormalY;
	thrust::device_vector<double> SurfaceNormalZ;

	thrust::device_vector<int> nndata1;
	thrust::device_vector<int> nndata2;
	thrust::device_vector<int> nndata3;
	thrust::device_vector<int> nndata4;
	thrust::device_vector<int> nndata5;
	thrust::device_vector<int> nndata6;
	thrust::device_vector<int> nndata7;
	thrust::device_vector<int> nndata8;
	thrust::device_vector<int> nndata9;
//	thrust::device_vector<int> nndata10;
//	thrust::device_vector<int> nndata11;
//	thrust::device_vector<int> nndata12;

	thrust::device_vector<bool> isNodeFixed;
	//GLOBAL COORDS
	// X,Y,Z, location, velocity and force of all nodes
	thrust::device_vector<double> prevNodeLocX;
	thrust::device_vector<double> prevNodeLocY;
	thrust::device_vector<double> prevNodeLocZ;

	thrust::device_vector<double> prevNodeForceX;
	thrust::device_vector<double> prevNodeForceY;
	thrust::device_vector<double> prevNodeForceZ;

 	thrust::device_vector<double> nodeLocX;
	thrust::device_vector<double> nodeLocY;
	thrust::device_vector<double> nodeLocZ;

	thrust::device_vector<double> nodeForceX;
	thrust::device_vector<double> nodeForceY;
	thrust::device_vector<double> nodeForceZ;

	thrust::device_vector<double> tempNodeForceXReduced;
	thrust::device_vector<double> tempNodeForceYReduced;
	thrust::device_vector<double> tempNodeForceZReduced;
	thrust::device_vector<double> tempNodeForceXUnreduced;
	thrust::device_vector<double> tempNodeForceYUnreduced;
	thrust::device_vector<double> tempNodeForceZUnreduced;

	//LOCAL COORDS
	//indices of each triangle
	int num_triangles;
	thrust::device_vector<int> triangles2Nodes_1;
	thrust::device_vector<int> triangles2Nodes_2;
	thrust::device_vector<int> triangles2Nodes_3;


	//indices of each edge
	int num_edges;
	thrust::device_vector<int> edges2Nodes_1;
	thrust::device_vector<int> edges2Nodes_2;

	//indices of 2 triangle on each edge
	thrust::device_vector<int> edges2Triangles_1;
	thrust::device_vector<int> edges2Triangles_2;

	//indices of edges on each triangle.
	thrust::device_vector<int> triangles2Edges_1;
	thrust::device_vector<int> triangles2Edges_2;
	thrust::device_vector<int> triangles2Edges_3;

	thrust::device_vector<int> triangles2Triangles_1;
	thrust::device_vector<int> triangles2Triangles_2;
	thrust::device_vector<int> triangles2Triangles_3;

};



//struct used for linking of nodes in network
struct AuxVecs {
	// bucket key means which bucket ID does a certain point fit into
	thrust::device_vector<int> id_bucket;	//bucket id
	// bucket value means global rank of a certain point
	thrust::device_vector<int> id_value;//node id
	// bucket key expanded means what are the bucket IDs are the neighbors of a certain point
	thrust::device_vector<int> id_bucket_expanded;
	// bucket value expanded means each point ( represented by its global rank) will have multiple copies
	thrust::device_vector<int> id_value_expanded;

	// begin position of a keys in id_bucket_expanded and id_value_expanded
	//entry keyBegin[bucketKey] returns start of indices to link
	thrust::device_vector<int> keyBegin;
	// end position of a keys in id_bucket_expanded and id_value_expanded
	thrust::device_vector<int> keyEnd;

	int endIndexid_bucket;
};



struct DomainParams {
	double minX;
	double maxX;
	double minY;
	double maxY;
	double minZ;
	double maxZ;
	double originMinX;
	double originMaxX;
	double originMinY;
	double originMaxY;
	double originMinZ;
	double originMaxZ;
	double gridSpacing = 1.5;//bucket scheme search distance, must be larger than cutoff for capsid
	int XBucketCount;
	int YBucketCount;
	int ZBucketCount;
	int totalBucketCount = 0;
};

struct LJInfoVecs{
	double LJ_PosX;
	double LJ_PosY;
	double LJ_PosZ;
	thrust::device_vector<double> LJ_PosX_all;
	thrust::device_vector<double> LJ_PosY_all;
	thrust::device_vector<double> LJ_PosZ_all;



	double Rmin_M=0.97;//1.0;
	double Rcutoff_M=0.97;//1.0;
	double Rmin_LJ;
	double Rcutoff_LJ;

	double epsilon_M=1.0;
	double epsilon_M_att1;
	double epsilon_M_att2;
	double epsilon_M_rep1;
	double epsilon_M_rep2;
	double epsilon_LJ;
	double epsilon_LJ_rep1;
	double epsilon_LJ_rep2;
	double spring_constant;

	thrust::device_vector<int> node_id_close;
	double lj_energy_M;
	double lj_energy_LJ;
	double forceX;
	double forceY;
	double forceZ;
	thrust::device_vector<double> forceX_all;
	thrust::device_vector<double> forceY_all;
	thrust::device_vector<double> forceZ_all;

};

struct AreaTriangleInfoVecs {
	double dummy;
	int factor = 3;//used for reduction
	double initial_area = 0.0048013;
	double spring_constant;
	double spring_constant_weak;

	double area_triangle_energy;

	thrust::device_vector<int> tempNodeIdUnreduced;
	thrust::device_vector<double> tempNodeForceXUnreduced;
	thrust::device_vector<double> tempNodeForceYUnreduced;
	thrust::device_vector<double> tempNodeForceZUnreduced;

	thrust::device_vector<int> tempNodeIdReduced;
	thrust::device_vector<double> tempNodeForceXReduced;
	thrust::device_vector<double> tempNodeForceYReduced;
	thrust::device_vector<double> tempNodeForceZReduced;

};

struct BendingTriangleInfoVecs {
	int numBendingSprings=0;

	int factor = 4;//used for reduction
	double spring_constant;
	double spring_constant_weak;
	double spring_constant_raft;
	double spring_constant_coat;
	double initial_angle = 0.0;//radians
	double initial_angle_raft;
	double initial_angle_coat;
	double initial_angle_bud;

	double bending_triangle_energy;

	thrust::device_vector<int> tempNodeIdUnreduced;
	thrust::device_vector<double> tempNodeForceXUnreduced;
	thrust::device_vector<double> tempNodeForceYUnreduced;
	thrust::device_vector<double> tempNodeForceZUnreduced;

	thrust::device_vector<int> tempNodeIdReduced;
	thrust::device_vector<double> tempNodeForceXReduced;
	thrust::device_vector<double> tempNodeForceYReduced;
	thrust::device_vector<double> tempNodeForceZReduced;
};
struct LinearSpringInfoVecs {

	int factor = 2;//used for reduction
	double spring_constant;
	double spring_constant_weak;
	double spring_constant_att1;
	double spring_constant_att2;
	double spring_constant_rep1; //This is the "D" in Morse potential
	double spring_constant_rep2; //This is the "a" in Morse potential

	double linear_spring_energy;
	double memrepulsion_energy;
	double scalar_edge_length;
	
	thrust::device_vector<double> edge_initial_length;

	thrust::device_vector<int> tempNodeIdUnreduced;
	thrust::device_vector<double> tempNodeForceXUnreduced;
	thrust::device_vector<double> tempNodeForceYUnreduced;
	thrust::device_vector<double> tempNodeForceZUnreduced;

	thrust::device_vector<int> tempNodeIdReduced;
	thrust::device_vector<double> tempNodeForceXReduced;
	thrust::device_vector<double> tempNodeForceYReduced;
	thrust::device_vector<double> tempNodeForceZReduced;
};

struct GeneralParams{
	std::vector<int> edge_undergoing_growth;
	bool nonuniform_wall_weakening;
	std::vector<int> triangle_undergoing_growth;
	double chemdiff_time_step_size;
	int current_total_sim_step;
	double chemdiff_max_step;
	double current_bud_area;
	double kT;
	double kT_growth;
	double tau;
	int solve_time=100;
	double Rmin = 1.0; //Current mesh minimum edge length, subject to change.
	double Rmin_growth;
	double abs_Rmin;
	int iteration = 0;
	int maxNodeCount;
	int maxNodeCountLJ;
	double length_scale;
	//parameters for advancing timestep and determining equilibrium

	double dt;
	double nodeMass = 1.0;
	int edges_in_upperhem_list_length;
	double growth_energy_scaling;
	thrust::device_vector<int> edge_to_ljparticle;
	thrust::device_vector<int> nodes_in_tip;
	thrust::device_vector<int> nodes_in_upperhem;
	thrust::device_vector<int> edges_in_upperhem;
	thrust::device_vector<int> edges_in_upperhem_list;
	thrust::device_vector<int> edges_in_tip;
	thrust::device_vector<int> triangles_in_tip;
	thrust::device_vector<int> triangles_in_upperhem;
	thrust::device_vector<int> boundaries_in_upperhem;
	thrust::device_vector<double> angle_per_edge;
	double centerX = 0.0;
	double centerY = 0.0;
	double centerZ = 0.0;

	double current_total_volume;
	double true_current_total_volume;
	double eq_total_volume;
	double volume_spring_constant;
	double volume_energy;
	double eq_total_boundary_length;
	double line_tension_energy;
	double line_tension_constant;
	double safeguardthreshold;
	

	//int num_of_triangles;
	//int num_of_edges;
	int true_num_edges;

	double insertion_energy_cost;
	double strain_threshold;

	int SCALE_TYPE;
	double scaling_pow;
	double gausssigma;
	double hilleqnconst;
	double hilleqnpow;

	thrust::device_vector<int> no_weakening;
	double septin_ring_z;
	double boundary_z;

};


class Storage;
class SystemBuilder;
struct HostSetInfoVecs;

typedef struct {
	int vertex[2];      //nodes of each edge
	int adj_ele[2];     //the elements who have the edge as the common edge
	double L0;          //initial length
	double cost0;       //cosine of the const of the theta0
	double sint0;
	double theta0;
}LINK;

typedef struct {
        int     Id;              // id of the element
	int     vertex[3];       //three nodes, anti-clockwise in out-normal direction ----- save index of nodes
	int     adj_ele[3];      //three adjacent element, rule of the name:  
                                 //adj_ele[0]: the adjacent triangle which has the common edge E(1,2) which has the 
	                         //node 1 and 2
	int		adj_ele_upperhem[3]; //Will indicate whether the adjacent elements are in upperhem or not (i.e. in budding region or not)
	int 	upperhem; //Will indicate whether the element is in upperhem or not.
	double  area,area0;         //area of each element 
        double  center[3];          // NEW::::coordinates of cell center. 
        double  length_side[3];
        CVector out_norm;           //unit out-normal vector
        int     N_ad_Mol;           // Number of adhesion receptor on the element.
        int     ad_on[MAX_AD_MV_ELE];     // For each MV, its status (whether adhere to substrate) is recorded here
                                          // MAX_AD_MV_ELE > N_ad_MV
        double  bond_len[MAX_AD_MV_ELE];  // length of each MV

        /*for DG */
        CVector  side_vector[3];   // side vectors v[i] = p[i+1]-p[i] in global coord.
        CVector  side_norm[3];     // side vectors v[i]\cdot side_norm[i] = 0. outward conormal in global coord. 
        double   aff_conorm[3][4];          // outward conormal in its own (current) local crds.
        double   aff_adj_ele_conorm[3][4];  // outward conormals of  adj elements in local crds of current element. 
        double   affine_pt[3][4], aff_cent[4]; // crds of vertices and centriod in local crds.

        double   dg_rho[RK_STEP][MAX_N_COEF];      
        double   dg_q[RK_STEP][2][MAX_N_COEF];            // LDG: flux q for diffusion problem.
		double	 b[RK_STEP][MAX_N_COEF]; 						//This array is for the inhibitor of u (represented as dg_rho).
		double   u_for_dg[1];

        // The following saved fluxes are used for every case.
        double   u_hat_flux_store[3][4][2];                        //side:point:saved_flux_vector
        bool     uhat_f_flag[3];
        double   q_hat_flux_store[3][4][2];                        //side:point:saved_flux_vector
        bool     qhat_f_flag[3];
        
        // The following saved fluxes are for solving CH eq. 
        double   r_hat_flux_store[3][4][2];                        //side:point:saved_flux_vector
        bool     rhat_f_flag[3];

        double   **mass_matrix, **mass_inv;               // for mass matrix on whole cell. same as tri->Bmass_matrix. tri->Bmass_inv. 
        double   **vec_vh_mass_matrix, **vec_vh_mass_inv; // for mass matrix of vector basis
        double   **affine_trans, **inv_affine_trans;      // affine transformation matrix.

        /*variables for LDG for solving Cahn-Hilliard*/
        double   dg_S[RK_STEP][2][MAX_N_COEF];   // LDG: variable S.
        double   dg_P[RK_STEP][2][MAX_N_COEF];   // LDG: variable P.
        double   dg_Q[RK_STEP][MAX_N_COEF];      // LDG: variable q = \gamma \nabla\cdot W.
        double   dg_W[RK_STEP][2][MAX_N_COEF];   // LDG: variable W.
        double   dg_r[RK_STEP][MAX_N_COEF];      // LDG: variable r: derivative of hydrophobic part of free energy.
}ELEMENT;

typedef struct {
	int NNum;  // number of nodes from biconcavefine.neu file
	int ENum;  // number of elements from biconcavefine.neu file
	
	std::list<ELEMENT>   *Curr_ele; // elements in interior of subdomain. 
                                   // cur_ele_it=rbc->Curr_ele->begin(). 
                                   //  NOTE::: cur_ele_it has its own memory storage and IS NOT SAME as ele[i].
	std::list<int>	*Curr_ele_index;

	std::list<ELEMENT>::iterator Curr_it;
	std::list<int>::iterator     Curr_ind_it;

        std::list<ELEMENT>::iterator **sv_it_adj_ele; // use to save iterator of adj. ele. 
		
			
	//ELEMENT*	Buf_ele[MAXD][2]; // save elements in buffer zone. 
	//int*	        Buf_ele_index[MAXD][2];
	//int 	        Buf_ele_num[MAXD][2];


	std::list<int>       *cur_node_index;
        std::list<int>       *cur_node_priority; // node inside unbuffered subdomain has value 0.
                                            // Receive node (in the buffer zone) has value 1???
	std::list<int>       *curr_link_index;
	int             *node_nb;

	std::list<int>	*EleNode;  /// Elements have EleNode[]; indices of the elements sharing the node
                                   /// EleNode = new list<int>[rbc->NNum];
	ELEMENT 	*ele;			 //element ele[ENum]. ele[i] 
	// CVector		  node[NNum];	 //node
	CVector 	*node;			 //node
        CVector         *normal_node_fitting; // unit normal at the node by fitting. 
	// double		  nw[NNum]; 	 //area of each element
	double		*nw;   //1/3 of area of triangles surrounding a node. See Comp_ele()
	CVector		*la[2]; // Lagrange multiplier
	LINK		*Link;	 //edge
	double		Area0, Vol0, Vol_target;  // Vol0: initial volume. Vol_target: final volume
        double          X0;                       // average link length
        double          *Mean_curv; // Total mean curvature at nodes, which is (k1+k2). 
        double          *Gauss_curv; // Gauss curvature at nodes 
        /*See paper Discretizing Laplace-Beltrami operator from differential quantities*/
        double          **geom_g; // determinant and derivative of determinant of g_[u,v]
        double          **inv_geom_uv; // inverse of g_[u,v]
        double          **der_inv_geom_uv; // derivative of entry of inverse of g_[u,v]
        //FUND_GEOM       *fund_geom;        

        double          *lambda_inex; /*Lagrange multiplier of local inextensibility*/
        double          lambda_vol;  /*Lagrange multiplier of Volume change. Use to simulate
                                      sphere deforming to biconcave red blood cell, during which
                                      the volume changes. */
	int LNum;		    // Number of edge

	int k;			 /// This only affects function void LoadMesh()

	double		tau,dx,dy,dz,dt,cs2,gra;

	CVector         Center;
	CVector         *nodeloc;
        
        CVector         *fn;
        double          *fMarker; // Mark the node where the force is added
        double          *chem;       // concentration on the node
        int             lNMark;
        int             rNMark;
        int             uNMark;
        int             bNMark;
        CVector         *node0;   // X^{n}
        CVector         *nodep;   // X^{n-1}
        CVector         *nodeVel; 

	//MESHP           mesh_para;
//	RECT_GRID	*wall_grid;
        int             CELL_ID; 
        char            Mesh_file_name[128]; 
        double          avg_radi; //average radius
} RBC;

struct _Front {
	POINTER *RBC_pointer;
	int Num_RBCs; 
};
typedef struct _Front Front;



class System {
public:
	std::weak_ptr<SystemBuilder> weak_bld_ptr;
	GeneralParams generalParams;
	DomainParams domainParams;
	AuxVecs auxVecs;
	CoordInfoVecs coordInfoVecs;

	CapsidInfoVecs capsidInfoVecs;
	LinearSpringInfoVecs linearSpringInfoVecs;
	BendingTriangleInfoVecs bendingTriangleInfoVecs;
	AreaTriangleInfoVecs areaTriangleInfoVecs;
	LJInfoVecs ljInfoVecs;
	RBC rbc;
	ELEMENT element;

	std::shared_ptr<Storage> storage;

	//gsl_vector* df;
	//gsl_vector* locations;
    



public:

	System();

	void set_weak_builder(std::weak_ptr<SystemBuilder> _weak_bld_ptr);

	void PrintForce();

	void initializeSystem(HostSetInfoVecs& hostSetInfoVecs);

	void assignStorage(std::shared_ptr<Storage> _storage);
 
	void solveSystem();

	void setExtras();

	void setBucketScheme();

	//void Solve_Forces(const gsl_vector* temp_locations);
	void Solve_Forces();
	//double Solve_Energy(const gsl_vector* temp_locations);

	//void dev_to_gsl_loc_update(gsl_vector* temp_locations);
	//void gsl_to_dev_loc_update(const gsl_vector* temp_locations);
	

};



#endif /*POLYMERSYSTEM_H_*/
