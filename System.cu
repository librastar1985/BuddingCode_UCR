#include "System.h"
#include "SystemStructures.h" 
#include "AreaTriangles.h"
#include "AreaTrianglesEnergy.h"
#include "BendingTriangles.h"
#include "BendingTrianglesEnergy.h"
#include "MemRepulsionSprings.h"
#include "MemRepulsionEnergy.h"
#include "LinearSprings.h"
#include "LinearSpringsEnergy.h"
#include "LJSprings.h"
#include "LJSprings_LJ.h"
#include "NodeAdvance.h"
#include "BucketScheme.h"
#include "Storage.h" 
#include "Edgeswap_test.h"
#include "SystemBuilder.h"
#include <vector>
#include "VolumeComp.h"
#include "VolumeSprings.h"
#include <bits/stdc++.h>
#include "LineTensionSprings.h"
#include "Growth.h"

///////////////////////////////////////////////////////////////////
///////////////////////// WARNING ////////////////////////////////
//////////////////REMEMBER TO CHANGE THE /////////////////////////
/////////////EQUILIBRIUM LENGTH OF EACH TRIANGLE EDGE /////////////
//////////////// IN THE VECTOR INITIALIZATION ////////////////////
//////////////////SECTION TOWARD THE END OF THE CODE /////////////
////////////////////////////////////////////////////////////////////

 //somehow the gradient is not being set in my version

//bool IsPos (int i){return (i>=0);}
int count_bigger(const std::vector<int>& elems) {
    return std::count_if(elems.begin(), elems.end(), [](int c){return c >= 0;});
}

System::System() {};

void System::Solve_Forces(){

	thrust::fill(coordInfoVecs.nodeForceX.begin(), coordInfoVecs.nodeForceX.end(), 0.0);
	thrust::fill(coordInfoVecs.nodeForceY.begin(), coordInfoVecs.nodeForceY.end(), 0.0);
	thrust::fill(coordInfoVecs.nodeForceZ.begin(), coordInfoVecs.nodeForceZ.end(), 0.0);
	
	//setBucketScheme();
	ComputeLinearSprings( 
		generalParams, 
		coordInfoVecs,
		linearSpringInfoVecs, 
		ljInfoVecs);
	
	ComputeAreaTriangleSprings(
		generalParams,
		coordInfoVecs,
		areaTriangleInfoVecs);

	ComputeCosTriangleSprings(
		generalParams,
		coordInfoVecs,  
		bendingTriangleInfoVecs); 
	
	ComputeMemRepulsionSprings(
		coordInfoVecs,
		linearSpringInfoVecs, 
		capsidInfoVecs,
		generalParams,
		auxVecs);

	ComputeVolume(
		generalParams,
		coordInfoVecs,
		linearSpringInfoVecs,
		ljInfoVecs);

	ComputeVolumeSprings(
		coordInfoVecs,
		linearSpringInfoVecs, 
		capsidInfoVecs,
		generalParams,
		auxVecs);

	ComputeLineTensionSprings(
		generalParams,
		coordInfoVecs,
		linearSpringInfoVecs);
		
};


void System::solveSystem() {
	
	generalParams.safeguardthreshold = 12;
	//safeguardthreshold is the maximum number of neighboring nodes a node can have.

	//////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////
	////////////////////////// PARAMETER SETTINGS ////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////

	double Max_Runtime = 0.0005;
	double Max_RunStep = Max_Runtime/generalParams.dt;
	std::cout<<"Max runtime = "<<Max_Runtime<<std::endl;
	std::cout<<"Max runstep = "<<Max_RunStep<<std::endl;
	bool runSim = true;
	int num_edge_loop;
	double initial_kT;
	initial_kT = generalParams.kT;//This is for the acceptance of change after looping through every edge within proximity.
	double SAMPLE_SIZE = 0.01;
	std::cout<<"Sample ratio: "<<SAMPLE_SIZE<<std::endl;
	//This determines the number of edges to test for bondflip remeshing

	auto edgeswap_ptr = std::make_shared<Edgeswap>(coordInfoVecs, generalParams);
	int RECORD_TIME = 100;//round(Max_RunStep/2);
	std::cout<<"Record frequency = "<<RECORD_TIME<<std::endl;
	//int GROWTH_TIME = 1;
	//std::cout<<"Growth frequency = "<<GROWTH_TIME<<std::endl;
	int NKBT = 10000; //The max number of edge-swap attempt per kBT value
	std::cout<<"Number of edge-swap per kBT value (or total number of edge-swap if kBT is fixed)= "<<NKBT<<std::endl;
	double min_kT = 0.00000000001;//0.21;
	std::cout<<"min kT for simulation termination = "<<min_kT<<std::endl;
	int WHEN = 0;
	double old_total_energy = 0.0;
	double new_total_energy = 0.0;
	double energy_gradient = 0.0;
	int Num_of_step_run = 0;
	auto build_ptr = weak_bld_ptr.lock();//upgrade weak builder to access host variables.
	//std::cout<<"initial LJ-x : "<< ljInfoVecs.LJ_PosX <<std::endl;
	//std::cout<<"initial LJ-y : "<< ljInfoVecs.LJ_PosY <<std::endl;
	//std::cout<<"initial LJ-z : "<< ljInfoVecs.LJ_PosZ <<std::endl;
		

    
	double min_energy;
	generalParams.true_num_edges = 0;
	for (int i = 0; i < coordInfoVecs.num_edges; i++){
		if (coordInfoVecs.edges2Nodes_1[i] != INT_MAX && coordInfoVecs.edges2Nodes_2[i] != INT_MAX){
			generalParams.true_num_edges += 1;
		}
	}
	
	/////////////////////////////////////////////////////////////////
	/////////////////////// MEMBRANE RELATED ////////////////////////
	/////////////////////////////////////////////////////////////////
	
	
	double VOLUME_FACTOR = 1.0;
	//VOLUME_FACTOR determines the target volume which equals to VOLUME_FACTOR*initial_volume.
	//double tip_depth = 0.5;
	//tip_depth is currently unused.
	
	double weakened = 6.0;
	//weakened determines the minimum height of the z-coordinate of the membrane node to be considered in the area of weakened mechanical properties.
	//double tip_base = 6.0;
	//tip_base currently unused.

	int translate_frequency = 5;
	//translate_frequency determines the frequency for the mesh to re-center and perform dynamical remeshing

	double EXPAN_THRESHOLD = 2.0;
	double EXPAN_THRESHOLD_weak = 1.5;
	std::cout<<"EXPANSION THRESHOLD = "<<EXPAN_THRESHOLD<<std::endl;
	int RULES_OF_EXPAN = 3;	//EXPAN_THRESHOLD is the yielding ratio where a pair of triangles will perform expansion.
	
	std::cout<<"EXPANSION RULE = "<<RULES_OF_EXPAN<<std::endl;
	//EXPAN_THRESHOLD_weak is the secondary yielding ratio.
	//RULES_OF_EXPAN controls how the EXPAN_THRESHOLD is applied:
	// 1:= Both trianglular areas must exceed the threshold value.
	// 2:= If one trianglular area exceeds the treshold value while the other exceeds the secondary threshold value.
	// 3:= If the combined area of the two triangles exceed 2*EXPAN_THRESHOLD.

	for (int i = 0; i < generalParams.maxNodeCount; i++){
		generalParams.centerX += coordInfoVecs.nodeLocX[i];
		generalParams.centerY += coordInfoVecs.nodeLocY[i];
		generalParams.centerZ += coordInfoVecs.nodeLocZ[i];
	}
	generalParams.centerX = generalParams.centerX/generalParams.maxNodeCount;
	generalParams.centerY = generalParams.centerY/generalParams.maxNodeCount;
	generalParams.centerZ = generalParams.centerZ/generalParams.maxNodeCount;
	double displacementX, displacementY, displacementZ;
	double newcenterX, newcenterY, newcenterZ;
	//centerX, centerY, centerZ is determined as the referenced origin for recentering of the mesh.

	std::vector<int> VectorShuffleForGrowthLoop;
	std::vector<int> VectorShuffleForFilamentLoop;

	//double max_height = coordInfoVecs.nodeLocZ[35];
	//double min_height = coordInfoVecs.nodeLocZ[38];
	//Max and min height of the membrane nodes, these have to be changed if the mesh used is changed.

	generalParams.Rmin = 1.0;//0.1505;
	//Equilibrium length of an edge of the triangle.
	generalParams.abs_Rmin = 1.0;//0.1505;//0.586955;
	//Equilibrium distance between membrane node for volume exclusion.
	areaTriangleInfoVecs.initial_area =0.433;// 0.0098;//0.0048013;
	//Equilibrium triangular area.
	ljInfoVecs.Rmin_M = 0.0;
	//Equilibrium distance between the nucleus particle and membrane.
	ljInfoVecs.Rcutoff_M = 0.0;
	//Maximal interaction range between the nucleus and membrane.
	ljInfoVecs.Rmin_LJ = 0.0;//3.0//1.0;
	//Equilibrium distance between nuclei.
	ljInfoVecs.Rcutoff_LJ = 0.0;//3.0;//1.0;
	//Maximal interaction range between the nuclei.
	ljInfoVecs.epsilon_M_att1 = 0.0;//6.0;//16.0;
	ljInfoVecs.epsilon_M_att2 = 0.0;//1.0;//1.0;
	std::cout<<"Morse_NM_D_att = "<<ljInfoVecs.epsilon_M_att1<<std::endl;
	std::cout<<"Morse_NM_a_att = "<<ljInfoVecs.epsilon_M_att2<<std::endl;
	//Coefficient for the attractive interaction between nuclei and membrane.
	ljInfoVecs.epsilon_M_rep1 = 0.0;//12.5;//16.0;
	ljInfoVecs.epsilon_M_rep2 = 0.0;//0.5;//1.0;
	std::cout<<"Morse_NM_D_rep = "<<ljInfoVecs.epsilon_M_rep1<<std::endl;
	std::cout<<"Morse_NM_a_rep = "<<ljInfoVecs.epsilon_M_rep2<<std::endl;
	//Coefficient for the repulsive interaction between nuclei and membrane.
	
	ljInfoVecs.epsilon_LJ_rep1 = 0.0;//10.0;//0.5;// 0.06;//7.5;
	ljInfoVecs.epsilon_LJ_rep2 = 0.0;//0.5;//1.0;//1.0;//1.0;
	std::cout<<"Morse_NN_D = "<<ljInfoVecs.epsilon_LJ_rep1<<std::endl;
	std::cout<<"Morse_NN_a = "<<ljInfoVecs.epsilon_LJ_rep2<<std::endl;
	//Coefficient of the interaction between nuclei.

	linearSpringInfoVecs.spring_constant_rep1 = 0.01;
	linearSpringInfoVecs.spring_constant_rep2 = 9.0;
	std::cout<<"Membrane volume exclusion Morse D = "<<linearSpringInfoVecs.spring_constant_rep1<<std::endl;
	std::cout<<"Membrane volume exclusion Morse a = "<<linearSpringInfoVecs.spring_constant_rep2<<std::endl;
	//The coefficient used for non-neighboring membrane node volume exclusion.
	//rep1 is the "D" and rep2 is the "alpha" in the standard form of Morse potential.

	generalParams.volume_spring_constant = 52.5;
	std::cout<<"spring constant for volume conservation = "<<generalParams.volume_spring_constant<<std::endl;
	generalParams.line_tension_constant = 200.0;
	std::cout<<"spring constant for the septin ring = "<<generalParams.line_tension_constant<<std::endl;
	generalParams.length_scale = 0.8;//0.1577;//1.0*generalParams.Rmin;// 0.8333;
	std::cout<<"equilibrium length of each segment of the septin ring = "<<generalParams.length_scale<<std::endl;

	double scale_linear = 75.0/9.0;//75.0/15.0;
	double scale_bend = 75.0/20.0;//75.0/7.5;
	double scale_area = 75.0/5.0;//75.0/15.0;
	std::cout<<"scaling of different region linear = "<<scale_linear<<std::endl;
	std::cout<<"scaling of different region bend = "<<scale_bend<<std::endl;
	std::cout<<"scaling of different region area = "<<scale_area<<std::endl;
	linearSpringInfoVecs.spring_constant_weak = linearSpringInfoVecs.spring_constant/scale_linear;
	bendingTriangleInfoVecs.spring_constant_weak = bendingTriangleInfoVecs.spring_constant/scale_bend;
	areaTriangleInfoVecs.spring_constant_weak = areaTriangleInfoVecs.spring_constant/scale_area;
	//Scaling of the weakend mechanical properties.

	bendingTriangleInfoVecs.initial_angle =  0.0906;//0.17549;//0.15;//0.0906;
	bendingTriangleInfoVecs.initial_angle_raft =  0.0906;//0.17549;//0.15;
	bendingTriangleInfoVecs.initial_angle_coat = 0.0906;//0.17549;//0.15;//0.167448079;
	std::cout<<"equilibrium bending angle of the membrane = "<<bendingTriangleInfoVecs.initial_angle<<std::endl;
	//raft and coat are current unused due to the assumption of uniform preferred curvature.
	
	bendingTriangleInfoVecs.spring_constant_raft = 0.0;//bendingTriangleInfoVecs.spring_constant;
	bendingTriangleInfoVecs.spring_constant_coat = 0.0;//bendingTriangleInfoVecs.spring_constant;
	bendingTriangleInfoVecs.spring_constant = bendingTriangleInfoVecs.spring_constant*(2/sqrt(3));
	bendingTriangleInfoVecs.spring_constant_raft = bendingTriangleInfoVecs.spring_constant_raft*(2/sqrt(3));
	bendingTriangleInfoVecs.spring_constant_coat = bendingTriangleInfoVecs.spring_constant_coat*(2/sqrt(3));
	std::cout<<"Effective bending coefficient is calculated by multiplying 2/sqrt(3)"<<std::endl;
	std::cout<<"effective bending coefficient of the membrane = "<<bendingTriangleInfoVecs.spring_constant<<std::endl;
	std::cout<<"effective bending coefficient of the membrane raft = "<<bendingTriangleInfoVecs.spring_constant_raft<<std::endl;
	std::cout<<"effective bending coefficient of the membrane coat = "<<bendingTriangleInfoVecs.spring_constant_coat<<std::endl;

	

	/////////////////////////////////////////////////////////////////
	////////////////// END OF MEMBRANE RELATED //////////////////////
	/////////////////////////////////////////////////////////////////


	/////////////////////////////////////////////////////////////////
	//////////////////////// NULCEUS RELATED ////////////////////////
	/////////////////////////////////////////////////////////////////
	double beta1 = 0.0;
	double beta2 = 0.0;
	std::cout<<"manual push speed for the nucleus tip = "<<beta1<<std::endl;
	std::cout<<"manual push speed for the remainder of the nucleus = "<<beta2<<std::endl;
	//beta1 is the vertical speed (0, 0, beta1) applied to the nucleus tip.
	//beta2 is the vertical speed (0, 0, beta2) applied to the remainder of the nucleus.

	std::vector<double> V1 = {-0.0};/*, 0.0  ,  0.1966  ,  0.5547 ,  -0.4689 ,   0.2422 ,  -0.2229,
							   -0.4312 ,  -0.0185 ,   0.2887 ,   0.3187 ,   0.7140 ,  
								0.2231 ,  -0.1921 ,	  -0.5541 ,   -0.1542 ,   -0.1689 ,    0.4391 ,
							   -0.6661 ,  -0.6381 ,   0.6256 ,   0.0466 ,  -0.0610 ,   0.5134};
								*/
	std::vector<double> V2 = {0.0};/*, 0.0 ,  -0.4595 ,  -0.4129 ,   0.0954 ,   0.1764 ,   0.4186 ,
							  -0.5602 ,  -0.6082 ,  -0.5318 ,   0.3561 ,   0.0753 ,
							  -0.0917 ,  -0.2596 , 0.2871 ,  -0.3918 ,   0.5195 ,   0.5579 ,
							  -0.2805 ,   0.0133  , -0.0073 ,   0.7426 ,   0.0614 ,  -0.1506};
								*/
	std::vector<double> V3 = { 0.6390};/*, 0.0 ,  -0.5511 ,   0.0267 ,  -0.5240  , -0.4004 ,   0.2850 ,
							   0.2032 ,  -0.1771 ,   0.4048 ,   0.3461 ,  -0.2034 ,
							   0.5041 ,  -0.4535 ,	-0.1241 ,   0.5722 ,  -0.3748 ,  -0.1335 ,
							   -0.0851 ,   0.3213 ,   0.2389 ,   0.0044 ,  -0.7424 ,  -0.7450};
							   */
	//V1, V2, and V3 are the (x,y,z)-coordinate of the nucleus particles.

	for (int i = 0; i < V1.size(); i++){
		ljInfoVecs.LJ_PosX_all.push_back(V1[i]); 
		ljInfoVecs.LJ_PosY_all.push_back(V2[i]);
		ljInfoVecs.LJ_PosZ_all.push_back(V3[i]);
	}  
	
	double NUCLEUS_UPPERHEM_BASE = 0.5;
	double NUCLEUS_LOWERHEM_BASE = -0.6;
	//These values defines the z-coordinate requirement for nucleus particles to be considered tip-region or base-region. This is used to 
	// determine where to apply spring or constant force.

	//////////////////////////////////////////////////////////////////
	///////////////// END OF NUCLEUS RELATED /////////////////////////
	//////////////////////////////////////////////////////////////////

	/*std::vector<int> filament_base(generalParams.maxNodeCountLJ, -1); //= {0,1,2,3,4,5,6,7,8,9,10,11};//{35, 21, 38, etc if we need more points}
	double filament_strength = 0.0;
	double filament_strength_pull = 1.0*filament_strength;
	double filament_Rmin = ((max_height - min_height)/4.0);
	std::cout<<"filament_strength = "<<filament_strength<<std::endl;
	std::cout<<"filament_strength for vertical pull = "<<filament_strength_pull<<std::endl;
	std::cout<<"filament_Rmin = "<<filament_Rmin<<std::endl;
	
	//First, determine the initial membrane nodes having filament bridges
	//with the nuclei particles
	for (int i = 0; i < generalParams.maxNodeCountLJ; i++){
		if (i == 0){
			filament_base[i] = 35;
			continue;
		}
		for (int j = 0; j < generalParams.maxNodeCount; j++){
			double xsquared = (ljInfoVecs.LJ_PosX_all[i] - coordInfoVecs.nodeLocX[j])*
								(ljInfoVecs.LJ_PosX_all[i] - coordInfoVecs.nodeLocX[j]);
			double ysquared = (ljInfoVecs.LJ_PosY_all[i] - coordInfoVecs.nodeLocY[j])*
								(ljInfoVecs.LJ_PosY_all[i] - coordInfoVecs.nodeLocY[j]);
			double zsquared = (ljInfoVecs.LJ_PosZ_all[i] - coordInfoVecs.nodeLocZ[j])*
								(ljInfoVecs.LJ_PosZ_all[i] - coordInfoVecs.nodeLocZ[j]);
			double R = sqrt(xsquared + ysquared + zsquared);
			if (R < filament_Rmin*1.1 && j != 35){
				filament_base[i] = j;
				break;
			}
		}
	}*/
	
	//std::vector<double> filament_Rmin;
	//for (int i = 0; i < V3.size();i++){
	//	filament_Rmin.push_back(sqrt((V3[i] - coordInfoVecs.nodeLocZ[38])*(V3[i] - coordInfoVecs.nodeLocZ[38])));
	//}
	//double filament_Rmin = sqrt((V3.back() - coordInfoVecs.nodeLocZ[38])*(V3.back() - coordInfoVecs.nodeLocZ[38]));
	//This part calculates the filament connecting the minimum point (in terms of z-coordinate) to the base of the nuclei cluster.


	//////////////////////////////////////////////////////////////////
	/////////// IDENTIFYING REGIONS WITH DIFFERENT MECH PROP /////////
	//////////////////////////////////////////////////////////////////

	/*ljInfoVecs.forceX_all.reserve(ljInfoVecs.LJ_PosX_all.size());
	ljInfoVecs.forceY_all.reserve(ljInfoVecs.LJ_PosX_all.size());
	ljInfoVecs.forceZ_all.reserve(ljInfoVecs.LJ_PosX_all.size());

	generalParams.maxNodeCountLJ = ljInfoVecs.LJ_PosX_all.size();
	std::vector<int> nucleus_in_upperhem(generalParams.maxNodeCountLJ, -1);
	std::vector<int> nucleus_in_lowerhem(generalParams.maxNodeCountLJ, -1);
	for (int i = 0; i < generalParams.maxNodeCountLJ; i++){
		if (ljInfoVecs.LJ_PosZ_all[i] > NUCLEUS_UPPERHEM_BASE){
			nucleus_in_upperhem[i] = 1;
		}
		if (ljInfoVecs.LJ_PosZ_all[i] < NUCLEUS_LOWERHEM_BASE){
			nucleus_in_lowerhem[i] = 1;
		}
	}*/
	

	std::vector<int> out;
	//int ALPHA;

	std::vector<bool> boundary_edges;
	boundary_edges.reserve(coordInfoVecs.num_edges);
	for (int i = 0; i < coordInfoVecs.num_edges; i++){
		if (coordInfoVecs.edges2Triangles_1[i] == coordInfoVecs.edges2Triangles_2[i]){
			boundary_edges.push_back(true);
		}
		else {
			boundary_edges.push_back(false);
		}
	}

	std::vector<int> edgeIndices;
	edgeIndices.reserve(coordInfoVecs.num_edges);
	for (int i = 0; i < coordInfoVecs.num_edges; ++i){
		//edgeIndices.push_back(edge_to_ljparticle[i]);
		if (boundary_edges[i] == false){
			edgeIndices.push_back(i);
		}
		else {
			edgeIndices.push_back(-1);
		}
	}



	auto it = remove_if(edgeIndices.begin(), edgeIndices.end(),  [](const int i) {return i < 0; });
	edgeIndices.erase(it, edgeIndices.end());
	
	//std::vector<int> row2 = {35 ,   76 ,   79 ,  111 ,  113 ,  151 ,  153 ,  360 ,  361 ,  362 ,  363 ,  364 ,  365 ,  505 ,  506 ,  515 ,  516 ,  593 ,  632};
	//std::vector<int> nodes_to_center;
	generalParams.nodes_in_upperhem.resize(generalParams.maxNodeCount,-1);

	for (int i = 0; i < generalParams.maxNodeCount; i++){
		if (coordInfoVecs.nodeLocZ[i] > (generalParams.centerZ + weakened)){
			generalParams.nodes_in_upperhem[i] = 1;
		}
		else{
			generalParams.nodes_in_upperhem[i] = -1;
		}
	//	std::cout<<"nodes "<<i<<" "<<generalParams.nodes_in_upperhem[i]<<std::endl;		
	}

	//std::vector<int> nodes_to_center;
	//std::vector<int> nodes_in_tip;
	//nodes_in_tip.resize(generalParams.maxNodeCount);
	//for (int i = 0; i < generalParams.maxNodeCount; i++){
	//	if (coordInfoVecs.nodeLocZ[i] > (generalParams.centerZ + tip_base)){
	//		nodes_in_tip[i] = 1;
	//	}
	//	else{
	//		nodes_in_tip[i] = -1;
	//	}
	//	std::cout<<"nodes "<<i<<" "<<generalParams.nodes_in_upperhem[i]<<std::endl;		
	//}

	generalParams.triangles_in_upperhem.resize(coordInfoVecs.num_triangles);
	for (int i = 0; i < coordInfoVecs.num_triangles; i++){
		int aaa = generalParams.nodes_in_upperhem[coordInfoVecs.triangles2Nodes_1[i]];
		//std::cout<<aaa<<std::endl;
		int bbb = generalParams.nodes_in_upperhem[coordInfoVecs.triangles2Nodes_2[i]];
		//std::cout<<bbb<<std::endl;
		int ccc = generalParams.nodes_in_upperhem[coordInfoVecs.triangles2Nodes_3[i]];
		//std::cout<<ccc<<std::endl;
		if ((aaa+bbb+ccc)==3){
			generalParams.triangles_in_upperhem[i] = 1;
			//triangles_in_upperhem.push_back(i);
		}
		else if ((aaa+bbb+ccc)==1){
			generalParams.triangles_in_upperhem[i] = 0;
			//triangles_in_upperhem.push_back(i);
		}
		else{
			generalParams.triangles_in_upperhem[i] = -1;
		}
	//	std::cout<<"triangle "<<i<<" "<<generalParams.triangles_in_upperhem[i]<<std::endl;		
	}

	std::vector<int> edges_in_upperhem;
	generalParams.edges_in_upperhem.resize(coordInfoVecs.num_edges);
	for (int i = 0; i < coordInfoVecs.num_edges; i++){
		int aaa = generalParams.triangles_in_upperhem[coordInfoVecs.edges2Triangles_1[i]];//generalParams.nodes_in_upperhem[coordInfoVecs.edges2Nodes_1[i]];
		int bbb = generalParams.triangles_in_upperhem[coordInfoVecs.edges2Triangles_2[i]];//generalParams.nodes_in_upperhem[coordInfoVecs.edges2Nodes_2[i]];
		if (aaa == 1 && bbb == 1){
			generalParams.edges_in_upperhem[i] = 1;
			edges_in_upperhem.push_back(i);
		}
		else if (aaa == 1 || bbb == 1){
			generalParams.edges_in_upperhem[i] = 0;
		}
		else{
			generalParams.edges_in_upperhem[i] = -1;
		}
		
	}
	

	//Find the boundary of the nodes_in_upperhem region
	generalParams.boundaries_in_upperhem.resize(coordInfoVecs.num_edges);
	for (int i = 0; i < coordInfoVecs.num_edges; i++){
		double T1 = coordInfoVecs.edges2Triangles_1[i];
		double T2 = coordInfoVecs.edges2Triangles_2[i];
		if (generalParams.triangles_in_upperhem[T1] == 1 && generalParams.triangles_in_upperhem[T2] == 0){
			generalParams.boundaries_in_upperhem[i] = 1;
			//double bdry_node1 = coordInfoVecs.edges2Nodes_1[i];
			//double bdry_node2 = coordInfoVecs.edges2Nodes_2[i];
			//coordInfoVecs.isNodeFixed[bdry_node1] = true;
			//coordInfoVecs.isNodeFixed[bdry_node2] = true;
		}
		else if (generalParams.triangles_in_upperhem[T1] == 0 && generalParams.triangles_in_upperhem[T2] == 1){
			generalParams.boundaries_in_upperhem[i] = 1;
			//double bdry_node1 = coordInfoVecs.edges2Nodes_1[i];
			//double bdry_node2 = coordInfoVecs.edges2Nodes_2[i];
			//coordInfoVecs.isNodeFixed[bdry_node1] = true;
			//coordInfoVecs.isNodeFixed[bdry_node2] = true;
		}
		else {
			generalParams.boundaries_in_upperhem[i] = -1;
		}
	}
	generalParams.eq_total_boundary_length = generalParams.boundaries_in_upperhem.size()*generalParams.Rmin;
	
	

	int true_num_edges_in_upperhem = 0;
	for (int i = 0; i < edges_in_upperhem.size(); i++){
		if (edges_in_upperhem[i] != INT_MAX && edges_in_upperhem[i] >= 0){
		true_num_edges_in_upperhem += 1;
		}
	}
	

	//std::vector<int> edge_to_ljparticle;
	//generalParams.edge_to_ljparticle.reserve(coordInfoVecs.num_edges);
	for (int i = 0; i < coordInfoVecs.num_edges; i++){
		generalParams.edge_to_ljparticle.push_back(-1);
	};
	/////////////////////////////////////////////////////////////////////
	////////////// END OF IDENTIFYING REG. WITH DIFF. MECH PROP /////////
	/////////////////////////////////////////////////////////////////////

	
	ComputeVolume(
		generalParams,
		coordInfoVecs,
		linearSpringInfoVecs,
		ljInfoVecs
	);
	double initial_volume;
	initial_volume = generalParams.true_current_total_volume;
	generalParams.eq_total_volume = generalParams.true_current_total_volume*VOLUME_FACTOR;//This is for setting different equilibrium volume to mimic growth or shirnkage.
	std::cout<<"true_current_total_volume = "<<generalParams.true_current_total_volume<<std::endl;
	std::cout<<"eq_total_volume = "<<generalParams.eq_total_volume<<std::endl;

	//////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////// START OF ACTUAL SIMULATION /////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////
	
	while (runSim == true){
		//WHEN += 1;
		double current_time = 0.0;

		//generalParams.kT = 1.0;//reset kT before simulations starts.
		//Max_Runtime = 0.0;//2.5;
		int translate_counter = 0;
			while (current_time < (Max_Runtime)){
					translate_counter += 1;
					Solve_Forces();
				
					double energy_rep =
					ComputeMemRepulsionEnergy(
						coordInfoVecs,
						linearSpringInfoVecs, 
						capsidInfoVecs,
						generalParams,
						auxVecs);

					//now forces are computed, move nodes.
					
					

					/*for (int i = 0; i < ljInfoVecs.LJ_PosX_all.size(); i++){

						ljInfoVecs.LJ_PosX = ljInfoVecs.LJ_PosX_all[i];
				//		std::cout<<"LJ_PosX = "<<ljInfoVecs.LJ_PosX<<std::endl;
						ljInfoVecs.LJ_PosY = ljInfoVecs.LJ_PosY_all[i];
				//		std::cout<<"LJ_PosY = "<<ljInfoVecs.LJ_PosY<<std::endl;
						ljInfoVecs.LJ_PosZ = ljInfoVecs.LJ_PosZ_all[i];
						
						
					
						ComputeLJSprings(
							coordInfoVecs,
							ljInfoVecs,
							generalParams);
						ljInfoVecs.forceX_all[i] =  ljInfoVecs.forceX;
						ljInfoVecs.forceY_all[i] =  ljInfoVecs.forceY;
						ljInfoVecs.forceZ_all[i] =  ljInfoVecs.forceZ;						

						ComputeLJSprings_LJ(
							coordInfoVecs,
							ljInfoVecs,
							generalParams);
						ljInfoVecs.forceX_all[i] +=  ljInfoVecs.forceX;
						ljInfoVecs.forceY_all[i] +=  ljInfoVecs.forceY;
						ljInfoVecs.forceZ_all[i] +=  ljInfoVecs.forceZ;	
						
						for (int w = 0; w < nodes_in_tip.size(); w++){

							if (nodes_in_tip[i] == 1 && nucleus_in_upperhem[i] == 1){
								double R = sqrt( 
									(ljInfoVecs.LJ_PosX - coordInfoVecs.nodeLocX[w]) * (ljInfoVecs.LJ_PosX - coordInfoVecs.nodeLocX[w]) + 
									(ljInfoVecs.LJ_PosY - coordInfoVecs.nodeLocY[w]) * (ljInfoVecs.LJ_PosY - coordInfoVecs.nodeLocY[w]) + 
									(ljInfoVecs.LJ_PosZ - coordInfoVecs.nodeLocZ[w]) * (ljInfoVecs.LJ_PosZ - coordInfoVecs.nodeLocZ[w]) );
								double magnitude = -2.0*(pull_strength/2.0)*(R - ljInfoVecs.Rmin_M)*(1.0/R);
													//2*40.0*(1-exp(-1.5*(R-ljInfoVecs.Rmin_M)))*
													//(-exp(-1.5*(R-ljInfoVecs.Rmin_M)))*
													//(1.0/R);
								double forceX = -magnitude*(ljInfoVecs.LJ_PosX - coordInfoVecs.nodeLocX[w]);//xLoc_LR;
								double forceY = -magnitude*(ljInfoVecs.LJ_PosY - coordInfoVecs.nodeLocY[w]);//yLoc_LR;
								double forceZ = -magnitude*(ljInfoVecs.LJ_PosZ - coordInfoVecs.nodeLocZ[w]);//zLoc_LR;
								ljInfoVecs.forceX_all[i] +=  -forceX;
								ljInfoVecs.forceY_all[i] +=  -forceY;
								ljInfoVecs.forceZ_all[i] +=  -forceZ;
								//coordInfoVecs.nodeForceX[35] += forceX;
								//coordInfoVecs.nodeForceY[35] += forceY;
								//coordInfoVecs.nodeForceZ[35] += forceZ;
							}
						}
						//for (int i = 0; i < filament_base.size(); i++){
						if (filament_base[i] != -1){
							double R = sqrt( 
								(ljInfoVecs.LJ_PosX - coordInfoVecs.nodeLocX[filament_base[i]]) * (ljInfoVecs.LJ_PosX - coordInfoVecs.nodeLocX[filament_base[i]]) + 
								(ljInfoVecs.LJ_PosY - coordInfoVecs.nodeLocY[filament_base[i]]) * (ljInfoVecs.LJ_PosY - coordInfoVecs.nodeLocY[filament_base[i]]) + 
								(ljInfoVecs.LJ_PosZ - coordInfoVecs.nodeLocZ[filament_base[i]]) * (ljInfoVecs.LJ_PosZ - coordInfoVecs.nodeLocZ[filament_base[i]]) );
								double magnitude;
								if (i == 0){
									magnitude = -2.0*(filament_strength_pull/2.0)*(R - filament_Rmin)*(1.0/R);
								}
								else{
									magnitude = -2.0*(filament_strength/2.0)*(R - filament_Rmin)*(1.0/R);
								}
							double forceX = -magnitude*(ljInfoVecs.LJ_PosX - coordInfoVecs.nodeLocX[filament_base[i]]);  
							double forceY = -magnitude*(ljInfoVecs.LJ_PosY - coordInfoVecs.nodeLocY[filament_base[i]]);  
							double forceZ = -magnitude*(ljInfoVecs.LJ_PosZ - coordInfoVecs.nodeLocZ[filament_base[i]]);  
							ljInfoVecs.forceX_all[i] +=  -forceX;
							ljInfoVecs.forceY_all[i] +=  -forceY;
							ljInfoVecs.forceZ_all[i] +=  -forceZ;
							coordInfoVecs.nodeForceX[filament_base[i]] += forceX;
							coordInfoVecs.nodeForceY[filament_base[i]] += forceY;
							coordInfoVecs.nodeForceZ[filament_base[i]] += forceZ;	
						//}	
							}
						
					}*/

					double beta;
					/*bool target_height_reached = false;
					if (ljInfoVecs.LJ_PosZ_all[0] >= (min_height + 0.875*(max_height - min_height))){
						target_height_reached = true;
					}*/

					//std::random_device rand_dev;
					//std::mt19937 generator_noise(rand_dev());
					//std::normal_distribution<double> distribution_noise(0.0, 0.01);

					/*for (int i = 0; i < ljInfoVecs.LJ_PosX_all.size(); i++){
						beta = 0.0;
						if(nucleus_in_upperhem[i] == 1){
							beta = beta1;
						}
						else{
							beta = beta2;
						} 

						if (target_height_reached == true){
							beta = 0.0;
						}

						//std::random_device rand_dev;
						//std::mt19937 generator_noise(rand_dev());
						//std::normal_distribution<double> distribution_noise(0.0, 0.01);
						double noise1 = 0.0;//distribution_noise(generator_noise);
						double noise2 = 0.0;//distribution_noise(generator_noise);
						double noise3 = 0.0;//distribution_noise(generator_noise);

						ljInfoVecs.LJ_PosX_all[i] = ljInfoVecs.LJ_PosX_all[i] + generalParams.dt * (ljInfoVecs.forceX_all[i] + noise1);
						ljInfoVecs.LJ_PosY_all[i] = ljInfoVecs.LJ_PosY_all[i] + generalParams.dt * (ljInfoVecs.forceY_all[i] + noise2);
						ljInfoVecs.LJ_PosZ_all[i] = ljInfoVecs.LJ_PosZ_all[i] + generalParams.dt * (ljInfoVecs.forceZ_all[i] + beta + noise3);
					
					}*/
				
				/*for (int i = 0; i < generalParams.maxNodeCount; i++){
					if (coordInfoVecs.isNodeFixed[i] == true){
						coordInfoVecs.nodeForceX[i] = 0.0;
						coordInfoVecs.nodeForceY[i] = 0.0;
						coordInfoVecs.nodeForceZ[i] = 0.0;
					}
				}*/

				AdvancePositions(
					coordInfoVecs,
					generalParams,
					domainParams);

				if (translate_counter % translate_frequency == 1){

					newcenterX = 0.0;
					newcenterY = 0.0;
					newcenterZ = 0.0;
					for (int i = 0; i < generalParams.maxNodeCount; i++){//for (int i = 0; i < coordInfoVecs.nodeLocX.size(); i++){
						newcenterX += coordInfoVecs.nodeLocX[i];
						newcenterY += coordInfoVecs.nodeLocY[i];
						newcenterZ += coordInfoVecs.nodeLocZ[i];
					}
					newcenterX = newcenterX/generalParams.maxNodeCount; //coordInfoVecs.nodeLocX.size();
					newcenterY = newcenterY/generalParams.maxNodeCount; //coordInfoVecs.nodeLocX.size();
					newcenterZ = newcenterZ/generalParams.maxNodeCount; //coordInfoVecs.nodeLocX.size();
					displacementX = newcenterX - generalParams.centerX;
					displacementY = newcenterY - generalParams.centerY;
					displacementZ = newcenterZ - generalParams.centerZ;

					for (int i = 0; i < generalParams.maxNodeCount; i++){
						coordInfoVecs.nodeLocX[i] += -displacementX;
						coordInfoVecs.nodeLocY[i] += -displacementY;
						coordInfoVecs.nodeLocZ[i] += -displacementZ;
					}
					for (int i = 0; i < ljInfoVecs.LJ_PosX_all.size(); i++){
						ljInfoVecs.LJ_PosX_all[i] += -displacementX;
						ljInfoVecs.LJ_PosY_all[i] += -displacementY;
						ljInfoVecs.LJ_PosZ_all[i] += -displacementZ;
					}

					//Here we re-establish the new filament base according to the current location of nuclei nodes
					/* int maxElementIndex = std::max_element(coordInfoVecs.nodeLocZ.begin(),coordInfoVecs.nodeLocZ.end()) - coordInfoVecs.nodeLocZ.begin();
					for (int i = 0; i < generalParams.maxNodeCountLJ; i++){
						if (i == 0){
							filament_base[i] = maxElementIndex;
							continue;
						}
						for (int j = 0; j < generalParams.maxNodeCount; j++){
							double xsquared = (ljInfoVecs.LJ_PosX_all[i] - coordInfoVecs.nodeLocX[j])*
												(ljInfoVecs.LJ_PosX_all[i] - coordInfoVecs.nodeLocX[j]);
							double ysquared = (ljInfoVecs.LJ_PosY_all[i] - coordInfoVecs.nodeLocY[j])*
												(ljInfoVecs.LJ_PosY_all[i] - coordInfoVecs.nodeLocY[j]);
							double zsquared = (ljInfoVecs.LJ_PosZ_all[i] - coordInfoVecs.nodeLocZ[j])*
												(ljInfoVecs.LJ_PosZ_all[i] - coordInfoVecs.nodeLocZ[j]);
							double R = sqrt(xsquared + ysquared + zsquared);
							if (R < (max_height - min_height)/2.0 && j != maxElementIndex){
								filament_base[i] = j;
								break;
							}
							else{filament_base[i] = -1;}
						}
					} */
				}
							
					new_total_energy = linearSpringInfoVecs.linear_spring_energy + 
						areaTriangleInfoVecs.area_triangle_energy + 
						bendingTriangleInfoVecs.bending_triangle_energy + 
						0.5*energy_rep + 
						//ljInfoVecs.lj_energy_M +
						//ljInfoVecs.lj_energy_LJ +
						generalParams.volume_energy;
//std::cout<<"LINEAR ENERGY = "<<linearSpringInfoVecs.linear_spring_energy<<std::endl;
//std::cout<<"BEND ENERGY = "<<bendingTriangleInfoVecs.bending_triangle_energy<<std::endl;
//std::cout<<"REPULSION ENERGY = "<<energy_rep<<std::endl;
//std::cout<<"VOLUME ENERGY = "<<generalParams.volume_energy<<std::endl;
				energy_gradient = sqrt((new_total_energy - old_total_energy)*(new_total_energy - old_total_energy));
				old_total_energy = new_total_energy;
				current_time+=generalParams.dt;
				

			}
			
		   
			/*max_height = -10000.0;
			min_height = 10000.0;
			for (int k = 0; k < generalParams.maxNodeCount; k++){
				if (coordInfoVecs. nodeLocZ[k] >= max_height){
					max_height = coordInfoVecs. nodeLocZ[k];
				}
				if (coordInfoVecs.nodeLocZ[k] <= min_height){
					min_height = coordInfoVecs.nodeLocZ[k];
				}
			}*/

		std::cout<<"current time (1st iter before edgeswap): "<< current_time << std::endl;
		std::cout<<"current total energy (1st iter before edgeswap) = "<<new_total_energy<<std::endl;
		std::cout<<"true_current_total_volume = "<<generalParams.true_current_total_volume<<std::endl;
		std::cout<<"eq_total_volume = "<<generalParams.eq_total_volume<<std::endl;
		std::cout<<"current KBT = "<<generalParams.kT<<std::endl;
		if (isnan(new_total_energy)==1){
			std::cout<<"Nan or Inf position update !!!!"<<std::endl;
			runSim = false;
			break;
		}
	
		//edgeswap_ptr->transferDtoH(coordInfoVecs, build_ptr->hostSetInfoVecs);//Currently this is treated as a backup of coordInfoVecs
		////storage->print_VTK_File();
		////storage->storeVariables();
		//runSim = false;
		//break;

		int edgeswap_iteration = 0;
		//double preswap_energy = new_total_energy;
		//double postswap_energy;
		//double Ediff = 0.0;
		//initial_kT = generalParams.kT;
		num_edge_loop = round(true_num_edges_in_upperhem*SAMPLE_SIZE);
		if (num_edge_loop == 0){
			num_edge_loop = 1;
		}	
		std::cout<<"num_edge_loop = "<<num_edge_loop<<std::endl;
	
 		while (initial_kT > 0){
 					////////////////////NOW RELAX THE ATTEMPTED EDGESWAP//////////////////////
					 current_time = 0.0;
					 translate_counter = 0;
					 double VOLUME_RATIO = generalParams.true_current_total_volume/generalParams.eq_total_volume;
					if (VOLUME_RATIO > 0.75 && VOLUME_FACTOR <= 4.0){
						VOLUME_FACTOR += 0.1;
						generalParams.eq_total_volume = initial_volume*VOLUME_FACTOR;
					};

 					while (current_time < Max_Runtime){
						 translate_counter += 1;
						 //std::cout<<"ERROR BEFORE RELAXATION"<<std::endl;
						 Solve_Forces();

 						double energy_rep =
 						ComputeMemRepulsionEnergy(
 							coordInfoVecs,
 							linearSpringInfoVecs, 
 							capsidInfoVecs,
 							generalParams,
							 auxVecs);
					if ((generalParams.true_current_total_volume/initial_volume) > 2.0){
runSim = false;
initial_kT = -0.01;
break;
}
					
				
 						/*for (int i = 0; i < ljInfoVecs.LJ_PosX_all.size(); i++){
 							ljInfoVecs.LJ_PosX = ljInfoVecs.LJ_PosX_all[i];
 							ljInfoVecs.LJ_PosY = ljInfoVecs.LJ_PosY_all[i];
							 ljInfoVecs.LJ_PosZ = ljInfoVecs.LJ_PosZ_all[i];
							 
							
					
 							ComputeLJSprings(
 								coordInfoVecs,
 								ljInfoVecs,
 								generalParams);
 							ljInfoVecs.forceX_all[i] =  ljInfoVecs.forceX;
 							ljInfoVecs.forceY_all[i] =  ljInfoVecs.forceY;
							 ljInfoVecs.forceZ_all[i] =  ljInfoVecs.forceZ;

 							ComputeLJSprings_LJ(
 								coordInfoVecs,
 								ljInfoVecs,
 								generalParams);
 							ljInfoVecs.forceX_all[i] +=  ljInfoVecs.forceX;
 							ljInfoVecs.forceY_all[i] +=  ljInfoVecs.forceY;
							 ljInfoVecs.forceZ_all[i] +=  ljInfoVecs.forceZ;
							 
							//  for (int w = 0; w < nodes_in_tip.size(); w++){
							// 	if (nodes_in_tip[w] == 1 && nucleus_in_upperhem[i] == 1){
							// 		double R = sqrt( 
							// 			(ljInfoVecs.LJ_PosX - coordInfoVecs.nodeLocX[w]) * (ljInfoVecs.LJ_PosX - coordInfoVecs.nodeLocX[w]) + 
							// 			(ljInfoVecs.LJ_PosY - coordInfoVecs.nodeLocY[w]) * (ljInfoVecs.LJ_PosY - coordInfoVecs.nodeLocY[w]) + 
							// 			(ljInfoVecs.LJ_PosZ - coordInfoVecs.nodeLocZ[w]) * (ljInfoVecs.LJ_PosZ - coordInfoVecs.nodeLocZ[w]) );
							// 		double magnitude = -2.0*(pull_strength/2.0)*(R - ljInfoVecs.Rmin_M)*(1.0/R);
							// 							//2*40.0*(1-exp(-1.5*(R-ljInfoVecs.Rmin_M)))*
							// 							//(-exp(-1.5*(R-ljInfoVecs.Rmin_M)))*
							// 							//(1.0/R);
							// 		double forceX = -magnitude*(ljInfoVecs.LJ_PosX - coordInfoVecs.nodeLocX[w]);//xLoc_LR;
							// 		double forceY = -magnitude*(ljInfoVecs.LJ_PosY - coordInfoVecs.nodeLocY[w]);//yLoc_LR;
							// 		double forceZ = -magnitude*(ljInfoVecs.LJ_PosZ - coordInfoVecs.nodeLocZ[w]);//zLoc_LR;
							// 		ljInfoVecs.forceX_all[i] +=  -forceX;
							// 		ljInfoVecs.forceY_all[i] +=  -forceY;
							// 		ljInfoVecs.forceZ_all[i] +=  -forceZ;
							// 		//coordInfoVecs.nodeForceX[35] += forceX;
							// 		//coordInfoVecs.nodeForceY[35] += forceY;
							// 		//coordInfoVecs.nodeForceZ[35] += forceZ;
							// 	}
							// }
							//for (int k = 0; k < filament_base.size(); k++){
							if (filament_base[i] != -1){
								double R = sqrt( 
									(ljInfoVecs.LJ_PosX - coordInfoVecs.nodeLocX[filament_base[i]]) * (ljInfoVecs.LJ_PosX - coordInfoVecs.nodeLocX[filament_base[i]]) + 
									(ljInfoVecs.LJ_PosY - coordInfoVecs.nodeLocY[filament_base[i]]) * (ljInfoVecs.LJ_PosY - coordInfoVecs.nodeLocY[filament_base[i]]) + 
									(ljInfoVecs.LJ_PosZ - coordInfoVecs.nodeLocZ[filament_base[i]]) * (ljInfoVecs.LJ_PosZ - coordInfoVecs.nodeLocZ[filament_base[i]]) );
								//if (R > (filament_Rmin + 0.25)){
								//	continue;
								//}
								double magnitude;
								if (i == 0){
									magnitude = -2.0*(filament_strength_pull/2.0)*(R - filament_Rmin)*(1.0/R);
								}
								else{
									magnitude = -2.0*(filament_strength/2.0)*(R - filament_Rmin)*(1.0/R);
								}
								double forceX = -magnitude*(ljInfoVecs.LJ_PosX - coordInfoVecs.nodeLocX[filament_base[i]]);  
								double forceY = -magnitude*(ljInfoVecs.LJ_PosY - coordInfoVecs.nodeLocY[filament_base[i]]);  
								double forceZ = -magnitude*(ljInfoVecs.LJ_PosZ - coordInfoVecs.nodeLocZ[filament_base[i]]);  
								ljInfoVecs.forceX_all[i] +=  -forceX;
								ljInfoVecs.forceY_all[i] +=  -forceY;
								ljInfoVecs.forceZ_all[i] +=  -forceZ;
								coordInfoVecs.nodeForceX[filament_base[i]] += forceX;
								coordInfoVecs.nodeForceY[filament_base[i]] += forceY;
								coordInfoVecs.nodeForceZ[filament_base[i]] += forceZ;	
							//}	
								}
							 
 						}*/
					
 						//now forces are computed, move nodes.
						 double beta;
						/* bool target_height_reached = false;
						if (ljInfoVecs.LJ_PosZ_all[0] >= (min_height + 0.75*(max_height - min_height))){
							target_height_reached = true;
						}*/

						//std::random_device rand_dev;
						//std::mt19937 generator_noise(rand_dev());
						//std::normal_distribution<double> distribution_noise(0.0, 0.01);
						/*for (int i = 0; i < ljInfoVecs.LJ_PosX_all.size(); i++){
							beta = 0.0;
							if(nucleus_in_upperhem[i] == 1){
							   beta = beta1;
						   }
						   else{
							   beta = beta2;
						   } 

						   if (target_height_reached == true){
							   beta = 0.0;
						   }

						   if (nucleus_in_lowerhem[i] == 1){ //(i == (generalParams.maxNodeCountLJ-1)){
							double R = sqrt( 
								(ljInfoVecs.LJ_PosX_all[i] - coordInfoVecs.nodeLocX[38]) * (ljInfoVecs.LJ_PosX_all[i] - coordInfoVecs.nodeLocX[38]) + 
								(ljInfoVecs.LJ_PosY_all[i] - coordInfoVecs.nodeLocY[38]) * (ljInfoVecs.LJ_PosY_all[i] - coordInfoVecs.nodeLocY[38]) + 
								(ljInfoVecs.LJ_PosZ_all[i] - coordInfoVecs.nodeLocZ[38]) * (ljInfoVecs.LJ_PosZ_all[i] - coordInfoVecs.nodeLocZ[38]) );
							double magnitude = -2.0*(filament_strength/2.0)*(R - filament_Rmin[i])*(1.0/R);
							double forceX = -magnitude*(ljInfoVecs.LJ_PosX_all[i] - coordInfoVecs.nodeLocX[38]);  
							double forceY = -magnitude*(ljInfoVecs.LJ_PosY_all[i] - coordInfoVecs.nodeLocY[38]);  
							double forceZ = -magnitude*(ljInfoVecs.LJ_PosZ_all[i] - coordInfoVecs.nodeLocZ[38]);  
							ljInfoVecs.forceX_all[i] +=  -forceX;
							ljInfoVecs.forceY_all[i] +=  -forceY;
							ljInfoVecs.forceZ_all[i] +=  -forceZ;	
							
						}
						   //std::random_device rand_dev;
						   //std::mt19937 generator_noise(rand_dev());
						   //std::normal_distribution<double> distribution_noise(0.0, 0.01);
						   double noise1 = 0.0;//distribution_noise(generator_noise);
						   double noise2 = 0.0;//distribution_noise(generator_noise);
						   double noise3 = 0.0;//distribution_noise(generator_noise);
   
						   ljInfoVecs.LJ_PosX_all[i] = ljInfoVecs.LJ_PosX_all[i] + generalParams.dt * (ljInfoVecs.forceX_all[i] + noise1);
						   ljInfoVecs.LJ_PosY_all[i] = ljInfoVecs.LJ_PosY_all[i] + generalParams.dt * (ljInfoVecs.forceY_all[i] + noise2);
						   ljInfoVecs.LJ_PosZ_all[i] = ljInfoVecs.LJ_PosZ_all[i] + generalParams.dt * (ljInfoVecs.forceZ_all[i] + beta + noise3);
					   
					   }*/

					   /*for (int i = 0; i < generalParams.maxNodeCount; i++){
						if (coordInfoVecs.isNodeFixed[i] == true){
							coordInfoVecs.nodeForceX[i] = 0.0;
							coordInfoVecs.nodeForceY[i] = 0.0;
							coordInfoVecs.nodeForceZ[i] = 0.0;
						}
					}*/
						 
 						AdvancePositions(
 							coordInfoVecs,
 							generalParams,
							 domainParams);
						
						if (translate_counter % translate_frequency == 1){
							newcenterX = 0.0;
							newcenterY = 0.0;
							newcenterZ = 0.0;
							for (int i = 0; i < generalParams.maxNodeCount; i++){//for (int i = 0; i < coordInfoVecs.nodeLocX.size(); i++){
								newcenterX += coordInfoVecs.nodeLocX[i];
								newcenterY += coordInfoVecs.nodeLocY[i];
								newcenterZ += coordInfoVecs.nodeLocZ[i];
							}
							newcenterX = newcenterX/generalParams.maxNodeCount; //coordInfoVecs.nodeLocX.size();
							newcenterY = newcenterY/generalParams.maxNodeCount; //coordInfoVecs.nodeLocX.size();
							newcenterZ = newcenterZ/generalParams.maxNodeCount; //coordInfoVecs.nodeLocX.size();
							displacementX = newcenterX - generalParams.centerX;
							displacementY = newcenterY - generalParams.centerY;
							displacementZ = newcenterZ - generalParams.centerZ;
			
							for (int i = 0; i < generalParams.maxNodeCount; i++){
							coordInfoVecs.nodeLocX[i] += -displacementX;
							coordInfoVecs.nodeLocY[i] += -displacementY;
							coordInfoVecs.nodeLocZ[i] += -displacementZ;
							}
							for (int i = 0; i < ljInfoVecs.LJ_PosX_all.size(); i++){
								ljInfoVecs.LJ_PosX_all[i] += -displacementX;
								ljInfoVecs.LJ_PosY_all[i] += -displacementY;
								ljInfoVecs.LJ_PosZ_all[i] += -displacementZ;
							}

							/* int maxElementIndex = std::max_element(coordInfoVecs.nodeLocZ.begin(),coordInfoVecs.nodeLocZ.end()) - coordInfoVecs.nodeLocZ.begin();
							for (int i = 0; i < generalParams.maxNodeCountLJ; i++){
								if (i == 0){
									filament_base[i] = maxElementIndex;
									continue;
								}
								for (int j = 0; j < generalParams.maxNodeCount; j++){
									double xsquared = (ljInfoVecs.LJ_PosX_all[i] - coordInfoVecs.nodeLocX[j])*
														(ljInfoVecs.LJ_PosX_all[i] - coordInfoVecs.nodeLocX[j]);
									double ysquared = (ljInfoVecs.LJ_PosY_all[i] - coordInfoVecs.nodeLocY[j])*
														(ljInfoVecs.LJ_PosY_all[i] - coordInfoVecs.nodeLocY[j]);
									double zsquared = (ljInfoVecs.LJ_PosZ_all[i] - coordInfoVecs.nodeLocZ[j])*
														(ljInfoVecs.LJ_PosZ_all[i] - coordInfoVecs.nodeLocZ[j]);
									double R = sqrt(xsquared + ysquared + zsquared);
									if (R < (max_height - min_height)/2.0 && j != maxElementIndex){
										filament_base[i] = j;
										break;
									}
									else{filament_base[i] = -1;}
								}
							} */
						
							ComputeVolume(
								generalParams,
								coordInfoVecs,
								linearSpringInfoVecs,
								ljInfoVecs);
							//std::cout<<"ERROR 1"<<std::endl;
							
							edgeswap_ptr->transferDtoH(coordInfoVecs, build_ptr->hostSetInfoVecs);
							//std::cout<<"ERROR 1.5"<<std::endl;
							for (int edge_loop = 0; edge_loop < num_edge_loop; edge_loop++) {
								//std::cout<<"edge_loop = "<<edge_loop<<std::endl;
								
								std::random_device rand_dev;
								std::mt19937 generator(rand_dev());
							   
							   std::uniform_int_distribution<int> distribution(1,edges_in_upperhem.size());
							   
							   int dice_roll = distribution(generator);
							   
							   int edge = edges_in_upperhem[dice_roll - 1];
							   //int edge = dice_roll -1;
							   while (generalParams.boundaries_in_upperhem[edge] == 1 || edge == INT_MAX){
									dice_roll = distribution(generator);
									
									edge =  edges_in_upperhem[dice_roll - 1];
									//edge = dice_roll -1;
								 }

								int ALPHA = edgeswap_ptr->edge_swap_host_vecs(
									edge,
									generalParams,
									build_ptr->hostSetInfoVecs,
									linearSpringInfoVecs,
									bendingTriangleInfoVecs,
									areaTriangleInfoVecs);
								
							}
							//NOTE: EDGESWAP ALGORITHM CURRENTLY IS WRITTEN TO ALLOW AT MOST 8 NEIGHBORING NODES PER NODE.
							//std::cout<<"ERROR 2"<<std::endl;
							edgeswap_ptr->transferHtoD(coordInfoVecs, build_ptr->hostSetInfoVecs);//Currently this is treated as a backup of coordInfoVecs
							//std::cout<<"ERROR 2.5"<<std::endl;
							
							
						}
						
 						new_total_energy = linearSpringInfoVecs.linear_spring_energy + 
 							areaTriangleInfoVecs.area_triangle_energy + 
 							bendingTriangleInfoVecs.bending_triangle_energy +
 							0.5*energy_rep +
 							//ljInfoVecs.lj_energy_M +  
							// ljInfoVecs.lj_energy_LJ +
							 generalParams.volume_energy;
 						//std::cout<<"new_total_energy = "<<new_total_energy<<std::endl;

 						energy_gradient = sqrt((new_total_energy - old_total_energy)*(new_total_energy - old_total_energy));
 						
 					old_total_energy = new_total_energy;
 					current_time+=generalParams.dt;
					 }
					 
					 
					/*max_height = -10000.0;
					min_height = 10000.0;
					for (int k = 0; k < generalParams.maxNodeCount; k++){
						if (coordInfoVecs. nodeLocZ[k] >= max_height){
							max_height = coordInfoVecs. nodeLocZ[k];
						}
						if (coordInfoVecs.nodeLocZ[k] <= min_height){
							min_height = coordInfoVecs.nodeLocZ[k];
						}
					}*/
											
			
 								
 					if (edgeswap_iteration % RECORD_TIME == 0){
						generalParams.true_num_edges = 0;
						for (int i = 0; i < coordInfoVecs.num_edges; i++){
							if (coordInfoVecs.edges2Nodes_1[i] != INT_MAX && coordInfoVecs.edges2Nodes_2[i] != INT_MAX){
								generalParams.true_num_edges += 1;
							}
						 }
						 storage->print_VTK_File();
						 //storage->storeVariables();
						 std::cout<<"current total energy = "<< new_total_energy<<std::endl;
						 std::cout<<"true current total volume = "<<generalParams.true_current_total_volume<<std::endl;
						std::cout<<"equilibrium total volume = "<<generalParams.eq_total_volume<<std::endl;
						// std::cout<<"LINEAR ENERGY = "<<linearSpringInfoVecs.linear_spring_energy<<std::endl;
						//std::cout<<"BEND ENRGY = "<<bendingTriangleInfoVecs.bending_triangle_energy<<std::endl;
						//std::cout<<"REPULSION ENERGY = "<<energy_rep<<std::endl;
						//std::cout<<"VOLUME ENERGY = "<<generalParams.volume_energy<<std::endl;
 					}
 					if (edgeswap_iteration == NKBT-1 ){
 						//storage->storeVariables();
					 }

					

					 edgeswap_iteration += 1;
					 
					/*if (edgeswap_iteration % GROWTH_TIME == 0){

						for (int i = 0; i < coordInfoVecs.nodeLocX.size(); i++){
							generalParams.centerX += coordInfoVecs.nodeLocX[i];
							generalParams.centerY += coordInfoVecs.nodeLocY[i];
							generalParams.centerZ += coordInfoVecs.nodeLocZ[i];
						}
						generalParams.centerX = generalParams.centerX/coordInfoVecs.nodeLocX.size();
						generalParams.centerY = generalParams.centerY/coordInfoVecs.nodeLocX.size();
						generalParams.centerZ = generalParams.centerZ/coordInfoVecs.nodeLocX.size();

						double x,y,z;
						std::random_device rand_dev0;
						std::mt19937 generator0(rand_dev0());
						std::uniform_real_distribution<double> guess(generalParams.centerX-1.0, generalParams.centerX+1.0);
						x = 5.0;//guess(generator0);
						y = 5.0;//guess(generator0);
						z = 5.0;//guess(generator0);
						bool goodchoice = false;
						double GAP;
						while (sqrt(x*x + y*y + z*z) > (2.0) && goodchoice == false){
							x = guess(generator0);
							y = guess(generator0);
							z = guess(generator0);
							if (sqrt(x*x + y*y + z*z) > 2.0){
								continue;
							}
							for (int i = 0; i < ljInfoVecs.LJ_PosX_all.size(); i++){
								GAP = sqrt((x-ljInfoVecs.LJ_PosX_all[i])*(x-ljInfoVecs.LJ_PosX_all[i]) +
											(y-ljInfoVecs.LJ_PosY_all[i])*(x-ljInfoVecs.LJ_PosY_all[i]) +
											(z-ljInfoVecs.LJ_PosZ_all[i])*(x-ljInfoVecs.LJ_PosZ_all[i]));
								if (GAP < 0.65){
									goodchoice = false;
									break;
								}
								else{goodchoice = true;}
							}
						}
						ljInfoVecs.LJ_PosX_all.push_back(x);
						ljInfoVecs.LJ_PosY_all.push_back(y);
						ljInfoVecs.LJ_PosZ_all.push_back(z);
						ljInfoVecs.forceX_all.resize(ljInfoVecs.LJ_PosX_all.size());
						ljInfoVecs.forceY_all.resize(ljInfoVecs.LJ_PosX_all.size());
						ljInfoVecs.forceZ_all.resize(ljInfoVecs.LJ_PosX_all.size());
						generalParams.maxNodeCountLJ = ljInfoVecs.LJ_PosX_all.size();
					}*/
 					//std::cout<<"edgeswap_iteration = "<<edgeswap_iteration<<std::endl;
 					if (edgeswap_iteration == NKBT){
 						generalParams.kT = -1.0;//generalParams.kT - 0.072;
 						std::cout<<"Current kBT = "<<generalParams.kT<<std::endl;
 						edgeswap_iteration = 0;
 					}
 					if (generalParams.kT < min_kT){
 						initial_kT = -1.0;
					runSim = false;
					break;
					 }

//std::cout<<"ERROR BEFORE GROWTH"<<std::endl;
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////// GROWTH OF THE CELL (MEMBRANE) ////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// GROWTH ALGORITHM IS CURRENTLY WRITTEN TO ALLOW NO MORE THEN 8 NEIGHBORING NODES PER NODE /////////////////////////////

VectorShuffleForGrowthLoop.clear();
for (int i = 0; i < coordInfoVecs.num_edges; i++){
	VectorShuffleForGrowthLoop.push_back(i);
}
std::random_device rand_dev;
std::mt19937 generator2(rand_dev());
std::shuffle(std::begin(VectorShuffleForGrowthLoop), std::end(VectorShuffleForGrowthLoop), generator2);

bool triggered = false;
//int triggered_counter = 0;
for (int p = 0; p < VectorShuffleForGrowthLoop.size(); p++){
	//std::cout<<"p = "<<p<<std::endl;
		int k = VectorShuffleForGrowthLoop[p];
		if (coordInfoVecs. edges2Nodes_1[k] == INT_MAX || coordInfoVecs. edges2Nodes_2[k] == INT_MAX){
			continue;
		}
		if (generalParams.boundaries_in_upperhem[k] == 1){
			continue;
		}
		if (generalParams.edges_in_upperhem[k] < 0 || generalParams.edges_in_upperhem[k] == INT_MAX){
			continue;
		}
		int iedge = k;
		int elem1 = coordInfoVecs.edges2Triangles_1[iedge];
		int elem2 = coordInfoVecs.edges2Triangles_2[iedge];
        int first_v = coordInfoVecs.triangles2Nodes_1[elem1];
        int second_v = coordInfoVecs.triangles2Nodes_2[elem1];
        int third_v = coordInfoVecs.triangles2Nodes_3[elem1];
        double v1x = coordInfoVecs.nodeLocX[second_v] - coordInfoVecs.nodeLocX[first_v];
        double v1y = coordInfoVecs.nodeLocY[second_v] - coordInfoVecs.nodeLocY[first_v];
        double v1z = coordInfoVecs.nodeLocZ[second_v] - coordInfoVecs.nodeLocZ[first_v];
        double v2x = coordInfoVecs.nodeLocX[third_v] - coordInfoVecs.nodeLocX[first_v];
        double v2y = coordInfoVecs.nodeLocY[third_v] - coordInfoVecs.nodeLocY[first_v];
        double v2z = coordInfoVecs.nodeLocZ[third_v] - coordInfoVecs.nodeLocZ[first_v];
		//std::cout<<"v1x "<<v1x<<std::endl;
		//std::cout<<"v1y "<<v1y<<std::endl;
		//std::cout<<"v1z "<<v1z<<std::endl;
		//std::cout<<"v2x "<<v2x<<std::endl;
		//std::cout<<"v2y "<<v2y<<std::endl;
		//std::cout<<"v2z "<<v2z<<std::endl;
		double This_area_v = sqrt((v1y*v2z - v2y*v1z)*(v1y*v2z - v2y*v1z) + 
                                ((-v1x*v2z) + (v2x*v1z))*((-v1x*v2z) + (v2x*v1z)) +
								(v1x*v2y - v2x*v1y)*(v1x*v2y - v2x*v1y))/2.0;
		//						std::cout<<"This_area_v = "<<This_area_v<<std::endl;
		
		int first_w = coordInfoVecs.triangles2Nodes_1[elem2];
        int second_w = coordInfoVecs.triangles2Nodes_2[elem2];
        int third_w = coordInfoVecs.triangles2Nodes_3[elem2];
		double w1x = coordInfoVecs.nodeLocX[second_w] - coordInfoVecs.nodeLocX[first_w];
        double w1y = coordInfoVecs.nodeLocY[second_w] - coordInfoVecs.nodeLocY[first_w];
        double w1z = coordInfoVecs.nodeLocZ[second_w] - coordInfoVecs.nodeLocZ[first_w];
        double w2x = coordInfoVecs.nodeLocX[third_w] - coordInfoVecs.nodeLocX[first_w];
        double w2y = coordInfoVecs.nodeLocY[third_w] - coordInfoVecs.nodeLocY[first_w];
		double w2z = coordInfoVecs.nodeLocZ[third_w] - coordInfoVecs.nodeLocZ[first_w];
		//std::cout<<"w1x "<<w1x<<std::endl;
		//std::cout<<"w1y "<<w1y<<std::endl;
		//std::cout<<"w1z "<<w1z<<std::endl;
		//std::cout<<"w2x "<<w2x<<std::endl;
		//std::cout<<"w2y "<<w2y<<std::endl;
		//std::cout<<"w2z "<<w2z<<std::endl;
        double This_area_w = sqrt((w1y*w2z - w2y*w1z)*(w1y*w2z - w2y*w1z) + 
                                ((-w1x*w2z) + (w2x*w1z))*((-w1x*w2z) + (w2x*w1z)) +
								(w1x*w2y - w2x*w1y)*(w1x*w2y - w2x*w1y))/2.0;
		//						std::cout<<"This_area_w = "<<This_area_w<<std::endl;
								
		int node1 = coordInfoVecs.edges2Nodes_1[iedge];
		int node2 = coordInfoVecs.edges2Nodes_2[iedge];
		double edge_dist = sqrt((coordInfoVecs.nodeLocX[node1] - coordInfoVecs.nodeLocX[node2])*(coordInfoVecs.nodeLocX[node1] - coordInfoVecs.nodeLocX[node2]) + 
					(coordInfoVecs.nodeLocY[node1] - coordInfoVecs.nodeLocY[node2])*(coordInfoVecs.nodeLocY[node1] - coordInfoVecs.nodeLocY[node2]) +
				(coordInfoVecs.nodeLocZ[node1] - coordInfoVecs.nodeLocZ[node2])*(coordInfoVecs.nodeLocZ[node1] - coordInfoVecs.nodeLocZ[node2]));
		if (RULES_OF_EXPAN == 1){
			if ((This_area_v/ areaTriangleInfoVecs.initial_area >= EXPAN_THRESHOLD) && (This_area_w/ areaTriangleInfoVecs.initial_area >= EXPAN_THRESHOLD)){
				triggered = true;
			}
			else{continue;}
		}
		else if (RULES_OF_EXPAN == 2){
			if ((This_area_v/ areaTriangleInfoVecs.initial_area >= EXPAN_THRESHOLD) && (This_area_w/ areaTriangleInfoVecs.initial_area >= EXPAN_THRESHOLD_weak)){
				triggered = true;
			}
			else if ((This_area_v/ areaTriangleInfoVecs.initial_area >= EXPAN_THRESHOLD_weak) && (This_area_w/ areaTriangleInfoVecs.initial_area >= EXPAN_THRESHOLD)){
				triggered = true;
			}
			//else if (edge_dist >= 2.0*generalParams.Rmin){
			//	triggered = true;
			//}
			else{continue;}
		}
		else if (RULES_OF_EXPAN == 3){
			//std::cout<<"This_area_v + This_area_w = "<<This_area_v + This_area_w<<std::endl;
			if ((This_area_v + This_area_w)/(2.0*areaTriangleInfoVecs.initial_area) >= EXPAN_THRESHOLD_weak){
				//std::cout<<4.0*areaTriangleInfoVecs.initial_area<<std::endl;
				triggered = true;
			}
			else{continue;}
		}
		
		////std::cout<<"GROWTH ERROR 2"<<std::endl;	
		int t1e1, t1e2, t1e3, t2e1, t2e2, t2e3;

		if (coordInfoVecs.triangles2Edges_1[elem1] == iedge){
			t1e1 = coordInfoVecs.triangles2Edges_2[elem1];
			t1e2 = coordInfoVecs.triangles2Edges_3[elem1];
			//t1e3 = coordInfoVecs.triangles2Edges_1[elem1];
		}
		else if (coordInfoVecs.triangles2Edges_2[elem1] == iedge){
			t1e1 = coordInfoVecs.triangles2Edges_3[elem1];
			t1e2 = coordInfoVecs.triangles2Edges_1[elem1];
			//t1e3 = coordInfoVecs.triangles2Edges_2[elem1];
		} 
		else if (coordInfoVecs.triangles2Edges_3[elem1] == iedge){
			t1e1 = coordInfoVecs.triangles2Edges_1[elem1];
			t1e2 = coordInfoVecs.triangles2Edges_2[elem1];
			//t1e3 = coordInfoVecs.triangles2Edges_3[elem1];
		}
		////std::cout<<"GROWTH ERROR 3"<<std::endl;	

		if (coordInfoVecs.triangles2Edges_1[elem2] == iedge){
			t2e1 = coordInfoVecs.triangles2Edges_2[elem2];
			t2e2 = coordInfoVecs.triangles2Edges_3[elem2];
			//t2e3 = coordInfoVecs.triangles2Edges_1[elem2];
		}
		else if (coordInfoVecs.triangles2Edges_2[elem2] == iedge){
			t2e1 = coordInfoVecs.triangles2Edges_3[elem2];
			t2e2 = coordInfoVecs.triangles2Edges_1[elem2];
			//t2e3 = coordInfoVecs.triangles2Edges_2[elem2];
		} 
		else if (coordInfoVecs.triangles2Edges_3[elem2] == iedge){
			t2e1 = coordInfoVecs.triangles2Edges_1[elem2];
			t2e2 = coordInfoVecs.triangles2Edges_2[elem2];
			//t2e3 = coordInfoVecs.triangles2Edges_3[elem2];
		}
		//Note that in the above assignment, t1e3 and t2e3 are the edges shared by the neighboring triangles T1 and T2.
		////std::cout<<"GROWTH ERROR 4"<<std::endl;	

		
		int n1, n2, n3, n4;
		
		if ((coordInfoVecs.edges2Nodes_1[t1e1] == coordInfoVecs. edges2Nodes_1[iedge]) || (coordInfoVecs.edges2Nodes_1[t1e1] == coordInfoVecs. edges2Nodes_2[iedge]) ){
			n1 = coordInfoVecs.edges2Nodes_1[t1e1];
			n2 = coordInfoVecs.edges2Nodes_2[t1e1];
			if (coordInfoVecs.edges2Nodes_1[t1e1] == coordInfoVecs. edges2Nodes_1[iedge]){
				n3 = coordInfoVecs.edges2Nodes_2[iedge];
			}
			else if (coordInfoVecs.edges2Nodes_1[t1e1] == coordInfoVecs. edges2Nodes_2[iedge]){
				n3 = coordInfoVecs.edges2Nodes_1[iedge];
			}
		}
		else if ((coordInfoVecs.edges2Nodes_2[t1e1] == coordInfoVecs. edges2Nodes_1[iedge]) || (coordInfoVecs.edges2Nodes_2[t1e1] == coordInfoVecs. edges2Nodes_2[iedge]) ){
			n1 = coordInfoVecs.edges2Nodes_2[t1e1];
			n2 = coordInfoVecs.edges2Nodes_1[t1e1];
			if (coordInfoVecs.edges2Nodes_2[t1e1] == coordInfoVecs. edges2Nodes_1[iedge]){
				n3 = coordInfoVecs.edges2Nodes_2[iedge];
			}
			else if (coordInfoVecs.edges2Nodes_2[t1e1] == coordInfoVecs. edges2Nodes_2[iedge]){
				n3 = coordInfoVecs.edges2Nodes_1[iedge];
			}
		}
		////std::cout<<"GROWTH ERROR 5"<<std::endl;	

		if (coordInfoVecs.edges2Nodes_1[t2e1] == coordInfoVecs.edges2Nodes_1[iedge] || coordInfoVecs.edges2Nodes_1[t2e1] == coordInfoVecs.edges2Nodes_2[iedge]){
			n4 = coordInfoVecs.edges2Nodes_2[t2e1];
		}
		else if (coordInfoVecs.edges2Nodes_2[t2e1] == coordInfoVecs.edges2Nodes_1[iedge] || coordInfoVecs.edges2Nodes_2[t2e1] == coordInfoVecs.edges2Nodes_2[iedge]){
			n4 = coordInfoVecs.edges2Nodes_1[t2e1];
		}
		int safe_growth1 = 0;
		int safe_growth2 = 0;
		if (coordInfoVecs. nndata1[n2] >= 0){safe_growth1 += 1;        }
        if (coordInfoVecs. nndata2[n2] >= 0){safe_growth1 += 1;        }
        if (coordInfoVecs. nndata3[n2] >= 0){safe_growth1 += 1;        }
        if (coordInfoVecs. nndata4[n2] >= 0){safe_growth1 += 1;        }
        if (coordInfoVecs. nndata5[n2] >= 0){safe_growth1 += 1;        }
        if (coordInfoVecs. nndata6[n2] >= 0){safe_growth1 += 1;        }
        if (coordInfoVecs. nndata7[n2] >= 0){safe_growth1 += 1;        }
        if (coordInfoVecs. nndata8[n2] >= 0){safe_growth1 += 1;        }
        if (coordInfoVecs. nndata9[n2] >= 0){safe_growth1 += 1;        }
        if (coordInfoVecs. nndata10[n2] >= 0){safe_growth1 += 1;        }
        if (coordInfoVecs. nndata11[n2] >= 0){safe_growth1 += 1;        }
		if (coordInfoVecs. nndata12[n2] >= 0){safe_growth1 += 1;        }
		if (coordInfoVecs. nndata1[n4] >= 0){safe_growth2 += 1;        }
        if (coordInfoVecs. nndata2[n4] >= 0){safe_growth2 += 1;        }
        if (coordInfoVecs. nndata3[n4] >= 0){safe_growth2 += 1;        }
        if (coordInfoVecs. nndata4[n4] >= 0){safe_growth2 += 1;        }
        if (coordInfoVecs. nndata5[n4] >= 0){safe_growth2 += 1;        }
        if (coordInfoVecs. nndata6[n4] >= 0){safe_growth2 += 1;        }
        if (coordInfoVecs. nndata7[n4] >= 0){safe_growth2 += 1;        }
        if (coordInfoVecs. nndata8[n4] >= 0){safe_growth2 += 1;        }
        if (coordInfoVecs. nndata9[n4] >= 0){safe_growth2 += 1;        }
        if (coordInfoVecs. nndata10[n4] >= 0){safe_growth2 += 1;        }
        if (coordInfoVecs. nndata11[n4] >= 0){safe_growth2 += 1;        }
		if (coordInfoVecs. nndata12[n4] >= 0){safe_growth2 += 1;        }
		
		if (safe_growth1 >= generalParams.safeguardthreshold || safe_growth2 >= generalParams.safeguardthreshold){
			continue;
		}

		//std::cout<<"n1 = "<<n1<<std::endl;
		//std::cout<<"n2 = "<<n2<<std::endl;
		//std::cout<<"n3 = "<<n3<<std::endl;
		//std::cout<<"n4 = "<<n4<<std::endl;
		//These extract the indices of vertices of the selected triangles "elem1" and "elem2". Now we have n1, n2, n3, n4 in the correct orientation (supposedly).

		////std::cout<<"GROWTH ERROR 6"<<std::endl;	
		int edgeindex, a, a1, a2, a3, temp1, temp2;
		//std::cout<<"maxNodeCount = "<< generalParams.maxNodeCount<<std::endl;
		double newx = (coordInfoVecs.nodeLocX[coordInfoVecs.edges2Nodes_1[iedge]] + coordInfoVecs.nodeLocX[coordInfoVecs.edges2Nodes_2[iedge]])/2.0;
		////std::cout<<"GROWTH ERROR 6.1"<<std::endl;	
		coordInfoVecs.nodeLocX[generalParams. maxNodeCount] = newx;
		////std::cout<<"GROWTH ERROR 6.2"<<std::endl;	
		double newy = (coordInfoVecs.nodeLocY[coordInfoVecs.edges2Nodes_1[iedge]] + coordInfoVecs.nodeLocY[coordInfoVecs.edges2Nodes_2[iedge]])/2.0;
		////std::cout<<"GROWTH ERROR 6.3"<<std::endl;	
		coordInfoVecs.nodeLocY[generalParams. maxNodeCount] = newy;
		////std::cout<<"GROWTH ERROR 6.4"<<std::endl;	
		double newz = (coordInfoVecs.nodeLocZ[coordInfoVecs.edges2Nodes_1[iedge]] + coordInfoVecs.nodeLocZ[coordInfoVecs.edges2Nodes_2[iedge]])/2.0;
		////std::cout<<"GROWTH ERROR 6.5"<<std::endl;	
		coordInfoVecs.nodeLocZ[generalParams. maxNodeCount] = newz;
		//These are the coordinate of the new vertex. Its index is "coordInfoVecs.nodeLocX.size()-1"

		//Before editing major data structures, we will update the nndata here since it is only affected by the addition of new nodes.

		//int NODESIZE= generalParams.maxNodeCount;//coordInfoVecs.nodeLocX.size();
		////std::cout<<"GROWTH ERROR 7"<<std::endl;			
		coordInfoVecs.triangles2Nodes_1[coordInfoVecs.num_triangles] = n1;
		coordInfoVecs.triangles2Nodes_2[coordInfoVecs.num_triangles] = n2;
		coordInfoVecs.triangles2Nodes_3[coordInfoVecs.num_triangles] = generalParams.maxNodeCount;
		coordInfoVecs.num_triangles += 1;
		//NOTE: What this +1 actually does is that it specifies the location to write
		//any new data. Here it points to the location to write new triangles information.
		//This is a new triangle associated with (tn1, tn2, newnode). Its index is "coordInfoVecs.triangles2Nodes_1.size()-4".
		////std::cout<<"GROWTH ERROR 8"<<std::endl;	
		coordInfoVecs.triangles2Nodes_1[coordInfoVecs.num_triangles] =(n2);
		coordInfoVecs.triangles2Nodes_2[coordInfoVecs.num_triangles] =(n3);
		coordInfoVecs.triangles2Nodes_3[coordInfoVecs.num_triangles] = generalParams.maxNodeCount;
		coordInfoVecs.num_triangles += 1;
		//This is a new triangle associated with (tn2, tn3, newnode). Its index is "coordInfoVecs.triangles2Nodes_1.size()-3".
		////std::cout<<"GROWTH ERROR 9"<<std::endl;	
		coordInfoVecs.triangles2Nodes_1[coordInfoVecs.num_triangles] =(n3);
		coordInfoVecs.triangles2Nodes_2[coordInfoVecs.num_triangles] =(n4);
		coordInfoVecs.triangles2Nodes_3[coordInfoVecs.num_triangles] = generalParams.maxNodeCount;
		coordInfoVecs.num_triangles += 1;
		//This is a new triangle associated with (tn3, tn1, newnode). Its index is "coordInfoVecs.triangles2Nodes_1.size()-2".
		////std::cout<<"GROWTH ERROR 10"<<std::endl;	
		coordInfoVecs.triangles2Nodes_1[coordInfoVecs.num_triangles] =(n4);
		coordInfoVecs.triangles2Nodes_2[coordInfoVecs.num_triangles] =(n1);
		coordInfoVecs.triangles2Nodes_3[coordInfoVecs.num_triangles] = generalParams.maxNodeCount;
		coordInfoVecs.num_triangles += 1;
		//This is a new triangle associated with (tn3, tn1, newnode). Its index is "coordInfoVecs.triangles2Nodes_1.size()-1".
		////std::cout<<"GROWTH ERROR 11"<<std::endl;	
		//Now we add new edges formed by the addition of the new node.
		coordInfoVecs.edges2Nodes_1[coordInfoVecs.num_edges] = (generalParams.maxNodeCount);
		coordInfoVecs.edges2Nodes_2[coordInfoVecs.num_edges] = (n1);
		coordInfoVecs.edges2Triangles_1[coordInfoVecs.num_edges] = coordInfoVecs.num_triangles - 4;
		coordInfoVecs.edges2Triangles_2[coordInfoVecs.num_edges] = coordInfoVecs.num_triangles - 1;
		coordInfoVecs.num_edges += 1;
		////std::cout<<"GROWTH ERROR 12"<<std::endl;	
		//This is a new edge associated with (newnode, tn1). Its index is "edges2Nodes_1.size()-4".
		coordInfoVecs.edges2Nodes_1[coordInfoVecs.num_edges] = (generalParams.maxNodeCount);
		coordInfoVecs.edges2Nodes_2[coordInfoVecs.num_edges] = (n2);
		coordInfoVecs.edges2Triangles_1[coordInfoVecs.num_edges] = coordInfoVecs.num_triangles - 3;
		coordInfoVecs.edges2Triangles_2[coordInfoVecs.num_edges] = coordInfoVecs.num_triangles - 4;
		coordInfoVecs.num_edges += 1;
		////std::cout<<"GROWTH ERROR 13"<<std::endl;	
		//This is a new edge associated with (newnode, tn2). Its index is "edges2Nodes_1.size()-3".
		coordInfoVecs.edges2Nodes_1[coordInfoVecs.num_edges] = (generalParams.maxNodeCount);
		coordInfoVecs.edges2Nodes_2[coordInfoVecs.num_edges] = (n3);
		coordInfoVecs.edges2Triangles_1[coordInfoVecs.num_edges] = coordInfoVecs.num_triangles - 2;
		coordInfoVecs.edges2Triangles_2[coordInfoVecs.num_edges] = coordInfoVecs.num_triangles - 3;
		coordInfoVecs.num_edges += 1;
		////std::cout<<"GROWTH ERROR 14"<<std::endl;	
		//This is a new edge associated with (newnode, tn3). Its index is "edges2Nodes_1.size()-2".
		coordInfoVecs.edges2Nodes_1[coordInfoVecs.num_edges] = (generalParams.maxNodeCount);
		coordInfoVecs.edges2Nodes_2[coordInfoVecs.num_edges] = (n4);
		coordInfoVecs.edges2Triangles_1[coordInfoVecs.num_edges] = coordInfoVecs.num_triangles - 1;
		coordInfoVecs.edges2Triangles_2[coordInfoVecs.num_edges] = coordInfoVecs.num_triangles - 2;
		coordInfoVecs.num_edges += 1;
		////std::cout<<"GROWTH ERROR 15"<<std::endl;	
		for (int j = 0; j < 4; j++){
		//	//std::cout<<"GROWTH ERROR 16"<<std::endl;				
			//Now we check to see if the order of update is correct, i.e. are edges2Triangles data in correct orientation.
			//This is crucial in the bendingspring computation.
			edgeindex = (coordInfoVecs.num_edges - (4-j));
			a = coordInfoVecs.edges2Triangles_1[edgeindex];
			if ((coordInfoVecs.triangles2Nodes_1[a] == coordInfoVecs.edges2Nodes_1[edgeindex]) && (coordInfoVecs.triangles2Nodes_2[a] == coordInfoVecs.edges2Nodes_2[edgeindex])){
				a1 = 1;
			}
			else{
				a1 = 0;
			}
			if ((coordInfoVecs.triangles2Nodes_2[a] == coordInfoVecs.edges2Nodes_1[edgeindex]) && (coordInfoVecs.triangles2Nodes_3[a] == coordInfoVecs.edges2Nodes_2[edgeindex])){
				a2 = 1;
			}
			else{
				a2 = 0;
			}
			if ((coordInfoVecs.triangles2Nodes_3[a] == coordInfoVecs.edges2Nodes_1[edgeindex]) && (coordInfoVecs.triangles2Nodes_1[a] == coordInfoVecs.edges2Nodes_2[edgeindex])){
				a3 = 1;
			}
			else{
				a3 = 0;
			}

			if ((a1+a2+a3) == 0){
				temp1 = coordInfoVecs.edges2Triangles_1[edgeindex];
				temp2 = coordInfoVecs.edges2Triangles_2[edgeindex];
				coordInfoVecs.edges2Triangles_1[edgeindex] = temp2;
				coordInfoVecs.edges2Triangles_2[edgeindex] = temp1;
			}
			else{}
			//This checks if the orientation is correct or not, if not, flip the ordering.
		}
		//This is a new edge associated with (newnode, tn3). Its index is "edges2Nodes_1.size()-1".
		generalParams.maxNodeCount += 1;

		coordInfoVecs.nndata1[generalParams.maxNodeCount-1] =  (n1);
		coordInfoVecs.nndata2[generalParams.maxNodeCount-1] =  (n2);
		coordInfoVecs.nndata3[generalParams.maxNodeCount-1] =  (n3);
		coordInfoVecs.nndata4[generalParams.maxNodeCount-1] =  (n4);
		coordInfoVecs.nndata5[generalParams.maxNodeCount-1] =  (-2);
		coordInfoVecs.nndata6[generalParams.maxNodeCount-1] =  (-2);
		coordInfoVecs.nndata7[generalParams.maxNodeCount-1] =  (-2);
		coordInfoVecs.nndata8[generalParams.maxNodeCount-1] =  (-2);
		coordInfoVecs.nndata9[generalParams.maxNodeCount-1] =  (-2);
		coordInfoVecs.nndata10[generalParams.maxNodeCount-1] = (-2);
		coordInfoVecs.nndata11[generalParams.maxNodeCount-1] = (-2);
		coordInfoVecs.nndata12[generalParams.maxNodeCount-1] = (-2);
		for (int j = 0; j < 2; j++){
			int nn, nnn, nnnn;
			if (j == 0){
				nn = n1;
				nnn = n3;
				nnnn = generalParams.maxNodeCount-1;
			}
			else if (j == 1){
				nn = n3;
				nnn = n1;
				nnnn = generalParams.maxNodeCount-1;
			}
			if (coordInfoVecs.nndata1[nn] == nnn){
				coordInfoVecs.nndata1[nn] = nnnn;
			}
			else if (coordInfoVecs.nndata2[nn] == nnn){
				coordInfoVecs.nndata2[nn] = nnnn;
			}
			else if (coordInfoVecs.nndata3[nn] == nnn){
				coordInfoVecs.nndata3[nn] = nnnn;
			}
			else if (coordInfoVecs.nndata4[nn] == nnn){
				coordInfoVecs.nndata4[nn] = nnnn;
			}
			else if (coordInfoVecs.nndata5[nn] == nnn){
				coordInfoVecs.nndata5[nn] = nnnn;
			}
			else if (coordInfoVecs.nndata6[nn] == nnn){
				coordInfoVecs.nndata6[nn] = nnnn;
			}
			else if (coordInfoVecs.nndata7[nn] == nnn){
				coordInfoVecs.nndata7[nn] = nnnn;
			}
			else if (coordInfoVecs.nndata8[nn] == nnn){
				coordInfoVecs.nndata8[nn] = nnnn;
			}
			else if (coordInfoVecs.nndata9[nn] == nnn){
				coordInfoVecs.nndata9[nn] = nnnn;
			}
			else if (coordInfoVecs.nndata10[nn] == nnn){
				coordInfoVecs.nndata10[nn] = nnnn;
			}
			else if (coordInfoVecs.nndata11[nn] == nnn){
				coordInfoVecs.nndata11[nn] = nnnn;
			}
			else if (coordInfoVecs.nndata12[nn] == nnn){
				coordInfoVecs.nndata12[nn] = nnnn;
			}
		}

		for (int j = 0; j < 2; j++){
			int nn, nnn;
			if (j == 0){
				nn = n2;
				nnn = generalParams.maxNodeCount-1;
			}
			else if (j == 1){
				nn = n4;
				nnn = generalParams.maxNodeCount-1;
			}
			if (coordInfoVecs.nndata1[nn] < 0){
				coordInfoVecs.nndata1[nn] = nnn;
			}
			else if (coordInfoVecs.nndata2[nn] < 0){
				coordInfoVecs.nndata2[nn] = nnn;
			}
			else if (coordInfoVecs.nndata3[nn] < 0){
				coordInfoVecs.nndata3[nn] = nnn;
			}
			else if (coordInfoVecs.nndata4[nn] < 0){
				coordInfoVecs.nndata4[nn] = nnn;
			}
			else if (coordInfoVecs.nndata5[nn] < 0){
				coordInfoVecs.nndata5[nn] = nnn;
			}
			else if (coordInfoVecs.nndata6[nn] < 0){
				coordInfoVecs.nndata6[nn] = nnn;
			}
			else if (coordInfoVecs.nndata7[nn] < 0){
				coordInfoVecs.nndata7[nn] = nnn;
			}
			else if (coordInfoVecs.nndata8[nn] < 0){
				coordInfoVecs.nndata8[nn] = nnn;
			}
			else if (coordInfoVecs.nndata9[nn] < 0){
				coordInfoVecs.nndata9[nn] = nnn;
			}
			else if (coordInfoVecs.nndata10[nn] < 0){
				coordInfoVecs.nndata10[nn] = nnn;
			}
			else if (coordInfoVecs.nndata11[nn] < 0){
				coordInfoVecs.nndata11[nn] = nnn;
			}
			else if (coordInfoVecs.nndata12[nn] < 0){
				coordInfoVecs.nndata12[nn] = nnn;
			}
		}
		//generalParams.num_of_nodes += 1;

		


		////std::cout<<"GROWTH ERROR 17"<<std::endl;	
		//Now we update the edges2Triangles data structure with new edges.
		//std::cout<<"elem 1 = "<<elem1<<std::endl;
		//std::cout<<"elem 2 = "<<elem2<<std::endl;
		for (int i = 0; i < coordInfoVecs.num_edges; i++){
		//	std::cout<<"edges2triangles"<<" "<< i <<" : "<<coordInfoVecs. edges2Triangles_1[i]<<" "<<coordInfoVecs. edges2Triangles_2[i]<<std::endl;
		}
		int TRIANGLESIZE = coordInfoVecs.num_triangles;//coordInfoVecs.triangles2Nodes_1.size();
		if (coordInfoVecs.edges2Triangles_1[t1e1] == elem1){
			coordInfoVecs.edges2Triangles_1[t1e1] = TRIANGLESIZE-4;
		}
		else if (coordInfoVecs.edges2Triangles_2[t1e1] == elem1){
			coordInfoVecs.edges2Triangles_2[t1e1] = TRIANGLESIZE-4;
		}
		else{}
		////std::cout<<"GROWTH ERROR 18"<<std::endl;	
		if (coordInfoVecs.edges2Triangles_1[t1e2] == elem1){
			coordInfoVecs.edges2Triangles_1[t1e2] = TRIANGLESIZE-3;
		}
		else if (coordInfoVecs.edges2Triangles_2[t1e2] == elem1){
			coordInfoVecs.edges2Triangles_2[t1e2] = TRIANGLESIZE-3;
		}
		else{}
		////std::cout<<"GROWTH ERROR 19"<<std::endl;	
		if (coordInfoVecs.edges2Triangles_1[t2e1] == elem2){
			coordInfoVecs.edges2Triangles_1[t2e1] = TRIANGLESIZE-2;
		}
		else if (coordInfoVecs.edges2Triangles_2[t2e1] == elem2){
			coordInfoVecs.edges2Triangles_2[t2e1] = TRIANGLESIZE-2;
		}
		else{}
		////std::cout<<"GROWTH ERROR 20"<<std::endl;	
		if (coordInfoVecs.edges2Triangles_1[t2e2] == elem2){
			coordInfoVecs.edges2Triangles_1[t2e2] = TRIANGLESIZE-1;
		}
		else if (coordInfoVecs.edges2Triangles_2[t2e2] == elem2){
			coordInfoVecs.edges2Triangles_2[t2e2] = TRIANGLESIZE-1;
		}
		else{}
		//std::cout<<"t1e1 "<<t1e1<<std::endl;
		//std::cout<<"t1e2 "<<t1e2<<std::endl;
		//std::cout<<"t1e3 "<<t1e3<<std::endl;
		//std::cout<<"t2e1 "<<t2e1<<std::endl;
		//std::cout<<"t2e2 "<<t2e2<<std::endl;
		//std::cout<<"t2e3 "<<t2e3<<std::endl;

		//for (int i = 0; i < coordInfoVecs.num_edges; i++){
		//	std::cout<<"edges2triangles"<<" "<< i <<" : "<<coordInfoVecs. edges2Triangles_1[i]<<" "<<coordInfoVecs. edges2Triangles_2[i]<<std::endl;
		//}
		//The above change the existing edges2Triangles data structure to accomodate new triangles added.

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		
		//Now we will take care of the last unedited data structure "triangles2Edges".
		//int aa, bb;
		int EDGESIZE = coordInfoVecs.num_edges;//coordInfoVecs.edges2Nodes_1.size();
		for (int j = 0; j < 4; j++){
		//	//std::cout<<"GROWTH ERROR 21"<<std::endl;	
			if (j == 0){
				coordInfoVecs.triangles2Edges_1[coordInfoVecs.num_triangles - 4] = (EDGESIZE-4);
				coordInfoVecs.triangles2Edges_2[coordInfoVecs.num_triangles - 4] = (t1e1);
				coordInfoVecs.triangles2Edges_3[coordInfoVecs.num_triangles - 4] = (EDGESIZE-3);   
			}
			else if (j == 1){
				coordInfoVecs.triangles2Edges_1[coordInfoVecs.num_triangles - 3] = (EDGESIZE-3);
				coordInfoVecs.triangles2Edges_2[coordInfoVecs.num_triangles - 3] = (t1e2);
				coordInfoVecs.triangles2Edges_3[coordInfoVecs.num_triangles - 3] = (EDGESIZE-2);   
			}
			else if (j ==2){
				coordInfoVecs.triangles2Edges_1[coordInfoVecs.num_triangles - 2] = (EDGESIZE-2);
				coordInfoVecs.triangles2Edges_2[coordInfoVecs.num_triangles - 2] = (t2e1);
				coordInfoVecs.triangles2Edges_3[coordInfoVecs.num_triangles - 2] = (EDGESIZE-1);   
			}
			else if (j ==3){
				coordInfoVecs.triangles2Edges_1[coordInfoVecs.num_triangles - 1] = (EDGESIZE-1);
				coordInfoVecs.triangles2Edges_2[coordInfoVecs.num_triangles - 1] = (t2e2);
				coordInfoVecs.triangles2Edges_3[coordInfoVecs.num_triangles - 1] = (EDGESIZE-4);   
			}
			
		}
	
		
		if (generalParams.nodes_in_upperhem[coordInfoVecs.edges2Nodes_1[iedge]] == 1 && generalParams.nodes_in_upperhem[coordInfoVecs.edges2Nodes_2[iedge]] == 1){
			generalParams.nodes_in_upperhem.push_back(1);
		}
		else{
			generalParams.nodes_in_upperhem.push_back(-1);
		}
		//Finally, we will fill the edge data chosen for growth (expansion) with INT_MAX so its data is no longer relevant to the computation
		////std::cout<<"GROWTH ERROR 22"<<std::endl;	
		coordInfoVecs.edges2Nodes_1[iedge] = INT_MAX;
		coordInfoVecs.edges2Nodes_2[iedge] = INT_MAX;
		for (int i = 0; i < coordInfoVecs.num_triangles; i++){
		//	//std::cout<<"GROWTH ERROR 23"<<std::endl;	
			if (coordInfoVecs.triangles2Edges_1[i] == iedge){
				coordInfoVecs.triangles2Edges_1[i] = INT_MAX;
			}
			if (coordInfoVecs.triangles2Edges_2[i] == iedge){
				coordInfoVecs.triangles2Edges_2[i] = INT_MAX;
			}
			if (coordInfoVecs.triangles2Edges_3[i] == iedge){
				coordInfoVecs.triangles2Edges_3[i] = INT_MAX;
			}
		}
		coordInfoVecs.edges2Triangles_1[iedge] = INT_MAX;
		coordInfoVecs.edges2Triangles_2[iedge] = INT_MAX;
		
		////std::cout<<"GROWTH ERROR 24"<<std::endl;	
		
			coordInfoVecs.triangles2Nodes_1[elem1] = INT_MAX;
			coordInfoVecs.triangles2Nodes_2[elem1] = INT_MAX;
			coordInfoVecs.triangles2Nodes_3[elem1] = INT_MAX;
			coordInfoVecs.triangles2Nodes_1[elem2] = INT_MAX;
			coordInfoVecs.triangles2Nodes_2[elem2] = INT_MAX;
			coordInfoVecs.triangles2Nodes_3[elem2] = INT_MAX;
			
			//Delete the associated vertices information of the selected triangle.
			//Since we delete the chosen triangles, any triangle indexed lower than the deleted one will have its index reduced (or moved up) by 1.
			//Hence, we need to sweep through all data structures using the triangle index to change the index accordingly.
			for (int i = 0; i < coordInfoVecs.num_edges; i++){
		//		//std::cout<<"GROWTH ERROR 25"<<std::endl;	
				if (coordInfoVecs.edges2Triangles_1[i] == elem1 || coordInfoVecs.edges2Triangles_1[i] == elem2){
					coordInfoVecs.edges2Triangles_1[i] = INT_MAX;
				}
				if (coordInfoVecs.edges2Triangles_2[i] == elem1 || coordInfoVecs.edges2Triangles_2[i] == elem2){
					coordInfoVecs.edges2Triangles_2[i] = INT_MAX;
				}
			if (coordInfoVecs.edges2Triangles_1[i] != INT_MAX && coordInfoVecs.edges2Triangles_2[i] == INT_MAX){
				std::cout<<"modified edges2Triangles "<<coordInfoVecs.edges2Triangles_1[i]<<" "<<coordInfoVecs.edges2Triangles_2[i]<<std::endl;
				}
				else if (coordInfoVecs.edges2Triangles_1[i] == INT_MAX && coordInfoVecs.edges2Triangles_2[i] != INT_MAX){
					std::cout<<"modified edges2Triangles "<<coordInfoVecs.edges2Triangles_1[i]<<" "<<coordInfoVecs.edges2Triangles_2[i]<<std::endl;
					}
			}
			//This completes the sweep. After this, the indices of triangle used in edges2Triangles data structure should be the correct one.
		//	//std::cout<<"GROWTH ERROR 26"<<std::endl;	
			coordInfoVecs.triangles2Edges_1[elem1] = INT_MAX;
			coordInfoVecs.triangles2Edges_2[elem1] = INT_MAX;
			coordInfoVecs.triangles2Edges_3[elem1] = INT_MAX;
			coordInfoVecs.triangles2Edges_1[elem2] = INT_MAX;
			coordInfoVecs.triangles2Edges_2[elem2] = INT_MAX;
			coordInfoVecs.triangles2Edges_3[elem2] = INT_MAX;
			for (int i = 0; i < coordInfoVecs.num_triangles; i++){
		//		//std::cout<<"GROWTH ERROR 27"<<std::endl;	
				if (coordInfoVecs.triangles2Edges_1[i] == iedge){
					coordInfoVecs.triangles2Edges_1[i] = INT_MAX;
				}
				if (coordInfoVecs.triangles2Edges_2[i] == iedge ){
					coordInfoVecs.triangles2Edges_2[i] = INT_MAX;
				}
				if (coordInfoVecs.triangles2Edges_3[i] == iedge ){
					coordInfoVecs.triangles2Edges_3[i] = INT_MAX;
				}
			}
		
		//Erase the edge infomation related to the deleted triangle. Note the deletion should always start with the largest index.

		//Before we delete the edge, determine whether the newly added node is part of nodes_in_upperhem or not.
		

		
						//Erase the edge infomation related to the deleted triangle.

						//Now we update the nodes_in_upperhem and edges_in_upperhem data structures.
						//This ensures that the newly created edges will have the correct associated spring constant.
//std::cout<<"ERROR HERE?"<<std::endl;
		generalParams.edges_in_upperhem[iedge] = INT_MAX;
		for (int i = 0; i < edges_in_upperhem.size(); i++){
			if (edges_in_upperhem[i] == iedge){
				edges_in_upperhem[i] == INT_MAX;
				//break;
			}
		}

		for (int q = 0; q < 4; q++){
		//	//std::cout<<"GROWTH ERROR 30"<<std::endl;	
			int nodeP = coordInfoVecs.triangles2Nodes_1[coordInfoVecs.num_triangles - (4-q)]; 
			int nodeQ = coordInfoVecs.triangles2Nodes_2[coordInfoVecs.num_triangles - (4-q)];
			int nodeR = coordInfoVecs.triangles2Nodes_3[coordInfoVecs.num_triangles - (4-q)];

			if (generalParams.nodes_in_upperhem[nodeP]==1 && generalParams.nodes_in_upperhem[nodeQ] ==1 && generalParams.nodes_in_upperhem[nodeR] ==1){
				generalParams.triangles_in_upperhem.push_back(1);
			}
			else if (generalParams.nodes_in_upperhem[nodeP]==1 && generalParams.nodes_in_upperhem[nodeQ] ==1){
				generalParams.triangles_in_upperhem.push_back(0);
			}
			else if (generalParams.nodes_in_upperhem[nodeP]==1 && generalParams.nodes_in_upperhem[nodeR] ==1){
				generalParams.triangles_in_upperhem.push_back(0);
			}
			else if (generalParams.nodes_in_upperhem[nodeQ] ==1 && generalParams.nodes_in_upperhem[nodeR] ==1){
				generalParams.triangles_in_upperhem.push_back(0);
			}
			else{
				generalParams.triangles_in_upperhem.push_back(INT_MAX);
			}
		}
		//std::cout<<"edges2Triangles size"<<""<<coordInfoVecs.edges2Triangles_1.size()<<" "<<coordInfoVecs.edges2Triangles_2.size()<<std::endl;
		//std::cout<<"triangles_in_upperhem size "<<generalParams.triangles_in_upperhem.size()<<std::endl;	
		//std::cout<<"GROWTH ERROR 29"<<std::endl;	
		for (int q = 0; q < 4; q++){
			//std::cout<<"GROWTH ERROR 31"<<std::endl;	
			int elem_1 = coordInfoVecs.edges2Triangles_1[coordInfoVecs.num_edges-(4 - q)];
			//std::cout<<coordInfoVecs.num_edges-(4 - q)<<std::endl;
			//std::cout<<"elem_1 "<<elem_1<<std::endl;
			//std::cout<<generalParams.nodes_in_upperhem[nodeP]<<std::endl;
			int elem_2 = coordInfoVecs.edges2Triangles_2[coordInfoVecs.num_edges-(4 - q)];
			//std::cout<<"elem_2"<<elem_2<<std::endl;
			//std::cout<<generalParams.nodes_in_upperhem[nodeQ]<<std::endl;
			//std::cout<<"GROWTH ERROR 31.5"<<std::endl;
			if (generalParams.triangles_in_upperhem[elem_1] == 1 && generalParams.triangles_in_upperhem[elem_2] == 1){
				
				generalParams.edges_in_upperhem.push_back(1);
				//generalParams.edges_in_upperhem_index.push_back(generalParams.num_of_edges - (4 - q));
				edges_in_upperhem.push_back(coordInfoVecs.num_edges - (4 - q));
			}
			
			else if (generalParams.triangles_in_upperhem[elem_1] == 1 || generalParams.triangles_in_upperhem[elem_2] == 1){
				
				generalParams.edges_in_upperhem.push_back(0);
				//generalParams.edges_in_upperhem_index.push_back(generalParams.num_of_edges - (4 - q));
			}
			else{
				
				generalParams.edges_in_upperhem.push_back(-1);
			}

			generalParams.boundaries_in_upperhem.push_back(-1);
			
			
		}
		generalParams.triangles_in_upperhem[elem1] = INT_MAX;
		generalParams.triangles_in_upperhem[elem2] = INT_MAX;

		//nodes_in_tip.push_back(-1);
		
						
						

						//Here we regenerate the edges that we will loop through in both edgeswap and growth, if any edge deletion actually happened.
						
						//if (triggered_counter > 0){
						//	edges_in_upperhem_for_loop.clear();
						//	for (int i = 0; i < generalParams.edges_in_upperhem_index.size(); i++){
						//		if (generalParams.edges_in_upperhem_index[i] != INT_MAX){
						//			edges_in_upperhem_for_loop.push_back(generalParams.edges_in_upperhem_index[i]);
						//		}
						//	}
						//}
						
						//This should completes the dreadful data structure update associated with cell (membrane) growth.
						//Have fun modifying it if you need more function!
						//if (triggered == true){
						//	break;
						//}
					}

					//Recalculate the nodes_in_tip data.
				if (triggered == true){	
					/*if (edgeswap_iteration >= 0){
					//generalParams.nodes_in_upperhem.clear();
					//generalParams.nodes_in_upperhem.resize(generalParams.maxNodeCount,-1);
					double max_z = -2000.0;
					for (int k = 0; k < generalParams.maxNodeCount; k++){
						if (coordInfoVecs. nodeLocZ[k] >= max_z){
							max_z = coordInfoVecs. nodeLocZ[k];
						}
					}

					for (int i = 0; i < generalParams.maxNodeCount; i++){
						if (coordInfoVecs.nodeLocZ[i] > (max_z - tip_depth)){
							generalParams.nodes_in_upperhem[i] = 1;
						}
						else{
							generalParams.nodes_in_upperhem[i] = -1;
						}
					//	std::cout<<"nodes "<<i<<" "<<generalParams.nodes_in_upperhem[i]<<std::endl;		
					}*/
				
					//std::vector<int> nodes_to_center;
				
					////////////////////////////////////////////////////////////////////////////////////////////
				
				
				
					////////////////////////////////////////////////////////////////////////////////////////////
					//generalParams.triangles_in_upperhem.resize(coordInfoVecs.num_triangles);
					/*for (int i = 0; i < coordInfoVecs.num_triangles; i++){
						if (coordInfoVecs.triangles2Nodes_1[i] == INT_MAX || coordInfoVecs.triangles2Nodes_2[i] == INT_MAX || coordInfoVecs.triangles2Nodes_3[i] == INT_MAX){
							generalParams.triangles_in_upperhem[i] = INT_MAX;
							continue;
						}
						int aaa = generalParams.nodes_in_upperhem[coordInfoVecs.triangles2Nodes_1[i]];
						//std::cout<<aaa<<std::endl;
						int bbb = generalParams.nodes_in_upperhem[coordInfoVecs.triangles2Nodes_2[i]];
						//std::cout<<bbb<<std::endl;
						int ccc = generalParams.nodes_in_upperhem[coordInfoVecs.triangles2Nodes_3[i]];
						//std::cout<<ccc<<std::endl;
						if ((aaa+bbb+ccc)==3){
							generalParams.triangles_in_upperhem[i] = 1;
							//triangles_in_upperhem.push_back(i);
						}
						else if ((aaa+bbb+ccc)==1){
							generalParams.triangles_in_upperhem[i] = 0;
							//triangles_in_upperhem.push_back(i);
						}
						else{
							generalParams.triangles_in_upperhem[i] = -1;
						}
					//	std::cout<<"triangle "<<i<<" "<<generalParams.triangles_in_upperhem[i]<<std::endl;		
					}
					//std::cout<<"WHERE iS THE PROBLEM 1"<<std::endl;
				
					//std::vector<int> edges_in_upperhem;
					edges_in_upperhem.clear();
					//generalParams.edges_in_upperhem.resize(coordInfoVecs.num_edges);
					for (int i = 0; i < coordInfoVecs.num_edges; i++){
						
						if (coordInfoVecs.edges2Nodes_1[i] == INT_MAX || coordInfoVecs.edges2Nodes_2[i] == INT_MAX){
							generalParams.edges_in_upperhem[i] = INT_MAX;
							continue;
						}
						
						//std::cout<<"edges2Triangles_1 = "<<coordInfoVecs.edges2Triangles_1[i]<<std::endl;
						//std::cout<<"edges2Triangles_2 = "<<coordInfoVecs.edges2Triangles_2[i]<<std::endl;
						int aaa = generalParams.triangles_in_upperhem[coordInfoVecs.edges2Triangles_1[i]];//generalParams.nodes_in_upperhem[coordInfoVecs.edges2Nodes_1[i]];
						int bbb = generalParams.triangles_in_upperhem[coordInfoVecs.edges2Triangles_2[i]];//generalParams.nodes_in_upperhem[coordInfoVecs.edges2Nodes_2[i]];
						
						if (aaa == 1 && bbb == 1){
							generalParams.edges_in_upperhem[i] = 1;
							edges_in_upperhem.push_back(i);
						}
						else if (aaa == 1 || bbb == 1){
							generalParams.edges_in_upperhem[i] = 1;//0;
						}
						else{
							generalParams.edges_in_upperhem[i] = -1;
						}
						
					}
				}*/
					//std::cout<<"WHERE iS THE PROBLEM 2"<<std::endl;
					
					
					/*for (int k = 0; k < generalParams.maxNodeCount; k++){
						if (coordInfoVecs. nodeLocZ[k] >= (max_z - tip_depth)){
							nodes_in_tip[k] == 1;
						}
						else{
							nodes_in_tip[k] == -1;
						}
					}*/
					true_num_edges_in_upperhem = 0;
					for (int i = 0; i < edges_in_upperhem.size(); i++){
						if (edges_in_upperhem[i] != INT_MAX && edges_in_upperhem[i] >= 0){
							true_num_edges_in_upperhem += 1;
							//break;
						}
					}
					//std::cout<<"WHERE iS THE PROBLEM 3"<<std::endl;
				}
			
			
			
			
//std::cout<<"GROWTH DONE!"<<std::endl;
 ////storage->print_VTK_File();
////storage->storeVariables();

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////// END OF GROWTH SECTION //////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
					//max_height = -10000.0;
                      //                  min_height = 10000.0;
                        //                for (int k = 0; k < generalParams.maxNodeCount; k++){
                          //                      if (coordInfoVecs. nodeLocZ[k] >= max_height){
                            //                            max_height = coordInfoVecs. nodeLocZ[k];
                              //                  }
                                //                if (coordInfoVecs.nodeLocZ[k] <= min_height){
                                  //                      min_height = coordInfoVecs.nodeLocZ[k];
                                    //            }
                                      //  }
					/*VectorShuffleForFilamentLoop.clear();
					for (int i = 0; i < generalParams.maxNodeCount; i++){
						VectorShuffleForFilamentLoop.push_back(i);
					}
					//std::random_device rand_dev;
					std::mt19937 generator3(rand_dev());
					std::shuffle(std::begin(VectorShuffleForFilamentLoop), std::end(VectorShuffleForFilamentLoop), generator3);
					int maxElementIndex = std::max_element(coordInfoVecs.nodeLocZ.begin(),coordInfoVecs.nodeLocZ.end()) - coordInfoVecs.nodeLocZ.begin();
					for (int i = 0; i < generalParams.maxNodeCountLJ; i++){
						if (i == 0){
							filament_base[i] = maxElementIndex;
							continue;
						}
						for (int j = 0; j < VectorShuffleForFilamentLoop.size(); j++){
							double xsquared = (ljInfoVecs.LJ_PosX_all[i] - coordInfoVecs.nodeLocX[j])*
												(ljInfoVecs.LJ_PosX_all[i] - coordInfoVecs.nodeLocX[j]);
							double ysquared = (ljInfoVecs.LJ_PosY_all[i] - coordInfoVecs.nodeLocY[j])*
												(ljInfoVecs.LJ_PosY_all[i] - coordInfoVecs.nodeLocY[j]);
							double zsquared = (ljInfoVecs.LJ_PosZ_all[i] - coordInfoVecs.nodeLocZ[j])*
												(ljInfoVecs.LJ_PosZ_all[i] - coordInfoVecs.nodeLocZ[j]);
							double R = sqrt(xsquared + ysquared + zsquared);
							if (R < filament_Rmin*1.1  && j != maxElementIndex){
								filament_base[i] = j;
								break;
							}
							else{filament_base[i] = -1;}
						}
					}*/ 

ComputeVolume(
	generalParams,
	coordInfoVecs,
	linearSpringInfoVecs,
	ljInfoVecs);
//std::cout<<"ERROR 1"<<std::endl;

/*edgeswap_ptr->transferDtoH(coordInfoVecs, build_ptr->hostSetInfoVecs);
//std::cout<<"ERROR 1.5"<<std::endl;
for (int edge_loop = 0; edge_loop < num_edge_loop; edge_loop++) {
	//std::cout<<"edge_loop = "<<edge_loop<<std::endl;
	
	std::random_device rand_dev;
	std::mt19937 generator(rand_dev());
   
   std::uniform_int_distribution<int> distribution(1,edges_in_upperhem.size());
   
   int dice_roll = distribution(generator);
   
   int edge = edges_in_upperhem[dice_roll - 1];
   
   while (generalParams.boundaries_in_upperhem[edge] == 1 || edge == INT_MAX){
		dice_roll = distribution(generator);
		
		edge =  edges_in_upperhem[dice_roll - 1];
	 }
   //std::cout<<"edge = "<<edge<<std::endl;
	int ALPHA = edgeswap_ptr->edge_swap_host_vecs(
		edge,
		generalParams,
		build_ptr->hostSetInfoVecs,
		linearSpringInfoVecs,
		bendingTriangleInfoVecs,
		areaTriangleInfoVecs);
	
}
//NOTE: EDGESWAP ALGORITHM CURRENTLY IS WRITTEN TO ALLOW AT MOST 8 NEIGHBORING NODES PER NODE.
//std::cout<<"ERROR 2"<<std::endl;
edgeswap_ptr->transferHtoD(coordInfoVecs, build_ptr->hostSetInfoVecs);//Currently this is treated as a backup of coordInfoVecs
*/


 					
					
 			}
		
		}
		

	};
	
	





void System::assignStorage(std::shared_ptr<Storage> _storage) {
	storage = _storage;
};
void System::set_weak_builder(std::weak_ptr<SystemBuilder> _weak_bld_ptr) {
	weak_bld_ptr = _weak_bld_ptr;
};



//initialize memory for thrust vectors and set coordInfoVecs vals from input. 
void System::initializeSystem(HostSetInfoVecs& hostSetInfoVecs) {
	std::cout<<"Initializing"<<std::endl;

	generalParams.maxNodeCount = hostSetInfoVecs.nodeLocX.size();
	coordInfoVecs.num_edges = hostSetInfoVecs.edges2Nodes_1.size();
	coordInfoVecs.num_triangles = hostSetInfoVecs.triangles2Nodes_1.size();

	std::cout<<"num nodes: "<< generalParams.maxNodeCount << std::endl;
	std::cout<<"num edges: "<< coordInfoVecs.num_edges << std::endl;
	std::cout<<"num elems: "<< coordInfoVecs.num_triangles << std::endl;
	//allocate memory
	int mem_prealloc = 4;
	coordInfoVecs.isNodeFixed.resize(mem_prealloc*hostSetInfoVecs.nodeLocX.size(),false);
	coordInfoVecs.prevNodeLocX.resize(mem_prealloc*hostSetInfoVecs.nodeLocX.size());
	coordInfoVecs.prevNodeLocY.resize(mem_prealloc*hostSetInfoVecs.nodeLocX.size());
	coordInfoVecs.prevNodeLocZ.resize(mem_prealloc*hostSetInfoVecs.nodeLocX.size());

	coordInfoVecs.prevNodeForceX.resize(mem_prealloc*hostSetInfoVecs.nodeLocX.size());
	coordInfoVecs.prevNodeForceY.resize(mem_prealloc*hostSetInfoVecs.nodeLocX.size());
	coordInfoVecs.prevNodeForceZ.resize(mem_prealloc*hostSetInfoVecs.nodeLocX.size());
	
	coordInfoVecs.nodeLocX.resize(mem_prealloc*hostSetInfoVecs.nodeLocX.size());
	coordInfoVecs.nodeLocY.resize(mem_prealloc*hostSetInfoVecs.nodeLocX.size());
	coordInfoVecs.nodeLocZ.resize(mem_prealloc*hostSetInfoVecs.nodeLocX.size());

	coordInfoVecs.nodeForceX.resize(mem_prealloc*hostSetInfoVecs.nodeLocX.size(), 0.0);
	coordInfoVecs.nodeForceY.resize(mem_prealloc*hostSetInfoVecs.nodeLocX.size(), 0.0);
	coordInfoVecs.nodeForceZ.resize(mem_prealloc*hostSetInfoVecs.nodeLocX.size(), 0.0);

	coordInfoVecs.triangles2Nodes_1.resize( mem_prealloc*coordInfoVecs.num_triangles );
	coordInfoVecs.triangles2Nodes_2.resize( mem_prealloc*coordInfoVecs.num_triangles );
	coordInfoVecs.triangles2Nodes_3.resize( mem_prealloc*coordInfoVecs.num_triangles );
	
	coordInfoVecs.triangles2Edges_1.resize( mem_prealloc*coordInfoVecs.num_triangles );
	coordInfoVecs.triangles2Edges_2.resize( mem_prealloc*coordInfoVecs.num_triangles );
	coordInfoVecs.triangles2Edges_3.resize( mem_prealloc*coordInfoVecs.num_triangles );

	coordInfoVecs.edges2Nodes_1.resize( mem_prealloc*coordInfoVecs.num_edges );
	coordInfoVecs.edges2Nodes_2.resize( mem_prealloc*coordInfoVecs.num_edges );
	
	coordInfoVecs.edges2Triangles_1.resize( mem_prealloc*coordInfoVecs.num_edges );
	coordInfoVecs.edges2Triangles_2.resize( mem_prealloc*coordInfoVecs.num_edges );

	coordInfoVecs.nndata1.resize( mem_prealloc*generalParams.maxNodeCount);
	coordInfoVecs.nndata2.resize( mem_prealloc*generalParams.maxNodeCount);
	coordInfoVecs.nndata3.resize( mem_prealloc*generalParams.maxNodeCount);
	coordInfoVecs.nndata4.resize( mem_prealloc*generalParams.maxNodeCount);
	coordInfoVecs.nndata5.resize( mem_prealloc*generalParams.maxNodeCount);
	coordInfoVecs.nndata6.resize( mem_prealloc*generalParams.maxNodeCount);
	coordInfoVecs.nndata7.resize( mem_prealloc*generalParams.maxNodeCount);
	coordInfoVecs.nndata8.resize( mem_prealloc*generalParams.maxNodeCount);
	coordInfoVecs.nndata9.resize( mem_prealloc*generalParams.maxNodeCount);
	coordInfoVecs.nndata10.resize( mem_prealloc*generalParams.maxNodeCount);
	coordInfoVecs.nndata11.resize( mem_prealloc*generalParams.maxNodeCount);
	coordInfoVecs.nndata12.resize( mem_prealloc*generalParams.maxNodeCount);



	//copy info to GPU
	std::cout<<"Copying"<<std::endl;
	thrust::copy(hostSetInfoVecs.isNodeFixed.begin(),hostSetInfoVecs.isNodeFixed.end(), coordInfoVecs.isNodeFixed.begin());
	
	std::cout<<"fixed_node_in_host: "<<std::endl;
	for (int k = 0; k < hostSetInfoVecs.isNodeFixed.size(); k++){
		//std::cout<<hostSetInfoVecs.isNodeFixed[k]<<std::endl;
	}
	std::cout<<"end_of_fixed_node_host_printout"<<std::endl;
	std::cout<<"fixed_node_in_device: "<<std::endl;
	for (int k = 0; k < coordInfoVecs.isNodeFixed.size(); k++){
		//std::cout<<coordInfoVecs.isNodeFixed[k]<<std::endl;
	}
	std::cout<<"end_of_fixed_node_device_printout"<<std::endl;
std::cout<<"size of host fixed "<< hostSetInfoVecs.isNodeFixed.size()<<std::endl;
std::cout<<"size of device fixed "<< coordInfoVecs.isNodeFixed.size()<<std::endl;

	/*for (int k = 0; k < coordInfoVecs.isNodeFixed.size(); k++){
		bool isFixedHost = hostSetInfoVecs.isNodeFixed[k];
		bool isFixedDevice = coordInfoVecs.isNodeFixed[k];
		if (isFixedDevice != isFixedHost){

			std::cout<<"pos "<< k << " dev val = " << coordInfoVecs.isNodeFixed[k]
				<< " host val = " <<  hostSetInfoVecs.isNodeFixed[k] <<std::endl;
		}
	}*/
	thrust::fill(coordInfoVecs.nodeForceX.begin(), coordInfoVecs.nodeForceX.end(), 0.0);
	thrust::fill(coordInfoVecs.nodeForceY.begin(), coordInfoVecs.nodeForceY.end(), 0.0);
	thrust::fill(coordInfoVecs.nodeForceZ.begin(), coordInfoVecs.nodeForceZ.end(), 0.0);

	thrust::fill(coordInfoVecs.prevNodeForceX.begin(), coordInfoVecs.prevNodeForceX.end(), 0.0);
	thrust::fill(coordInfoVecs.prevNodeForceY.begin(), coordInfoVecs.prevNodeForceY.end(), 0.0);
	thrust::fill(coordInfoVecs.prevNodeForceZ.begin(), coordInfoVecs.prevNodeForceZ.end(), 0.0);
	
	thrust::copy(hostSetInfoVecs.nodeLocX.begin(), hostSetInfoVecs.nodeLocX.end(), coordInfoVecs.prevNodeLocX.begin() );
	thrust::copy(hostSetInfoVecs.nodeLocY.begin(), hostSetInfoVecs.nodeLocY.end(), coordInfoVecs.prevNodeLocY.begin() );
	thrust::copy(hostSetInfoVecs.nodeLocZ.begin(), hostSetInfoVecs.nodeLocZ.end(), coordInfoVecs.prevNodeLocZ.begin() );
	
	thrust::copy(hostSetInfoVecs.nodeLocX.begin(), hostSetInfoVecs.nodeLocX.end(), coordInfoVecs.nodeLocX.begin() );
	thrust::copy(hostSetInfoVecs.nodeLocY.begin(), hostSetInfoVecs.nodeLocY.end(), coordInfoVecs.nodeLocY.begin() );
	thrust::copy(hostSetInfoVecs.nodeLocZ.begin(), hostSetInfoVecs.nodeLocZ.end(), coordInfoVecs.nodeLocZ.begin() );
	
	thrust::copy(hostSetInfoVecs.triangles2Nodes_1.begin(), hostSetInfoVecs.triangles2Nodes_1.end(), coordInfoVecs.triangles2Nodes_1.begin() );
	thrust::copy(hostSetInfoVecs.triangles2Nodes_2.begin(), hostSetInfoVecs.triangles2Nodes_2.end(), coordInfoVecs.triangles2Nodes_2.begin() );
	thrust::copy(hostSetInfoVecs.triangles2Nodes_3.begin(), hostSetInfoVecs.triangles2Nodes_3.end(), coordInfoVecs.triangles2Nodes_3.begin() );
	
	thrust::copy(hostSetInfoVecs.triangles2Edges_1.begin(), hostSetInfoVecs.triangles2Edges_1.end(), coordInfoVecs.triangles2Edges_1.begin() );
	thrust::copy(hostSetInfoVecs.triangles2Edges_2.begin(), hostSetInfoVecs.triangles2Edges_2.end(), coordInfoVecs.triangles2Edges_2.begin() );
	thrust::copy(hostSetInfoVecs.triangles2Edges_3.begin(), hostSetInfoVecs.triangles2Edges_3.end(), coordInfoVecs.triangles2Edges_3.begin() );

	thrust::copy(hostSetInfoVecs.edges2Nodes_1.begin(), hostSetInfoVecs.edges2Nodes_1.end(), coordInfoVecs.edges2Nodes_1.begin() );
	thrust::copy(hostSetInfoVecs.edges2Nodes_2.begin(), hostSetInfoVecs.edges2Nodes_2.end(), coordInfoVecs.edges2Nodes_2.begin() );
	
	thrust::copy(hostSetInfoVecs.edges2Triangles_1.begin(), hostSetInfoVecs.edges2Triangles_1.end(), coordInfoVecs.edges2Triangles_1.begin() );
	thrust::copy(hostSetInfoVecs.edges2Triangles_2.begin(), hostSetInfoVecs.edges2Triangles_2.end(), coordInfoVecs.edges2Triangles_2.begin() );

	thrust::copy(hostSetInfoVecs.nndata1.begin(), hostSetInfoVecs.nndata1.end(), coordInfoVecs.nndata1.begin() );
	thrust::copy(hostSetInfoVecs.nndata2.begin(), hostSetInfoVecs.nndata2.end(), coordInfoVecs.nndata2.begin() );
	thrust::copy(hostSetInfoVecs.nndata3.begin(), hostSetInfoVecs.nndata3.end(), coordInfoVecs.nndata3.begin() );
	thrust::copy(hostSetInfoVecs.nndata4.begin(), hostSetInfoVecs.nndata4.end(), coordInfoVecs.nndata4.begin() );
	thrust::copy(hostSetInfoVecs.nndata5.begin(), hostSetInfoVecs.nndata5.end(), coordInfoVecs.nndata5.begin() );
	thrust::copy(hostSetInfoVecs.nndata6.begin(), hostSetInfoVecs.nndata6.end(), coordInfoVecs.nndata6.begin() );
	thrust::copy(hostSetInfoVecs.nndata7.begin(), hostSetInfoVecs.nndata7.end(), coordInfoVecs.nndata7.begin() );
	thrust::copy(hostSetInfoVecs.nndata8.begin(), hostSetInfoVecs.nndata8.end(), coordInfoVecs.nndata8.begin() );
	thrust::copy(hostSetInfoVecs.nndata9.begin(), hostSetInfoVecs.nndata9.end(), coordInfoVecs.nndata9.begin() );
	thrust::copy(hostSetInfoVecs.nndata10.begin(), hostSetInfoVecs.nndata10.end(), coordInfoVecs.nndata10.begin() );
	thrust::copy(hostSetInfoVecs.nndata11.begin(), hostSetInfoVecs.nndata11.end(), coordInfoVecs.nndata11.begin() );
	thrust::copy(hostSetInfoVecs.nndata12.begin(), hostSetInfoVecs.nndata12.end(), coordInfoVecs.nndata12.begin() );


 
	//allocate memory for other data structures.   

	//area triangle info vec
	//number of area springs is the number of triangles
	std::cout<<"Mem"<<std::endl;
	areaTriangleInfoVecs.tempNodeIdUnreduced.resize(mem_prealloc*areaTriangleInfoVecs.factor * coordInfoVecs.num_triangles);
	areaTriangleInfoVecs.tempNodeForceXUnreduced.resize(mem_prealloc*areaTriangleInfoVecs.factor * coordInfoVecs.num_triangles);
	areaTriangleInfoVecs.tempNodeForceYUnreduced.resize(mem_prealloc*areaTriangleInfoVecs.factor * coordInfoVecs.num_triangles);
	areaTriangleInfoVecs.tempNodeForceZUnreduced.resize(mem_prealloc*areaTriangleInfoVecs.factor * coordInfoVecs.num_triangles);
	
	areaTriangleInfoVecs.tempNodeIdReduced.resize(mem_prealloc*areaTriangleInfoVecs.factor * coordInfoVecs.num_triangles);
	areaTriangleInfoVecs.tempNodeForceXReduced.resize(mem_prealloc*areaTriangleInfoVecs.factor * coordInfoVecs.num_triangles);
	areaTriangleInfoVecs.tempNodeForceYReduced.resize(mem_prealloc*areaTriangleInfoVecs.factor * coordInfoVecs.num_triangles);
	areaTriangleInfoVecs.tempNodeForceZReduced.resize(mem_prealloc*areaTriangleInfoVecs.factor * coordInfoVecs.num_triangles);

	//beinding triangle info vec
	//num bending springs is the number of times each edge is between two triangles. 
	bendingTriangleInfoVecs.numBendingSprings = coordInfoVecs.num_edges;//coordInfoVecs.edges2Triangles_1.size();

	bendingTriangleInfoVecs.tempNodeIdUnreduced.resize(mem_prealloc*bendingTriangleInfoVecs.factor * bendingTriangleInfoVecs.numBendingSprings);
	bendingTriangleInfoVecs.tempNodeForceXUnreduced.resize(mem_prealloc*bendingTriangleInfoVecs.factor * bendingTriangleInfoVecs.numBendingSprings);
	bendingTriangleInfoVecs.tempNodeForceYUnreduced.resize(mem_prealloc*bendingTriangleInfoVecs.factor * bendingTriangleInfoVecs.numBendingSprings);
	bendingTriangleInfoVecs.tempNodeForceZUnreduced.resize(mem_prealloc*bendingTriangleInfoVecs.factor * bendingTriangleInfoVecs.numBendingSprings);
	
	bendingTriangleInfoVecs.tempNodeIdReduced.resize(mem_prealloc*bendingTriangleInfoVecs.factor * bendingTriangleInfoVecs.numBendingSprings);
	bendingTriangleInfoVecs.tempNodeForceXReduced.resize(mem_prealloc*bendingTriangleInfoVecs.factor * bendingTriangleInfoVecs.numBendingSprings);
	bendingTriangleInfoVecs.tempNodeForceYReduced.resize(mem_prealloc*bendingTriangleInfoVecs.factor * bendingTriangleInfoVecs.numBendingSprings);
	bendingTriangleInfoVecs.tempNodeForceZReduced.resize(mem_prealloc*bendingTriangleInfoVecs.factor * bendingTriangleInfoVecs.numBendingSprings);

	//linear springs
	
	linearSpringInfoVecs.tempNodeIdUnreduced.resize(mem_prealloc*linearSpringInfoVecs.factor*coordInfoVecs.num_edges);
	linearSpringInfoVecs.tempNodeForceXUnreduced.resize(mem_prealloc*linearSpringInfoVecs.factor*coordInfoVecs.num_edges);
	linearSpringInfoVecs.tempNodeForceYUnreduced.resize(mem_prealloc*linearSpringInfoVecs.factor*coordInfoVecs.num_edges);
	linearSpringInfoVecs.tempNodeForceZUnreduced.resize(mem_prealloc*linearSpringInfoVecs.factor*coordInfoVecs.num_edges);
	
	linearSpringInfoVecs.tempNodeIdReduced.resize(mem_prealloc*linearSpringInfoVecs.factor * coordInfoVecs.num_edges);
	linearSpringInfoVecs.tempNodeForceXReduced.resize(mem_prealloc*linearSpringInfoVecs.factor * coordInfoVecs.num_edges);
	linearSpringInfoVecs.tempNodeForceYReduced.resize(mem_prealloc*linearSpringInfoVecs.factor * coordInfoVecs.num_edges);
	linearSpringInfoVecs.tempNodeForceZReduced.resize(mem_prealloc*linearSpringInfoVecs.factor * coordInfoVecs.num_edges);
	
	linearSpringInfoVecs.edge_initial_length.clear();
	linearSpringInfoVecs.edge_initial_length.resize(mem_prealloc*coordInfoVecs.num_edges,1.0);
	
	thrust::copy(hostSetInfoVecs.edge_initial_length.begin(), hostSetInfoVecs.edge_initial_length.end(), linearSpringInfoVecs.edge_initial_length.begin() );

	//Resize the hostSetInfoVecs so that we can copy data back and forth between hostSetinfoVecs and coordInfoVecs without problem.
	hostSetInfoVecs.isNodeFixed.resize(mem_prealloc*hostSetInfoVecs.nodeLocX.size());
	//hostSetInfoVecs.prevNodeLocX.resize(mem_prealloc*hostSetInfoVecs.nodeLocX.size());
	//hostSetInfoVecs.prevNodeLocY.resize(mem_prealloc*hostSetInfoVecs.nodeLocX.size());
	//hostSetInfoVecs.prevNodeLocZ.resize(mem_prealloc*hostSetInfoVecs.nodeLocX.size());

	//hostSetInfoVecs.prevNodeForceX.resize(mem_prealloc*hostSetInfoVecs.nodeLocX.size());
	//hostSetInfoVecs.prevNodeForceY.resize(mem_prealloc*hostSetInfoVecs.nodeLocX.size());
	//hostSetInfoVecs.prevNodeForceZ.resize(mem_prealloc*hostSetInfoVecs.nodeLocX.size());
	
	hostSetInfoVecs.nodeLocX.resize(coordInfoVecs.nodeLocX.size());
	hostSetInfoVecs.nodeLocY.resize(coordInfoVecs.nodeLocX.size());
	hostSetInfoVecs.nodeLocZ.resize(coordInfoVecs.nodeLocX.size());
	std::cout<<"Host_nodeLocX size = "<<hostSetInfoVecs.nodeLocX.size()<<std::endl;

	hostSetInfoVecs.nodeForceX.resize(coordInfoVecs.nodeLocX.size());
	hostSetInfoVecs.nodeForceY.resize(coordInfoVecs.nodeLocX.size());
	hostSetInfoVecs.nodeForceZ.resize(coordInfoVecs.nodeLocX.size());
	std::cout<<"Host_nodeForceX size = "<<hostSetInfoVecs.nodeLocX.size()<<std::endl;

	hostSetInfoVecs.triangles2Nodes_1.resize( coordInfoVecs.triangles2Nodes_1.size() );
	hostSetInfoVecs.triangles2Nodes_2.resize( coordInfoVecs.triangles2Nodes_2.size() );
	hostSetInfoVecs.triangles2Nodes_3.resize( coordInfoVecs.triangles2Nodes_3.size() );
	std::cout<<"Host_triangles2Nodes size = "<<hostSetInfoVecs.triangles2Nodes_1.size()<<std::endl;
	
	hostSetInfoVecs.triangles2Edges_1.resize( coordInfoVecs.triangles2Edges_1.size() );
	hostSetInfoVecs.triangles2Edges_2.resize( coordInfoVecs.triangles2Edges_2.size() );
	hostSetInfoVecs.triangles2Edges_3.resize( coordInfoVecs.triangles2Edges_3.size() );
	std::cout<<"Host_triangles2Edges size = "<<hostSetInfoVecs.triangles2Edges_1.size()<<std::endl;

	hostSetInfoVecs.edges2Nodes_1.resize( coordInfoVecs.edges2Nodes_1.size() );
	hostSetInfoVecs.edges2Nodes_2.resize( coordInfoVecs.edges2Nodes_2.size() );
	std::cout<<"Host_edges2Nodes size = "<<hostSetInfoVecs.edges2Nodes_1.size()<<std::endl;
	
	hostSetInfoVecs.edges2Triangles_1.resize( coordInfoVecs.edges2Triangles_1.size() );
	hostSetInfoVecs.edges2Triangles_2.resize( coordInfoVecs.edges2Triangles_2.size() );
	std::cout<<"Host_edges2Triangles size = "<<hostSetInfoVecs.edges2Triangles_1.size()<<std::endl;

	hostSetInfoVecs.nndata1.resize( mem_prealloc*generalParams.maxNodeCount);
	hostSetInfoVecs.nndata2.resize( mem_prealloc*generalParams.maxNodeCount);
	hostSetInfoVecs.nndata3.resize( mem_prealloc*generalParams.maxNodeCount);
	hostSetInfoVecs.nndata4.resize( mem_prealloc*generalParams.maxNodeCount);
	hostSetInfoVecs.nndata5.resize( mem_prealloc*generalParams.maxNodeCount);
	hostSetInfoVecs.nndata6.resize( mem_prealloc*generalParams.maxNodeCount);
	hostSetInfoVecs.nndata7.resize( mem_prealloc*generalParams.maxNodeCount);
	hostSetInfoVecs.nndata8.resize( mem_prealloc*generalParams.maxNodeCount);
	hostSetInfoVecs.nndata9.resize( mem_prealloc*generalParams.maxNodeCount);
	hostSetInfoVecs.nndata10.resize( mem_prealloc*generalParams.maxNodeCount);
	hostSetInfoVecs.nndata11.resize( mem_prealloc*generalParams.maxNodeCount);
	hostSetInfoVecs.nndata12.resize( mem_prealloc*generalParams.maxNodeCount);

	std::cout<<"initial lengths: "<< linearSpringInfoVecs.edge_initial_length.size()<<std::endl;

	std::cout<<"System Ready"<<std::endl;

	//Generate LJ particle list. and set LJ particle midpoint.
	//double maxX_lj = *(thrust::max_element(coordInfoVecs.nodeLocX.begin(),coordInfoVecs.nodeLocX.end()));
	//double minX_lj = *(thrust::min_element(coordInfoVecs.nodeLocX.begin(),coordInfoVecs.nodeLocX.end()));
	//double maxY_lj = *(thrust::max_element(coordInfoVecs.nodeLocY.begin(),coordInfoVecs.nodeLocY.end()));
	//double minY_lj = *(thrust::min_element(coordInfoVecs.nodeLocY.begin(),coordInfoVecs.nodeLocY.end()));
	
	//ljInfoVecs.LJ_PosX = (maxX_lj + minX_lj)/2.0;
	//ljInfoVecs.LJ_PosY = (maxY_lj + minY_lj)/2.0;


	//currently unused
	/*thrust::host_vector<int> tempIds;
	for (int i = 0; i < hostSetInfoVecs.nodeLocX.size(); i++ ) {
		double xLoc = hostSetInfoVecs.nodeLocX[i];
		double yLoc = hostSetInfoVecs.nodeLocY[i];
		double zLoc = hostSetInfoVecs.nodeLocZ[i];
		
		double xDist = ljInfoVecs.LJ_PosX - xLoc;
		double yDist = ljInfoVecs.LJ_PosY - yLoc;
		double zDist = ljInfoVecs.LJ_PosZ - zLoc;

		double dist = std::sqrt(xDist*xDist + yDist*yDist + zDist*zDist);
		//just test all poitns for now. Optimize later.
		if (dist < ljInfoVecs.Rcutoff) {
			tempIds.push_back(i);
		}
	}
	ljInfoVecs.node_id_close.resize( tempIds.size() );
	thrust::copy(tempIds.begin(), tempIds.end(), ljInfoVecs.node_id_close.begin());
	std::cout<<"lj nodes: "<< ljInfoVecs.node_id_close.size() << std::endl;*/






	//last, set memory foor buckets.
	auxVecs.id_bucket.resize(generalParams.maxNodeCount);
	auxVecs.id_value.resize(generalParams.maxNodeCount);
	auxVecs.id_bucket_expanded.resize(27 * (generalParams.maxNodeCount));
	auxVecs.id_value_expanded.resize(27 *( generalParams.maxNodeCount ));
 


};


