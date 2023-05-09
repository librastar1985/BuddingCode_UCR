#include <stdio.h>
#include "System.h"
#include "SystemStructures.h" 
#include "AreaTriangles.h"
//#include "AreaTrianglesEnergy.h"
#include "BendingTriangles.h"
// #include "BendingTrianglesEnergy.h"
#include "MemRepulsionSprings_universal.h"
#include "MemRepulsionSprings_local.h"
#include "MemRepulsionEnergy.h"
#include "LinearSprings.h"
//#include "LinearSpringsEnergy.h"
//#include "LJSprings.h"
//#include "LJSprings_LJ.h"
#include "NodeAdvance.h"
//#include "BucketScheme.h"
#include "Storage.h" 
#include "Utilities.h"
#include "SystemBuilder.h"
#include <vector>
#include "VolumeComp.h"
#include "VolumeSprings.h"
#include <bits/stdc++.h>
#include "LineTensionSprings.h"
//#include "Growth.h"
#include <math.h>
#include "vmalloc.h"
#include <list>
#include "TurgorForce.h"

//#include "SurfaceNormal.h"
//#include "Nodes2Triangles.h"



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
	// std::cout<<"ERROR 1"<<std::endl;
	ComputeAreaTriangleSprings(
		
		generalParams,
		coordInfoVecs,
		areaTriangleInfoVecs);
		// std::cout<<"ERROR 2"<<std::endl;
	ComputeTurgorSprings(
		generalParams,
		coordInfoVecs,
		areaTriangleInfoVecs
	);
	// std::cout<<"ERROR 3"<<std::endl;
	ComputeCosTriangleSprings(
		
		generalParams,
		coordInfoVecs,  
		bendingTriangleInfoVecs); 
		// std::cout<<"ERROR 4"<<std::endl;
	// ComputeMemRepulsionSprings_universal(
	// 	coordInfoVecs,
	// 	linearSpringInfoVecs, 
	// 	capsidInfoVecs,
	// 	generalParams,
	// 	auxVecs);

	ComputeMemRepulsionSprings_local(
		coordInfoVecs,
		linearSpringInfoVecs, 
		capsidInfoVecs,
		generalParams,
		auxVecs);
		// std::cout<<"ERROR 5"<<std::endl;

	ComputeVolume(
		generalParams,
		coordInfoVecs,
		linearSpringInfoVecs,
		ljInfoVecs);

		// std::cout<<"ERROR 6"<<std::endl;
	/*ComputeVolumeSprings(
		coordInfoVecs,
		linearSpringInfoVecs, 
		capsidInfoVecs,
		generalParams,
		auxVecs);*/

	/* if (generalParams.true_current_total_volume/initial_volume >= 1.25){
	ComputeLineTensionSprings(
		generalParams,
		coordInfoVecs,
		linearSpringInfoVecs);
	} */
		
};

void System::solveSystem(){

	coordInfoVecs.k_0 = 20.0;
	coordInfoVecs.k_1 = 25.0;
	coordInfoVecs.k_2 = 5.0;
	coordInfoVecs.k_3 = 5.0;
	coordInfoVecs.k_4 = 1.0;
	coordInfoVecs.k_ss = 12;//10.75;
	coordInfoVecs.beta = 1.0/1.0;///1.45;
	coordInfoVecs.gamma = 1.0;
	coordInfoVecs.q1 = 10.0;
	coordInfoVecs.h = 10.0;

	uint mem_prealloc = 4; //Make sure that this number is the same as set in System::initializeSystem found near the bottom of this script.

	double max_conc_scaler_for_material_insert = 0.8;
	std::cout<<"multiplier applied to max conc to determine location of material insertion = "<<max_conc_scaler_for_material_insert<<std::endl;
	double current_edge_to_tip_height_scale = 2.0;
	std::cout<<"current_edge_to_tip_height_scale = "<<current_edge_to_tip_height_scale<<std::endl;
	//Determines how far away from the tip can new material be inserted.
	// double current_edge_to_tip_dist_scale = 4.0;
	// std::cout<<"current_edge_to_tip_dist_scale = "<<current_edge_to_tip_dist_scale<<std::endl;
	double bdry_to_tip_height_scale = 0.0;
	std::cout<<"bdry_to_tip_height_scale = "<<bdry_to_tip_height_scale<<std::endl;


	// double u_scalingPower = 8;//6;//8;//4;//6.0;
	double u_scalingPower = 7; // For test-compilation purpose, please delete this line and reactive the line above.
	std::cout<<"u_scalingPower (power of 'u' in the chem diff) = "<<u_scalingPower<<std::endl;
	double max_u_scalingPower = 8;
	std::cout<<"max u_scalingPower = "<<max_u_scalingPower<<std::endl;
	double max_u_scalingPower_postPeak = 8;
	std::cout<<"max_u_scalingPower_postPeak = "<<max_u_scalingPower_postPeak<<std::endl;
	double timeFrame = 5e3;
	std::cout<<"timeFrame needed to reach max u_scalingPower = "<<timeFrame<<std::endl;
	double timeFrame_postPeak = 5e3;
	std::cout<<"timeFrame needed to reach max u_scalingPower_postPeak = "<<timeFrame_postPeak<<std::endl;
	double u_scalingPower_Progress = sqrt((max_u_scalingPower - u_scalingPower)*(max_u_scalingPower - u_scalingPower))/timeFrame;
	double u_scalingPower_postPeak_Progress = sqrt((max_u_scalingPower - max_u_scalingPower_postPeak)*(max_u_scalingPower - max_u_scalingPower_postPeak))/timeFrame_postPeak;
	bool polarization_postPeak = false;
	bool triggered = false;
	generalParams.current_total_sim_step = 0;
	int relax_max_steps_before_growth_and_edgeswap = 3e3;
	Front *front = new Front();
	RBC *rbc = new RBC();
	RBC *n_rbc = new RBC();
	auto utilities_ptr = std::make_shared<Utilities>(coordInfoVecs, generalParams);
    auto build_ptr = weak_bld_ptr.lock();
	std::cout<<"Declaration of rbc and n_rbc complete."<<std::endl;
	std::cout<<"Utilities_ptr declaration complete."<<std::endl;
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////Build the "u" vector representing the external or internal influencer for polarization /////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	double max_cell_triangle_height, min_cell_triangle_height, v1,v2,v3,v4, cell_height;
	max_cell_triangle_height = -10000.0;
	min_cell_triangle_height = 10000.0;

	for (int i = 0; i < coordInfoVecs.num_triangles; i++){
		if (coordInfoVecs.triangles2Nodes_1[i] >= (INT_MAX-100) || coordInfoVecs.triangles2Nodes_2[i] >= (INT_MAX-100) || coordInfoVecs.triangles2Nodes_3[i] >= (INT_MAX-100)){
			continue;
		}
		if (coordInfoVecs.triangles2Nodes_1[i] <= (-INT_MAX+100) || coordInfoVecs.triangles2Nodes_2[i] <= (-INT_MAX+100) || coordInfoVecs.triangles2Nodes_3[i] <= (-INT_MAX+100)){
			continue;
		}
		v1 = coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_1[i]];
		v2 = coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_2[i]];
		v3 = coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_3[i]];
		v4 = (v1+v2+v3)/3.0;
		if (v4 >= max_cell_triangle_height){
			max_cell_triangle_height = v4;
		}
		if (v4 <= min_cell_triangle_height){
			min_cell_triangle_height = v4;
		}
	}
	cell_height = (max_cell_triangle_height - min_cell_triangle_height);
	std::cout<<"Determination of max triangle and min cell height complete."<<std::endl;
	double max_u = 1.5;//1.1;
	double min_u = 0.5;//0.9;
	std::cout<<"max_u = "<<max_u<<std::endl;
	std::cout<<"min_u = "<<min_u<<std::endl;
	// vector<double> u;
	// coordInfoVecs.u.resize(3*coordInfoVecs.num_triangles);
	for (int i = 0; i < coordInfoVecs.num_triangles; i++){
		if (coordInfoVecs.triangles2Nodes_1[i] >= (INT_MAX-100) || coordInfoVecs.triangles2Nodes_2[i] >= (INT_MAX-100) || coordInfoVecs.triangles2Nodes_3[i] >= (INT_MAX-100)){
			continue;
		}
		if (coordInfoVecs.triangles2Nodes_1[i] <= (-INT_MAX+100) || coordInfoVecs.triangles2Nodes_2[i] <= (-INT_MAX+100) || coordInfoVecs.triangles2Nodes_3[i] <= (-INT_MAX+100)){
			continue;
		}
		
		v1 = coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_1[i]];
		v2 = coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_2[i]];
		v3 = coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_3[i]];
		v4 = (v1+v2+v3)/3.0;
		coordInfoVecs.u[i] = min_u + (max_u - min_u)*pow(((v4-min_cell_triangle_height)/cell_height),u_scalingPower); //Try exponential power? C_0*exp(-x/lambda), varying lambda
		// coordInfoVecs.u[i] = max_u - (max_u - min_u)*exp((-(v4-min_cell_triangle_height)/cell_height)/u_scalingPower);
	}
	std::cout<<"Construction of the 'u' vector complete."<<std::endl;
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	generalParams.nodeMass = 1.0;
	int GROWTH_COUNTER = 0;
	int min_num_edge_loop = 1;
	std::cout<<"min_num_edge_loop for edgeswap = "<<min_num_edge_loop<<std::endl;

	// std::random_device rand_dev;
	// // std::mt19937 generator2(rand_dev());
	// std::mt19937 generator_edgeswap(rand_dev());

	double MAX_VOLUME_RATIO = 3.0;
	double MAX_BUD_AREA_RATIO = 100.0;
	int MAX_GROWTH_NUMBER = 1;
	std::cout<<"MAX_GROWTH_NUMBER (# of edge to expand) = "<<MAX_GROWTH_NUMBER<<std::endl;
	int GROWTH_FREQUENCY = 25;//95;//70;//25*3;
	std::cout<<"GROWTH_FREQ (how many times Max_Runtime has to be reached to perform growth"<<GROWTH_FREQUENCY<<std::endl;
	double energy_gradient_threshold = 0.02;//0.01;
	std::cout<<"ENERGY_GRADIENT_THRESHOLD = "<<energy_gradient_threshold<<std::endl;

	generalParams.kT_growth = 1.0;
	generalParams.SCALE_TYPE = 3; 
	// 0:= Gaussian-like weakening
	// 1:= a1*(pow(x,b)) + a2*(1-pow(x,b)) type weakening
	// 2:= pure Gaussian weakening
	// 3:= isotropic
	// 4:= hill equation
	//Note that (3) is used in combination with sigma = INT_MAX;
	std::cout<<"SCALE TYPE = "<<generalParams.SCALE_TYPE<<std::endl;
	std::cout<<"0:= sigmoidal Gaussian-like weakening, 1:= a1*(pow(x,b)) + a2*(1-pow(x,b)) type weakening, 2:= pure Gaussian weakening, 3:= isotropic, 4:= hill equation"<<std::endl;
	generalParams.scaling_pow = 2.0;
	std::cout<<"scaling_pow (this is for SCALE_TYPE = 1 case) = "<<generalParams.scaling_pow<<std::endl;
	generalParams.gausssigma = 0.1;
	std::cout<<"gausssigma (this is for the SCALE_TYPE = 0 case) = "<<generalParams.gausssigma<<std::endl;
	//coordInfoVecs.scaling_per_edge.
	//generalParams.hilleqnconst = 0.9;
	//generalParams.hilleqnpow = 40.0;
	std::vector<int> nodes_in_growth;
	std::vector<int> triangles_in_growth;
	std::vector<int> edges_in_growth;
	double dtb; //dtb := distance to boundary
	double dtb_max; //dtb_max := the max distance used to calculate the distance ratio in the Hill equation.
	double sigma = 0.0;//INT_MAX; //if this is set to be INT_MAX then we assume isotropic weakening.
	double sigma_true = sqrt(0.5); //This is the variance used to calculate the scaling of the wall weakening.
	std::cout<<"initial sigma (for gradient distribution variance), based on initial distribution of Cdc42, if using true gaussian weakening = "<<sigma<<std::endl;
	std::cout<<"If sigma = INT_MAX, then we have isotropic weakening scenario"<<std::endl;
	std::cout<<"true sigma (for gaussian-related distribution variance) = "<<sigma_true<<std::endl;

	generalParams.insertion_energy_cost = -log(0.0025);
	std::cout<<"GROWTH: material insertion energy cost (dependent on local chemical concentration) = "<<generalParams.insertion_energy_cost<<std::endl;
	generalParams.strain_threshold = 0.05;//0.01;
	std::cout<<"GROWTH: critical strain threshold used for insertion probability calculation = "<<generalParams.strain_threshold<<std::endl;

	generalParams.growth_energy_scaling = 1.0;//0.01375;
	std::cout<<"GROWTH ENERGY SCALING FOR MATERIAL INSERTION PROBABILITY = "<<generalParams.growth_energy_scaling<<std::endl;
	generalParams.safeguardthreshold = 9;
	std::cout<<"NEIGHBOR SAFE GUARD THRESHOLD = "<<generalParams.safeguardthreshold<<std::endl;
	//safeguardthreshold is the maximum number of neighboring nodes a node can have.

	//////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////
	////////////////////////// PARAMETER SETTINGS ////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////

	int TOTAL_GROWTH_COUNTER = 0;
	int TOTAL_GROWTH_ATTEMPT = 0;
	double max_runTimePerRelaxation = 50;
	std::cout<<"Number of sim steps per relaxation event : "<<max_runTimePerRelaxation<<std::endl;

	double Max_Runtime = generalParams.dt*max_runTimePerRelaxation;//50.0;//50.0;
	double minimal_run_time_ratio = 1.0;
	double Max_RunStep = Max_Runtime/generalParams.dt;
	std::cout<<"Max runtime = "<<Max_Runtime<<std::endl;
	std::cout<<"Max runstep = "<<Max_RunStep<<std::endl;
	bool runSim = true;
	int num_edge_loop;
	double initial_kT;
	initial_kT = generalParams.kT;//This is for the acceptance of change after looping through every edge within proximity.
	double SAMPLE_SIZE = 0.05;//0.025;
	std::cout<<"Sample ratio: "<<SAMPLE_SIZE<<std::endl;
	std::cout<<"If the Sample raio is 0, it means we have chosen a fixed number of attempt throughout the simulation"<<std::endl;
	//This determines the number of edges to test for bondflip remeshing

	// auto utilities_ptr = std::make_shared<Edgeswap>(coordInfoVecs, generalParams);
	int RECORD_TIME = 1;//round(Max_RunStep/2);
	std::cout<<"Record frequency = "<<RECORD_TIME<<std::endl;
	//int GROWTH_TIME = 1;
	//std::cout<<"Growth frequency = "<<GROWTH_TIME<<std::endl;
	int translate_frequency = 10;
	std::cout<<"recentering of the model cell frequency = "<<translate_frequency<<std::endl;
	//translate_frequency determines the frequency for the mesh to re-center and perform dynamical remeshing
	int NUMBER_OF_GROWTH_EVENT = 2000;//1000;//1000*2;
	// std::cout<<"Number of maximally allowed growth event = "<<NUMBER_OF_GROWTH_EVENT<<" which used to terminate the simulation if not enough growth is encountered for a prolonged simulation."<<std::endl;
	int NUMBER_OF_TARGETED_GROWTH_EVENT = 1;
	int NKBT = GROWTH_FREQUENCY*NUMBER_OF_GROWTH_EVENT*2;//GROWTH_FREQUENCY*NUMBER_OF_GROWTH_EVENT;//10000;//7500;
	std::cout<<"Number of edge-swap per kBT value (or total number of edge-swap if kBT is fixed), NKBT = "<<NKBT<<std::endl;
	int GROWTH_FREQUENCY_SCALE = 4;
	std::cout<<"GROWTH FREQ SCALE: decides how many growth event must be checked before recording the result = "<<GROWTH_FREQUENCY_SCALE<<std::endl;
	double min_kT = -0.1;//0.21;
	std::cout<<"min kT for simulation termination = "<<min_kT<<std::endl;
	int WHEN = 0;
	
	std::cout<<"Total number of simulation steps: "<<NKBT*max_runTimePerRelaxation<<std::endl;
	double old_total_energy = 0.0;
	double new_total_energy = 0.0;
	double energy_gradient = 0.0;
	double energy_rep = 0.0;
	int Num_of_step_run = 0;
	// auto build_ptr = weak_bld_ptr.lock();//upgrade weak builder to access host variables.
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
	
	//double COMPRESS = 2.0227;
	// double COMPRESS2 = -2.0227;

	/////////////////////////////////////////////////////////////////
	/////////////////////// MEMBRANE RELATED ////////////////////////
	/////////////////////////////////////////////////////////////////
	
	std::vector<double> nodenormal_1(generalParams.maxNodeCount, 0.0);
	std::vector<double> nodenormal_2(generalParams.maxNodeCount, 0.0);
	std::vector<double> nodenormal_3(generalParams.maxNodeCount, 0.0);
	int reduce_counter = 0;

	double VOLUME_FACTOR = MAX_VOLUME_RATIO;//1.6;//2.25;
	//VOLUME_FACTOR determines the target volume which equals to VOLUME_FACTOR*initial_volume.
	//double tip_depth = 0.5;
	//tip_depth is currently unused.

	double LINE_TENSION_THRESHOLD = -10000.0;
	std::cout<<"LINE TENSION THRESHOLD for activation of line tension = "<<LINE_TENSION_THRESHOLD<<std::endl;
	double VOLUME_THRESHOLD = 0.0;
	std::cout<<"VOLUME THRESHOLD for activation of weakened membrane = "<<VOLUME_THRESHOLD<<std::endl;
	
	double weakened = 1.90;//6.0;
	//weakened determines the minimum height of the z-coordinate of the membrane node to be considered in the area of weakened mechanical properties.
	//double tip_base = 6.0;
	//tip_base currently unused.

	// double EXPAN_THRESHOLD = 0.0;
	// double EXPAN_THRESHOLD_weak = 0.0;//1.75;
	// std::cout<<"EXPANSION THRESHOLD = "<<EXPAN_THRESHOLD<<std::endl;
	// int RULES_OF_EXPAN = 1;	//EXPAN_THRESHOLD is the yielding ratio where a pair of triangles will perform expansion.
	
	// std::cout<<"EXPANSION RULE = "<<RULES_OF_EXPAN<<std::endl;
	// //EXPAN_THRESHOLD_weak is the secondary yielding ratio.
	// //RULES_OF_EXPAN controls how the EXPAN_THRESHOLD is applied:
	// // 1:= Both trianglular areas must exceed the threshold value.
	// // 2:= If one trianglular area exceeds the treshold value while the other exceeds the secondary threshold value.
	// // 3:= If the combined area of the two triangles exceed 2*EXPAN_THRESHOLD.
	// // 4:= If a selected edges exceed the threshold value, split the triangles associated with the edge.

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
	std::vector<int> VectorShuffleForEdgeswapLoop;

	double max_height = coordInfoVecs.nodeLocZ[35];
	double min_height = coordInfoVecs.nodeLocZ[38];
	int max_height_index = 35;
	/*double max_height = -10000.0;
	int max_height_index = -1;
	std::vector<int> Stiffness_gradient();
    for (int k = 0; k < generalParams.maxNodeCount; k++){
        if (coordInfoVecs. nodeLocZ[k] >= max_height){
			max_height = coordInfoVecs. nodeLocZ[k];
			max_height_index = k;
            }
	}*/
	//Max and min height of the membrane nodes, these have to be changed if the mesh used is changed.

	generalParams.Rmin = 0.3012;//0.15;
	//Equilibrium length of an edge of the triangle.
	//generalParams.Rmin_growth = 0.329;
	generalParams.abs_Rmin = generalParams.Rmin;//0.15;
	std::cout<<"abs_Rmin = "<<generalParams.abs_Rmin<<std::endl;
	//Equilibrium distance between membrane node for volume exclusion.
	areaTriangleInfoVecs.initial_area = 0.039;//2835;//0.009808;//0.039;//0.03927344;//0.009817;
	std::cout<<"equilibrium triangular area = "<<areaTriangleInfoVecs.initial_area<<std::endl;
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

	linearSpringInfoVecs.spring_constant_rep1 = 0.01;//0.023;
	linearSpringInfoVecs.spring_constant_rep2 = 9.0;//5.0;
	std::cout<<"Membrane volume exclusion Morse D = "<<linearSpringInfoVecs.spring_constant_rep1<<std::endl;
	std::cout<<"Membrane volume exclusion Morse a = "<<linearSpringInfoVecs.spring_constant_rep2<<std::endl;
	//The coefficient used for non-neighboring membrane node volume exclusion.
	//rep1 is the "D" and rep2 is the "alpha" in the standard form of Morse potential.

	generalParams.volume_spring_constant = 0.2;//(1.0/3.0)*areaTriangleInfoVecs.initial_area*1.0;
	std::cout<<"spring constant for surface normal expansion (pressure within the cell) = "<<generalParams.volume_spring_constant<<std::endl;
	generalParams.line_tension_constant = 0.0;//250.0;
	std::cout<<"spring constant for the septin ring (before budding) = "<<generalParams.line_tension_constant<<std::endl;
	generalParams.length_scale = 1.0;//0.85;//0.1577;//1.0*generalParams.Rmin;// 0.8333;
	//std::cout<<"equilibrium length of each segment of the septin ring = "<<generalParams.length_scale<<std::endl;

	// bendingTriangleInfoVecs.spring_constant = bendingTriangleInfoVecs.spring_constant*(2.0/sqrt(3));

	double scale_linear = linearSpringInfoVecs.spring_constant*1.0;//0.25;//25.0/2.5;//75.0/15.0;
	double scale_bend = bendingTriangleInfoVecs.spring_constant*1.0;//0.05;//10.0/1.0;//75.0/7.5;
	double scale_area = areaTriangleInfoVecs.spring_constant*1.0;//0.25;//50.0/5.0;//75.0/15.0;
	std::cout<<"weakened region linear (before budding) = "<<scale_linear<<std::endl;
	std::cout<<"weakened region bend (before budding) = "<<scale_bend<<std::endl;
	std::cout<<"weakened region area (before budding) = "<<scale_area<<std::endl;
	//linearSpringInfoVecs.spring_constant_weak = linearSpringInfoVecs.spring_constant/scale_linear;
	//bendingTriangleInfoVecs.spring_constant_weak = bendingTriangleInfoVecs.spring_constant/scale_bend;
	//areaTriangleInfoVecs.spring_constant_weak = areaTriangleInfoVecs.spring_constant/scale_area;
	linearSpringInfoVecs.spring_constant_weak = scale_linear;
	bendingTriangleInfoVecs.spring_constant_weak = scale_bend;
	areaTriangleInfoVecs.spring_constant_weak = scale_area;
	//Scaling of the weakend mechanical properties.

	bendingTriangleInfoVecs.initial_angle = 0.087165870975460;//0.087249;//0.04335;
	bendingTriangleInfoVecs.initial_angle_raft = 0.087165870975460;//0.087249;//0.04335;
	bendingTriangleInfoVecs.initial_angle_coat = 0.087165870975460;//0.087249;//0.04335;
	std::cout<<"equilibrium bending angle of the membrane = "<<bendingTriangleInfoVecs.initial_angle<<std::endl;
	//raft and coat are current unused due to the assumption of uniform preferred curvature.

	bendingTriangleInfoVecs.initial_angle_bud = bendingTriangleInfoVecs.initial_angle;
	std::cout<<"equilibrium bending angle of the bud = "<<bendingTriangleInfoVecs.initial_angle_bud<<std::endl;
	
	// bendingTriangleInfoVecs.spring_constant_raft = 0.0;//bendingTriangleInfoVecs.spring_constant;
	// bendingTriangleInfoVecs.spring_constant_coat = 0.0;//bendingTriangleInfoVecs.spring_constant;
	// bendingTriangleInfoVecs.spring_constant = bendingTriangleInfoVecs.spring_constant*(2.0/sqrt(3));
	// bendingTriangleInfoVecs.spring_constant_raft = bendingTriangleInfoVecs.spring_constant_raft*(2.0/sqrt(3));
	// bendingTriangleInfoVecs.spring_constant_coat = bendingTriangleInfoVecs.spring_constant_coat*(2.0/sqrt(3));
	// std::cout<<"Effective bending coefficient is calculated by multiplying 2/sqrt(3)"<<std::endl;
	// std::cout<<"effective bending coefficient of the membrane = "<<bendingTriangleInfoVecs.spring_constant<<std::endl;
	// std::cout<<"effective bending coefficient of the membrane raft = "<<bendingTriangleInfoVecs.spring_constant_raft<<std::endl;
	// std::cout<<"effective bending coefficient of the membrane coat = "<<bendingTriangleInfoVecs.spring_constant_coat<<std::endl;

	std::vector<int> pull_nodes_up;// = {35,    76,    79,   111,   113,   151,   153,   360,   361,   362,   363,   364,   365,   505,   506,   515,   516,   593,   632};//{35, 360,   361,   362,   363,   364,   365};
	std::vector<int> pull_nodes_down;// = {38,    86,    89,   121,   123,   144,   146,   378,   379,   380,   381,   382,   383,   535,   536,   545,   546,   602,   626};//{38, 378,   379,   380,   381,   382,   383};
	std::vector<int> push_nodes_down;
	std::vector<int> push_nodes_up;
	// for (int i = 0; i < generalParams.maxNodeCount; i++){
	// 	if (coordInfoVecs.nodeLocZ[i] >= 1.43026488631){
	// 		pull_nodes_up.push_back(i);
	// 	}
	// 	if (coordInfoVecs.nodeLocZ[i] <= -1.43026488631){
	// 		pull_nodes_down.push_back(i);
	// 	}
	// }

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
	
	std::vector<int> row2 = {35 ,   76 ,   79 ,  111 ,  113 ,  151 ,  153 ,  360 ,  361 ,  362 ,  363 ,  364 ,  365 ,  505 ,  506 ,  515 ,  516 ,  593 ,  632};
	// std::vector<int> row2 = {35,76,79,111,113,151,153,360,361,362,363,364,365,505,506,515,516,593,632,840,841,842,
	//    843,844,845,1087,1090,1091,1105,1108,1109,1297,1299,1301,1309,1311,1313,1537,1539,1541,1549,1551,1553,2196,
	//   2197,2198,2199,2200,2201,2202,2203, 2204, 2205,2206,2207, 2208,2209,2210,2211,2212,2213};
	//std::vector<int> nodes_to_center;
	//generalParams.nodes_in_upperhem.resize(generalParams.maxNodeCount,-1);

	for (int i = 0; i < generalParams.maxNodeCount; i++){
		generalParams.nodes_in_upperhem[i] = -1;
		// generalParams.nodes_in_upperhem[i] = 1;
	}

	for (int i = 0; i < row2.size(); i++){
		generalParams.nodes_in_upperhem[row2[i]] = 1;
	//	std::cout<<"nodes "<<i<<" "<<generalParams.nodes_in_upperhem[i]<<std::endl;		
	}
	// for (int i = 0; i < generalParams.maxNodeCount; i++){
	// 	if (coordInfoVecs.nodeLocZ[i] > (generalParams.centerZ + weakened)){
	// 		generalParams.nodes_in_upperhem[i] = 1;
	// 	}
	// 	else{
	// 		generalParams.nodes_in_upperhem[i] = -1;
	// 	}
	// //	std::cout<<"nodes "<<i<<" "<<generalParams.nodes_in_upperhem[i]<<std::endl;		
	// }

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

	//generalParams.triangles_in_upperhem.resize(coordInfoVecs.num_triangles);
	for (int i = 0; i < coordInfoVecs.num_triangles; i++){
		if (coordInfoVecs.triangles2Nodes_1[i] >= (INT_MAX-1000) || coordInfoVecs.triangles2Nodes_1[i] < 0){
			generalParams.triangles_in_upperhem[i] = -1;
			continue;
		}
		else if (coordInfoVecs.triangles2Nodes_2[i] >= (INT_MAX-1000) || coordInfoVecs.triangles2Nodes_2[i] < 0){
			generalParams.triangles_in_upperhem[i] = -1;
			continue;
		}
		else if (coordInfoVecs.triangles2Nodes_3[i] >= (INT_MAX-1000) || coordInfoVecs.triangles2Nodes_3[i] < 0){
			generalParams.triangles_in_upperhem[i] = -1;
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
		//else if ((aaa+bbb+ccc)==1){
		//	generalParams.triangles_in_upperhem[i] = 0;
			//triangles_in_upperhem.push_back(i);
		//}
		else{
			generalParams.triangles_in_upperhem[i] = -1;
		}
	//	std::cout<<"triangle "<<i<<" "<<generalParams.triangles_in_upperhem[i]<<std::endl;		
	}

	//std::vector<int> edges_in_upperhem;
//	generalParams.edges_in_upperhem.resize(coordInfoVecs.num_edges);
	int edges_in_upperhem_COUNT = 0;
	for (int i = 0; i < coordInfoVecs.num_edges; i++){
		if (coordInfoVecs.edges2Triangles_1[i] >= (INT_MAX-1000) || coordInfoVecs.edges2Triangles_1[i] < 0){
			generalParams.edges_in_upperhem[i] = -1;
			generalParams.edges_in_upperhem_list[i] = -INT_MAX;
			// std::cout<<"Filter1, edge = "<<i<<std::endl;
			continue;
		}
		else if (coordInfoVecs.edges2Triangles_2[i] >= (INT_MAX-1000) || coordInfoVecs.edges2Triangles_2[i] < 0){
			generalParams.edges_in_upperhem[i] = -1;
			generalParams.edges_in_upperhem_list[i] = -INT_MAX;
			// std::cout<<"Filter2, edge = "<<i<<std::endl;
			continue;
		}
		else{
			// std::cout<<"Filter 3, edge = "<<i<<std::endl;
			int aaa = generalParams.triangles_in_upperhem[coordInfoVecs.edges2Triangles_1[i]];//generalParams.nodes_in_upperhem[coordInfoVecs.edges2Nodes_1[i]];
			int bbb = generalParams.triangles_in_upperhem[coordInfoVecs.edges2Triangles_2[i]];//generalParams.nodes_in_upperhem[coordInfoVecs.edges2Nodes_2[i]];
			if (aaa == 1 && bbb == 1){
				generalParams.edges_in_upperhem[i] = 1;
				//generalParams.edges_in_upperhem_list.push_back(i);
				generalParams.edges_in_upperhem_list[i] = i;
				edges_in_upperhem_COUNT += 1;
			}
			else if (aaa == 1 || bbb == 1){
				generalParams.edges_in_upperhem[i] = 1;
				generalParams.edges_in_upperhem_list[i] = -INT_MAX;
				edges_in_upperhem_COUNT += 1;
			}
			else{
				generalParams.edges_in_upperhem[i] = -1;
				generalParams.edges_in_upperhem_list[i] = -INT_MAX;
			}
		}
		// std::cout<<edges_in_upperhem_COUNT<<std::endl;
		
	}
	std::cout<<"INITIAL EDGES IN UPPERHEM = "<<edges_in_upperhem_COUNT<<std::endl;

	int COUNTING_EDGE = 0;
	for (int y = 0; y < coordInfoVecs.num_edges; y++){
		if (generalParams.edges_in_upperhem_list[y] >= 0){
			COUNTING_EDGE += 1;
		}
		generalParams.edges_in_upperhem_list_length = COUNTING_EDGE;
	}
	

	//Find the boundary of the nodes_in_upperhem region
	//generalParams.boundaries_in_upperhem.resize(coordInfoVecs.num_edges);
	std::vector<int> boundary_node_list;
	std::vector<int> boundary_edge_list;
	for (int i = 0; i < coordInfoVecs.num_edges; i++){
		double T1 = coordInfoVecs.edges2Triangles_1[i];
		double T2 = coordInfoVecs.edges2Triangles_2[i];
		if (T1 >= (INT_MAX - 1000) || T1 < 0 || T2 >= (INT_MAX-1000) || T2 < 0){
			continue;
		}
		if (generalParams.triangles_in_upperhem[T1] == 1 && generalParams.triangles_in_upperhem[T2] != 1){
			generalParams.boundaries_in_upperhem[i] = 1;
			//std::cout<<generalParams.boundaries_in_upperhem[i]<<std::endl;
		//	generalParams.triangles_in_upperhem[T1] = 0;
		//	generalParams.triangles_in_upperhem[T2] = 0;
			double bdry_node1 = coordInfoVecs.edges2Nodes_1[i];
			double bdry_node2 = coordInfoVecs.edges2Nodes_2[i];
			boundary_node_list.push_back(bdry_node1);
			boundary_node_list.push_back(bdry_node2);
			boundary_edge_list.push_back(i);
			//generalParams.nodes_in_upperhem[bdry_node1] = 0;
			//generalParams.nodes_in_upperhem[bdry_node2] = 0;
			// coordInfoVecs.isNodeFixed[bdry_node1] = true;
			// coordInfoVecs.isNodeFixed[bdry_node2] = true;
		}
		else if (generalParams.triangles_in_upperhem[T1] != 1 && generalParams.triangles_in_upperhem[T2] == 1){
			generalParams.boundaries_in_upperhem[i] = 1;
			//std::cout<<generalParams.boundaries_in_upperhem[i]<<std::endl;
		//	generalParams.triangles_in_upperhem[T1] = 0;
		//	generalParams.triangles_in_upperhem[T2] = 0;
			double bdry_node1 = coordInfoVecs.edges2Nodes_1[i];
			double bdry_node2 = coordInfoVecs.edges2Nodes_2[i];
			boundary_node_list.push_back(bdry_node1);
			boundary_node_list.push_back(bdry_node2);
			boundary_edge_list.push_back(i);
			//generalParams.nodes_in_upperhem[bdry_node1] = 0;
			//generalParams.nodes_in_upperhem[bdry_node2] = 0;
			// coordInfoVecs.isNodeFixed[bdry_node1] = true;
			// coordInfoVecs.isNodeFixed[bdry_node2] = true;
		}
		else {
			generalParams.boundaries_in_upperhem[i] = -1;
			//std::cout<<generalParams.boundaries_in_upperhem[i]<<std::endl;
		}
	}
	std::cout<<"size of boundary_node_list (this is double-counted) = "<<boundary_node_list.size()<<std::endl;
	//generalParams.eq_total_boundary_length = generalParams.boundaries_in_upperhem.size()*generalParams.Rmin;

	/*for (int i = 0; i < coordInfoVecs.num_edges; i++){
		int aaa = coordInfoVecs.edges2Nodes_1[i];
		int bbb = coordInfoVecs.edges2Nodes_2[i];
		if (aaa == 1 && bbb == 1){
			generalParams.edges_in_upperhem[i] = 1;
			generalParams.edges_in_upperhem_list.push_back(i);
		}
		else if (aaa == 1 || bbb == 1){
			generalParams.edges_in_upperhem[i] = 0;
		}
		else{
			generalParams.edges_in_upperhem[i] = -1;
		}
		
	}*/
	
	

	int true_num_edges_in_upperhem = 0;
	for (int i = 0; i < coordInfoVecs.num_edges; i++){
		if (generalParams.edges_in_upperhem_list[i] != INT_MAX && generalParams.edges_in_upperhem_list[i] >= 0){
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


	//std::cout<<"ERROR HERE?"<<std::endl;
	ComputeVolume(
		generalParams,
		coordInfoVecs,
		linearSpringInfoVecs,
		ljInfoVecs
	);
	//std::cout<<"ERROR HERE 2?"<<std::endl;
	double initial_volume;
	// initial_volume = generalParams.true_current_total_volume;
	// generalParams.eq_total_volume = generalParams.true_current_total_volume*VOLUME_FACTOR;//This is for setting different equilibrium volume to mimic growth or shirnkage.
	// std::cout<<"true current total volume = "<<generalParams.true_current_total_volume<<std::endl;
	// std::cout<<"eq total volume = "<<generalParams.eq_total_volume<<std::endl;

	//////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////// START OF ACTUAL SIMULATION /////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////

	/* Build the initial gradient weakend scale */
	dtb = 0.0;//dtb := distance to boundary
	generalParams.septin_ring_z = 0.0;
	generalParams.boundary_z = 0.0;
	//for (int k = 0; k < boundary_edge_list.size(); k++){
	for (int k = 0; k < boundary_node_list.size(); k++){
		double n1 = boundary_node_list[k];//coordInfoVecs.edges2Nodes_1[boundary_edge_list[k]];
		//double n2 = coordInfoVecs.edges2Nodes_2[boundary_edge_list[k]];
		//double cent_of_edge_x = (coordInfoVecs.nodeLocX[n1] + coordInfoVecs.nodeLocX[n2])/2.0;
		//double cent_of_edge_y = (coordInfoVecs.nodeLocY[n1] + coordInfoVecs.nodeLocY[n2])/2.0;
		//double cent_of_edge_z = (coordInfoVecs.nodeLocZ[n1] + coordInfoVecs.nodeLocZ[n2])/2.0;
		double dist_x = coordInfoVecs.nodeLocX[max_height_index] - coordInfoVecs.nodeLocX[n1];//cent_of_edge_x;
		double dist_y = coordInfoVecs.nodeLocY[max_height_index] - coordInfoVecs.nodeLocY[n1];//cent_of_edge_y;
		double dist_z = coordInfoVecs.nodeLocZ[max_height_index] - coordInfoVecs.nodeLocZ[n1];//cent_of_edge_z;
		// double temp_dist = sqrt((coordInfoVecs.nodeLocX[max_height_index] - coordInfoVecs.nodeLocX[n1])*(coordInfoVecs.nodeLocX[max_height_index] - coordInfoVecs.nodeLocX[n1]) +
		// (coordInfoVecs.nodeLocY[max_height_index] - coordInfoVecs.nodeLocY[n1])*(coordInfoVecs.nodeLocY[max_height_index] - coordInfoVecs.nodeLocY[n1]) +
		// 	(coordInfoVecs.nodeLocZ[max_height_index] - coordInfoVecs.nodeLocZ[n1])*(coordInfoVecs.nodeLocZ[max_height_index] - coordInfoVecs.nodeLocZ[n1]));
		double temp_dist = sqrt(dist_x*dist_x + dist_y*dist_y + dist_z*dist_z);
		generalParams.septin_ring_z += coordInfoVecs.nodeLocZ[n1];
		if (temp_dist >= dtb){
			dtb = temp_dist;
			/* "dtb" will be used to identify where the septin ring is located, and used to determine the Hill coefficient*/
		}
	}
	std::cout<<"dtb = "<<dtb<<std::endl;
	//generalParams.septin_ring_z = generalParams.septin_ring_z/boundary_node_list.size();
	//generalParams.boundary_z = generalParams.septin_ring_z - generalParams.Rmin;
	/* dtb will be only calculated once so we can effectively keep the Hill eqn curve consistent with only horizontal shift */
	dtb_max = dtb + (generalParams.Rmin);
	
	std::cout<<"initial distance between cell tip and the boundary of weakened area = "<<dtb<<std::endl;
	std::cout<<"Notice that here, the distance from the tip to the boundary is slightly extended by half of the equilibrium length of an edge"<<std::endl;
	//std::cout<<"If this message is present, we are forcing a fixed portion of the bud tip to be occupied by the max concentration"<<std::endl;
	//generalParams.hilleqnconst = (dtb + generalParams.Rmin/4.0)/dtb_max;
	generalParams.hilleqnconst = dtb/dtb_max;
	generalParams.hilleqnpow = 70.0;
	// std::cout<<"hill equation constant K = "<<generalParams.hilleqnconst<<std::endl;
	// std::cout<<"hill (equation) coefficient = "<<generalParams.hilleqnpow<<std::endl;
	// std::cout<<"NOTE: IN THIS SIMULATION, THE LOCATION WHERE 50% WEAKENING IS EXPERIENCED IS LOCATED SLIGHTLY AWAY FROM THE SEPTIN RING, "<<std::endl;
	// std::cout<<"THIS IS DUE TO THE FACT THAT IN ISOTROPIC CASE, SEPTIN RING LOCATION MUST BE SUFFICIENTLY WEAKENED TO INDUCE BUDDING"<<std::endl;
	// std::cout<<" "<<std::endl;
	// std::cout<<" "<<std::endl;
	// std::cout<<" "<<std::endl;
	// std::cout<<" "<<std::endl;
	// std::cout<<" "<<std::endl;
	// std::cout<<" "<<std::endl;
	// std::cout<<" "<<std::endl;
	// std::cout<<" "<<std::endl;


	utilities_ptr->transferDtoH(generalParams, coordInfoVecs, build_ptr->hostSetInfoVecs);
	std::cout<<"HA?"<<std::endl;
	utilities_ptr->gradient_weakening_update_host_vecs(sigma,
		//max_height_index,
		coordInfoVecs.nodeLocX[max_height_index],
		coordInfoVecs.nodeLocY[max_height_index],
		coordInfoVecs.nodeLocZ[max_height_index],
		dtb,
		dtb_max,
		generalParams,
		coordInfoVecs,
		build_ptr->hostSetInfoVecs);
		std::cout<<"HAHA?"<<std::endl;
	for (int u = 0; u < generalParams.maxNodeCount; u++){
		// std::cout<<"u = "<<u<<std::endl;
		int BETA = utilities_ptr->nodes2Triangles_host_vecs(
			u,
			build_ptr->hostSetInfoVecs,
			coordInfoVecs,
			generalParams,
			auxVecs);
	}
	utilities_ptr->transferHtoD(generalParams, coordInfoVecs, build_ptr->hostSetInfoVecs);
	std::cout<<"STARTING THE ACTUAL SIMULATION"<<std::endl;
	while (runSim == true){
		
		double current_time = 0.0;

		int translate_counter = 0;
		
        while (current_time < relax_max_steps_before_growth_and_edgeswap*(Max_Runtime)){
            translate_counter += 1;
            Solve_Forces();
            double beta;
                
            AdvancePositions(
                coordInfoVecs,
                generalParams,
                domainParams);
                        
            new_total_energy = linearSpringInfoVecs.linear_spring_energy + 
                areaTriangleInfoVecs.area_triangle_energy + 
                bendingTriangleInfoVecs.bending_triangle_energy;// + 
                0.5*energy_rep;// + 
                //ljInfoVecs.lj_energy_M +
                //ljInfoVecs.lj_energy_LJ +
                //generalParams.volume_energy;

            // energy_gradient = sqrt((new_total_energy - old_total_energy)*(new_total_energy - old_total_energy))/old_total_energy;
            // if (current_time >= Max_Runtime*minimal_run_time_ratio && (energy_gradient/generalParams.dt) < energy_gradient_threshold){
            //     break;
            //     }
            old_total_energy = new_total_energy;
            current_time+=generalParams.dt;
        }

		std::cout<<"Time used for 'steady state' initial condition before growth and edge swaps = "<<current_time<<std::endl;
		std::cout<<"current total energy (before growth and edge swaps) = "<<new_total_energy<<std::endl;
		std::cout<<"LINEAR ENERGY = "<<linearSpringInfoVecs.linear_spring_energy<<std::endl;
		std::cout<<"BEND ENERGY = "<<bendingTriangleInfoVecs.bending_triangle_energy<<std::endl;
		std::cout<<"AREA ENERGY = "<<areaTriangleInfoVecs.area_triangle_energy<<std::endl;
		//std::cout<<"REPULSION ENERGY = "<<energy_rep<<std::endl;
		std::cout<<"VOLUME ENERGY = "<<generalParams.volume_energy<<std::endl;
		std::cout<<"true_current_total_volume (before growth and edge swaps) = "<<generalParams.true_current_total_volume<<std::endl;
		// std::cout<<"eq_total_volume = "<<generalParams.eq_total_volume<<std::endl;
		std::cout<<"current KBT = "<<generalParams.kT<<std::endl;
		if (isnan(new_total_energy)==1){
			std::cout<<"Nan or Inf position update !!!!"<<std::endl;
			runSim = false;
			break;
		}
		double current_bud_area = 0.0;
		for (int k = 0; k < coordInfoVecs.num_triangles; k++){
			if (coordInfoVecs.triangles2Nodes_1[k] >= (INT_MAX - 1000.0) || coordInfoVecs.triangles2Nodes_1[k] <= (-INT_MAX + 1000.0) ||
				coordInfoVecs.triangles2Nodes_2[k] >= (INT_MAX - 1000.0) || coordInfoVecs.triangles2Nodes_2[k] <= (-INT_MAX + 1000.0) ||
				coordInfoVecs.triangles2Nodes_3[k] >= (INT_MAX - 1000.0) || coordInfoVecs.triangles2Nodes_3[k] <= (-INT_MAX + 1000.0)){
						continue;
					}
			else{
				if (generalParams.triangles_in_upperhem[k] == 1){
					double r1x = coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_1[k]];
					double r1y = coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_1[k]];
					double r1z = coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_1[k]];
					double r2x = coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_2[k]];
					double r2y = coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_2[k]];
					double r2z = coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_2[k]];
					double r3x = coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_3[k]];
					double r3y = coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_3[k]];
					double r3z = coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_3[k]];
					double norm_r1r2 = sqrt((r2x-r1x)*(r2x-r1x) + (r2y-r1y)*(r2y-r1y) + (r2z-r1z)*(r2z-r1z));
					double norm_r2r3 = sqrt((r3x-r2x)*(r3x-r2x) + (r3y-r2y)*(r3y-r2y) + (r3z-r2z)*(r3z-r2z));
					double norm_r3r1 = sqrt((r3x-r1x)*(r3x-r1x) + (r3y-r1y)*(r3y-r1y) + (r3z-r1z)*(r3z-r1z));
					double s = (norm_r1r2 + norm_r2r3 + norm_r3r1)/2.0;
					double area = sqrt(s*(s-norm_r1r2)*(s-norm_r2r3)*(s-norm_r3r1));
					current_bud_area += area;
				}
			}
		}
		double Initial_Bud_Area = current_bud_area;
		std::cout<<"Initial bud surface area (before growth and edge swaps) = "<<Initial_Bud_Area<<std::endl;

		generalParams.volume_spring_constant = 0.2;//(1.0/3.0)*areaTriangleInfoVecs.initial_area*1.0;
		std::cout<<"spring constant for surface normal expansion (pressure within the cell) = "<<generalParams.volume_spring_constant<<std::endl;
		generalParams.line_tension_constant = 50.0;//250.0;
		std::cout<<"spring constant for the septin ring = "<<generalParams.line_tension_constant<<std::endl;
		generalParams.length_scale = 1.0;//0.85;//0.1577;//1.0*generalParams.Rmin;// 0.8333;
		//std::cout<<"equilibrium length of each segment of the septin ring = "<<generalParams.length_scale<<std::endl;

		double scale_linear = linearSpringInfoVecs.spring_constant*0.75;//0.25;//25.0/2.5;//75.0/15.0;
		double scale_bend = bendingTriangleInfoVecs.spring_constant*0.135;//0.05;//10.0/1.0;//75.0/7.5;
		double scale_area = areaTriangleInfoVecs.spring_constant*0.75;//0.25;//50.0/5.0;//75.0/15.0;
		std::cout<<"weakened region linear = "<<scale_linear<<std::endl;
		std::cout<<"weakened region bend = "<<scale_bend<<std::endl;
		std::cout<<"weakened region area = "<<scale_area<<std::endl;
		//linearSpringInfoVecs.spring_constant_weak = linearSpringInfoVecs.spring_constant/scale_linear;
		//bendingTriangleInfoVecs.spring_constant_weak = bendingTriangleInfoVecs.spring_constant/scale_bend;
		//areaTriangleInfoVecs.spring_constant_weak = areaTriangleInfoVecs.spring_constant/scale_area;
		linearSpringInfoVecs.spring_constant_weak = scale_linear;
		bendingTriangleInfoVecs.spring_constant_weak = scale_bend;
		areaTriangleInfoVecs.spring_constant_weak = scale_area;
		//Scaling of the weakend mechanical properties.
		initial_volume = generalParams.true_current_total_volume;
		generalParams.eq_total_volume = generalParams.true_current_total_volume*VOLUME_FACTOR;//This is for setting different equilibrium volume to mimic growth or shirnkage.
		std::cout<<"true current total volume = "<<generalParams.true_current_total_volume<<std::endl;
		std::cout<<"eq total volume = "<<generalParams.eq_total_volume<<std::endl;
	
		storage->print_VTK_File();
		storage->storeVariables();

		utilities_ptr->LDG_Surface_Diffusion_Initialize(
			coordInfoVecs,
			generalParams,
			build_ptr->hostSetInfoVecs,
			auxVecs,
			rbc,
			n_rbc,
			mem_prealloc
		);

		int edgeswap_iteration = 0;
		num_edge_loop = 0;//round(true_num_edges_in_upperhem*SAMPLE_SIZE);	

		int LINE_TENSION_START = 0;
		
		bool WEAKENED_START = false;
		bool EDGESWAP_ALGORITHM_TRIGGERED;
		bool needToRebuildDiffStructAfterEdgeSwap = false;
		int number_of_simulation_step = 0;
 		while (initial_kT > 0){
			if (edgeswap_iteration >= NKBT){
				runSim = false;
				initial_kT = -1;
				break;
			}
			////////////////////NOW RELAX THE ATTEMPTED EDGESWAP//////////////////////
				current_time = 0.0;
				translate_counter = 0;
				double VOLUME_RATIO = generalParams.true_current_total_volume/generalParams.eq_total_volume;
			
			if (generalParams.true_current_total_volume/initial_volume >= LINE_TENSION_THRESHOLD && edgeswap_iteration == 0){
			// 	if (LINE_TENSION_START < 1){
					double DIST = 0.0;
					double COUNT = 0.0;
					for (int t = 0; t < coordInfoVecs.num_edges; t++){
						if (generalParams.boundaries_in_upperhem[t] == 1){
							COUNT += 1.0;
							int node1 = coordInfoVecs.edges2Nodes_1[t];
							int node2 = coordInfoVecs.edges2Nodes_2[t];
							DIST += sqrt((coordInfoVecs.nodeLocX[node2] - coordInfoVecs.nodeLocX[node1])*(coordInfoVecs.nodeLocX[node2] - coordInfoVecs.nodeLocX[node1]) +
							(coordInfoVecs.nodeLocY[node2] - coordInfoVecs.nodeLocY[node1])*(coordInfoVecs.nodeLocY[node2] - coordInfoVecs.nodeLocY[node1]) + 
							(coordInfoVecs.nodeLocZ[node2] - coordInfoVecs.nodeLocZ[node1])*(coordInfoVecs.nodeLocZ[node2] - coordInfoVecs.nodeLocZ[node1]));
						}
					}
					// for (int t = 0; t < coordInfoVecs.num_edges; t++){
					// 	if (generalParams.boundaries_in_upperhem[t] == 1){
					// 		COUNT += 1.0;
					// 	}
					// }
					generalParams.length_scale = (DIST/COUNT)/generalParams.Rmin;
					std::cout<<"equilibrium length of each segment of the septin ring = "<<generalParams.length_scale*generalParams.Rmin<<std::endl;
					generalParams.eq_total_boundary_length = COUNT*generalParams.length_scale* generalParams.Rmin;
					std::cout<<"equilibrium length of the septin ring = "<<generalParams.eq_total_boundary_length<<std::endl;
					LINE_TENSION_START += 1;
			// 	}
				
			}
			//std::cout<<"start relaxation step"<<std::endl;
			EDGESWAP_ALGORITHM_TRIGGERED = false;
			bool end_of_relaxation = false;
			while (current_time < Max_Runtime){
				number_of_simulation_step += 1;
				if (Max_Runtime <= 0.0){
					std::cout<<"Max_Runtime is set to be 0 or negative! "<<std::endl;
					break;
				}
					
				Solve_Forces();
					
				if (LINE_TENSION_START >= 1){
					ComputeLineTensionSprings(
						generalParams,
						coordInfoVecs,
						linearSpringInfoVecs);
					}
				//std::cout<<"STOPPED BEFORE MemRepul"<<std::endl;
				/*energy_rep =
				ComputeMemRepulsionEnergy(
					coordInfoVecs,
					linearSpringInfoVecs, 
					capsidInfoVecs,
					generalParams,
						auxVecs);*/						

				AdvancePositions(
					coordInfoVecs,
					generalParams,
						domainParams);

				new_total_energy = linearSpringInfoVecs.linear_spring_energy + 
			areaTriangleInfoVecs.area_triangle_energy + 
			bendingTriangleInfoVecs.bending_triangle_energy;// +
			0.5*energy_rep;

		energy_gradient = sqrt((new_total_energy - old_total_energy)*(new_total_energy - old_total_energy))/old_total_energy;
		old_total_energy = new_total_energy;
		current_time+=generalParams.dt;	
	
				if (translate_counter % translate_frequency == 0){
					//	std::cout<<"SIMULATIONs TRIGGER REPOSITIONING AND EDGESWAP?"<<std::endl;

					newcenterX = 0.0;
					newcenterY = 0.0;
					newcenterZ = 0.0;
				//	std::cout<<"HERE?"<<std::endl;
					
					for (int i = 0; i < generalParams.maxNodeCount; i++){//for (int i = 0; i < coordInfoVecs.nodeLocX.size(); i++){
						//std::cout<<i<<std::endl;
						newcenterX += coordInfoVecs.nodeLocX[i];
						//std::cout<<newcenterX<<std::endl;
						newcenterY += coordInfoVecs.nodeLocY[i];
						//std::cout<<newcenterY<<std::endl;
						newcenterZ += coordInfoVecs.nodeLocZ[i];
						//std::cout<<newcenterZ<<std::endl;
					}
				//	std::cout<<"HERE2?"<<std::endl;
					newcenterX = newcenterX/generalParams.maxNodeCount; //coordInfoVecs.nodeLocX.size();
					newcenterY = newcenterY/generalParams.maxNodeCount; //coordInfoVecs.nodeLocX.size();
					newcenterZ = newcenterZ/generalParams.maxNodeCount; //coordInfoVecs.nodeLocX.size();
					displacementX = newcenterX - generalParams.centerX;
					displacementY = newcenterY - generalParams.centerY;
					displacementZ = newcenterZ - generalParams.centerZ;
					
				//	std::cout<<"HERE3?"<<std::endl;
					for (int i = 0; i < generalParams.maxNodeCount; i++){
					coordInfoVecs.nodeLocX[i] += -displacementX;
					coordInfoVecs.nodeLocY[i] += -displacementY;
					coordInfoVecs.nodeLocZ[i] += -displacementZ;
					}
				//	std::cout<<"HERE4?"<<std::endl;
					for (int i = 0; i < ljInfoVecs.LJ_PosX_all.size(); i++){
						ljInfoVecs.LJ_PosX_all[i] += -displacementX;
						ljInfoVecs.LJ_PosY_all[i] += -displacementY;
						ljInfoVecs.LJ_PosZ_all[i] += -displacementZ;
					}

					ComputeVolume(
						generalParams,
						coordInfoVecs,
						linearSpringInfoVecs,
						ljInfoVecs);

				}

		}

		end_of_relaxation = true;

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////Recalculate "u" since the system is relaxed, potentially becoming different ////////////////////
		////////////////////from the initial configuration /////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		max_cell_triangle_height = -10000.0;
		min_cell_triangle_height = 10000.0;
		for (int i = 0; i < coordInfoVecs.num_triangles; i++){
			if (coordInfoVecs.triangles2Nodes_1[i] >= (INT_MAX-100) || coordInfoVecs.triangles2Nodes_2[i] >= (INT_MAX-100) || coordInfoVecs.triangles2Nodes_3[i] >= (INT_MAX-100)){
				continue;
			}
			if (coordInfoVecs.triangles2Nodes_1[i] < 0 || coordInfoVecs.triangles2Nodes_2[i] < 0 || coordInfoVecs.triangles2Nodes_3[i] < 0){
				continue;
			}
			v1 = coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_1[i]];
			v2 = coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_2[i]];
			v3 = coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_3[i]];
			v4 = (v1+v2+v3)/3.0;
			if (v4 >= max_cell_triangle_height){
				max_cell_triangle_height = v4;
			}
			if (v4 <= min_cell_triangle_height){
				min_cell_triangle_height = v4;
			}
		}
		cell_height = (max_cell_triangle_height - min_cell_triangle_height);
		
		// double max_u = 3.0;//1.1;
		// double min_u = 0.05;//0.9;
		// vector<double> u;
		// coordInfoVecs.u.resize(3*coordInfoVecs.num_triangles);
		for (int i = 0; i < coordInfoVecs.num_triangles; i++){
			if (coordInfoVecs.triangles2Nodes_1[i] >= (INT_MAX-100) || coordInfoVecs.triangles2Nodes_2[i] >= (INT_MAX-100) || coordInfoVecs.triangles2Nodes_3[i] >= (INT_MAX-100)){
				
				continue;
			}
			if (coordInfoVecs.triangles2Nodes_1[i] < 0 || coordInfoVecs.triangles2Nodes_2[i] < 0 || coordInfoVecs.triangles2Nodes_3[i] < 0){
				continue;
			}
			v1 = coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_1[i]];
			v2 = coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_2[i]];
			v3 = coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_3[i]];
			v4 = (v1+v2+v3)/3.0;
			coordInfoVecs.u[i] = min_u + (max_u - min_u)*pow(((v4-min_cell_triangle_height)/cell_height),u_scalingPower); //Try exponential power? C_0*exp(-x/lambda), varying lambda
			// coordInfoVecs.u[i] = max_u - (max_u - min_u)*exp((-(v4-min_cell_triangle_height)/cell_height)/u_scalingPower);
		}
		
		if (polarization_postPeak == false){
			u_scalingPower = u_scalingPower+ u_scalingPower_Progress;
			if (u_scalingPower >= max_u_scalingPower){
				polarization_postPeak = true;
				u_scalingPower = max_u_scalingPower;
			}
		}
		else{
			u_scalingPower = u_scalingPower - u_scalingPower_postPeak_Progress;
			if (u_scalingPower <= max_u_scalingPower_postPeak){
				u_scalingPower = max_u_scalingPower_postPeak;
			}
		}
		
		// coordInfoVecs.beta = 1.0/(min_u + ((/*cell_height - */4.0*generalParams.Rmin)/cell_height)*(max_u - min_u));
	
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////Run chemical diffusion//////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// generalParams.chemdiff_maxruntime = 2e3;
		int MAX_LOOP = 1;
		// if (edgeswap_iteration < 1){
			// utilities_ptr->LDG_Surface_Diffusion_Initialize(
			// 		coordInfoVecs,
			// 		generalParams,
			// 		build_ptr->hostSetInfoVecs,
			// 		auxVecs,
			// 		rbc,
			// 		n_rbc
			// );
		// }
		// if (needToRebuildDiffStructAfterEdgeSwap == true){
		// if (triggered == true){
			// utilities_ptr->LDG_Surface_Diffusion_Structure_Rebuild(coordInfoVecs,
			// 													generalParams,
			// 													build_ptr->hostSetInfoVecs,
			// 													auxVecs,
			// 													rbc,
			// 													n_rbc);
			// needToRebuildDiffStructAfterEdgeSwap = false;
		// }
		// if (triggered == true || edgeswap_iteration == 0){
			if (triggered == true){
				std::cout<<"Growth is triggered"<<std::endl;
			}
		if (edgeswap_iteration == 0 || edgeswap_iteration%GROWTH_FREQUENCY==0){
			std::cout<<"After rebuilding 'u', the new highest triangle is at the height of "<<max_cell_triangle_height<<std::endl;
			std::cout<<"new 1/beta calculated: "<< 1.0/coordInfoVecs.beta<<std::endl;
			utilities_ptr->LDG_Surface_Diffusion_Structure_Rebuild(coordInfoVecs,
																generalParams,
																build_ptr->hostSetInfoVecs,
																auxVecs,
																rbc,
																n_rbc);
			needToRebuildDiffStructAfterEdgeSwap = false;
			std::cout<<"Start solving PDE"<<std::endl;
			for (int i = 0; i < MAX_LOOP; i++){
				utilities_ptr->LDG_Surface_Diffusion_Solve(
									edgeswap_iteration,
									coordInfoVecs,
									generalParams,
									build_ptr->hostSetInfoVecs,
									auxVecs,
									rbc,
									n_rbc);
				// storage->storeVariables();
			}
			triggered = false;
		}
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			
		if (end_of_relaxation == true){
			std::random_device rand_dev;
			// std::mt19937 generator2(rand_dev());
			std::mt19937 generator_edgeswap(rand_dev());
				ComputeVolume(
					generalParams,
					coordInfoVecs,
					linearSpringInfoVecs,
					ljInfoVecs);

				if ((generalParams.true_current_total_volume/initial_volume) < 0.6 || generalParams.true_current_total_volume/initial_volume >= MAX_VOLUME_RATIO){
					generalParams.true_num_edges = 0;
					for (int i = 0; i < coordInfoVecs.num_edges; i++){
						if (coordInfoVecs.edges2Nodes_1[i] != INT_MAX && coordInfoVecs.edges2Nodes_2[i] != INT_MAX){
							generalParams.true_num_edges += 1;
						}
					}
					storage-> print_VTK_File();
					storage-> storeVariables();
					
					if (generalParams.true_current_total_volume/initial_volume < 0.6){
						std::cout<<"Cell over compression 60%"<<std::endl;
					}
					else if (generalParams.true_current_total_volume/initial_volume >= MAX_VOLUME_RATIO){
						std::cout<<"Target volume ratio exceeded. Current volume ratio = "<<generalParams.true_current_total_volume/initial_volume<<std::endl;
					}
					std::cout<<"Current number of edgeswap iteration performed at volume-related termination = "<<edgeswap_iteration<<std::endl;
					std::cout<<"Current number of simulation step at volume-related termination = "<<number_of_simulation_step<<std::endl;

					Max_Runtime = 0.0;
					runSim = false;
					initial_kT = -1;
					break;
				}
					double current_bud_area = 0.0;
					for (int k = 0; k < coordInfoVecs.num_triangles; k++){
						if (coordInfoVecs.triangles2Nodes_1[k] >= (INT_MAX - 1000.0) || coordInfoVecs.triangles2Nodes_1[k] <= (-INT_MAX + 1000.0) ||
							coordInfoVecs.triangles2Nodes_2[k] >= (INT_MAX - 1000.0) || coordInfoVecs.triangles2Nodes_2[k] <= (-INT_MAX + 1000.0) ||
							coordInfoVecs.triangles2Nodes_3[k] >= (INT_MAX - 1000.0) || coordInfoVecs.triangles2Nodes_3[k] <= (-INT_MAX + 1000.0)){
									continue;
								}
						else{
							if (generalParams.triangles_in_upperhem[k] == 1){
								double r1x = coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_1[k]];
								double r1y = coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_1[k]];
								double r1z = coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_1[k]];
								double r2x = coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_2[k]];
								double r2y = coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_2[k]];
								double r2z = coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_2[k]];
								double r3x = coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_3[k]];
								double r3y = coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_3[k]];
								double r3z = coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_3[k]];
								double norm_r1r2 = sqrt((r2x-r1x)*(r2x-r1x) + (r2y-r1y)*(r2y-r1y) + (r2z-r1z)*(r2z-r1z));
								double norm_r2r3 = sqrt((r3x-r2x)*(r3x-r2x) + (r3y-r2y)*(r3y-r2y) + (r3z-r2z)*(r3z-r2z));
								double norm_r3r1 = sqrt((r3x-r1x)*(r3x-r1x) + (r3y-r1y)*(r3y-r1y) + (r3z-r1z)*(r3z-r1z));
								double s = (norm_r1r2 + norm_r2r3 + norm_r3r1)/2.0;
								double area = sqrt(s*(s-norm_r1r2)*(s-norm_r2r3)*(s-norm_r3r1));
								current_bud_area += area;
							}
						}
					}
					// std::cout<<"Current bud surface area = "<<current_bud_area<<std::endl;
					if (current_bud_area/Initial_Bud_Area >= MAX_BUD_AREA_RATIO){
						std::cout<<"Target bud surface area ratio exceeded. Current bud surface area ratio = "<<current_bud_area/Initial_Bud_Area<<std::endl;
						std::cout<<"Current number of edgeswap iteration performed at area-related termination = "<<edgeswap_iteration<<std::endl;
						std::cout<<"Current number of simulation step at area-related termination = "<<number_of_simulation_step<<std::endl;
						generalParams.true_num_edges = 0;
						for (int i = 0; i < coordInfoVecs.num_edges; i++){
							if (coordInfoVecs.edges2Nodes_1[i] != INT_MAX && coordInfoVecs.edges2Nodes_2[i] != INT_MAX){
								generalParams.true_num_edges += 1;
							}
						}
						storage-> print_VTK_File();
						storage-> storeVariables();
						Max_Runtime = 0.0;
						runSim = false;
						initial_kT = -1;
						break;
					}

				// std::cout<<"entering edge swap algorithm"<<std::endl;							
				utilities_ptr->transferDtoH(generalParams, coordInfoVecs, build_ptr->hostSetInfoVecs);
				
				// double max_conc = -INT_MAX;
				// double min_conc = INT_MAX;
				// double temp_max_conc, temp_min_conc;
				// for (int i = 0; i < coordInfoVecs.num_triangles; i++){
				// 	if (coordInfoVecs.triangles2Nodes_1[i] >= (INT_MAX - 1000) || coordInfoVecs.triangles2Nodes_1[i] <= (-INT_MAX + 1000)){
				// 		continue;
				// 	}
				// 	else if (coordInfoVecs.triangles2Nodes_2[i] >= (INT_MAX - 1000) || coordInfoVecs.triangles2Nodes_2[i] <= (-INT_MAX + 1000)){
				// 		continue;
				// 	}
				// 	else if (coordInfoVecs.triangles2Nodes_3[i] >= (INT_MAX - 1000) || coordInfoVecs.triangles2Nodes_3[i] <= (-INT_MAX + 1000)){
				// 		continue;
				// 	}

				// 	temp_max_conc = coordInfoVecs.soln_per_triangle[i];
				// 	if (temp_max_conc > max_conc){
				// 		max_conc = temp_max_conc;
				// 	}
				// 	if (temp_min_conc < min_conc){
				// 		min_conc = temp_min_conc;
				// 	}
				// }

				VectorShuffleForEdgeswapLoop.clear();
				for (int i = 0; i < coordInfoVecs.num_edges; i++){
					if (generalParams.edges_in_upperhem_list[i] >= 0 && 
						generalParams.edges_in_upperhem_list[i] != INT_MAX &&
						//generalParams.edges_in_upperhem[i] < coordInfoVecs.num_edges &&
						//generalParams.edges_in_upperhem[i] != -INT_MAX &&
						generalParams.boundaries_in_upperhem[i] != 1){
						// VectorShuffleForEdgeswapLoop.push_back(generalParams.edges_in_upperhem_list[i]);
						if (coordInfoVecs.edges2Nodes_1[i] < 0 || coordInfoVecs.edges2Nodes_1[i] >= (INT_MAX-1000)){
							continue;
						}
						else if (coordInfoVecs.edges2Nodes_2[i] < 0 || coordInfoVecs.edges2Nodes_2[i] >= (INT_MAX-1000)){
							continue;
						}
						// double tri1 = coordInfoVecs.edges2Triangles_1[i];
						// double tri2 = coordInfoVecs.edges2Triangles_2[i];
						// double avg_conc = (coordInfoVecs.soln_per_triangle[tri1] + coordInfoVecs.soln_per_triangle[tri2])/2.0;
						// if (avg_conc > max_conc_scaler_for_material_insert*max_conc){
						// 
							// VectorShuffleForGrowthLoop.push_back(y);
							// VectorShuffleForGrowthLoop_COUNT += 1;
							// VectorShuffleForEdgeswapLoop.push_back(i);
						// }
						VectorShuffleForEdgeswapLoop.push_back(i);
					}
				}	
				// std::cout<<"edge swap candidate vector size = "<<VectorShuffleForEdgeswapLoop.size()<<std::endl;
				num_edge_loop = round(true_num_edges_in_upperhem*SAMPLE_SIZE);
				if (num_edge_loop <= min_num_edge_loop){
					num_edge_loop = min_num_edge_loop;
				}
				
				std::shuffle(std::begin(VectorShuffleForEdgeswapLoop), std::end(VectorShuffleForEdgeswapLoop), generator_edgeswap);
				for (int edge_loop = 0; edge_loop < num_edge_loop; edge_loop++) {
													
					std::uniform_int_distribution<int> distribution(1,VectorShuffleForEdgeswapLoop.size());
					
					int dice_roll = distribution(generator_edgeswap);
					
					int edge = VectorShuffleForEdgeswapLoop[dice_roll - 1];
					//int edge = dice_roll -1;
					while (generalParams.boundaries_in_upperhem[edge] == 1 || edge == INT_MAX || edge < 0){
						dice_roll = distribution(generator_edgeswap);
						
						int edge = VectorShuffleForEdgeswapLoop[dice_roll - 1];
						//edge =  generalParams.edges_in_upperhem_list[dice_roll - 1];
						//edge = dice_roll -1;
						}
					//int edge = generalParams.edges_in_upperhem_list[edge_loop];
					//int edge = VectorShuffleForEdgeswapLoop[edge_loop];
					// std::cout<<"edge = "<<edge<<std::endl;
					if (edge < 0 || edge == INT_MAX){
						continue;
					}

					int ALPHA = utilities_ptr->edge_swap_host_vecs(
						edge,
						generalParams,
						build_ptr->hostSetInfoVecs,
						linearSpringInfoVecs,
						bendingTriangleInfoVecs,
						areaTriangleInfoVecs);
					if (ALPHA != 0){
						needToRebuildDiffStructAfterEdgeSwap = true;
					}
				}
			
			
			for (int i = 0; i < coordInfoVecs.num_triangles; i++){
				utilities_ptr->triangles2Triangles_host_vecs(i, build_ptr->hostSetInfoVecs,coordInfoVecs,generalParams, auxVecs);
			}
		
			utilities_ptr->transferHtoD(generalParams, coordInfoVecs, build_ptr->hostSetInfoVecs);//Currently this is treated as a backup of coordInfoVecs
	
			EDGESWAP_ALGORITHM_TRIGGERED = true;
			edgeswap_iteration += 1;
			translate_counter += 1;
		}
		
	
		if (EDGESWAP_ALGORITHM_TRIGGERED == false){
			//std::cout<<"current_time = "<<current_time<<std::endl;
			std::cout<<"EDGE_SWAP IS TRIGGERED BECAUSE PREVIOUS RELAXATION STEPS SOMEHOW FAIL TO TRIGGER EDGESWAP NORMALLY. PLEASE INVESTIGATE."<<std::endl;
			runSim = false;
			initial_kT = -1;
			Max_Runtime = 0.0;
			break;
		}
		
		if (edgeswap_iteration % (GROWTH_FREQUENCY*GROWTH_FREQUENCY_SCALE) == 0){
				for (int v = 0; v < coordInfoVecs.num_edges; v++){
				double ev1 = coordInfoVecs.edges2Nodes_1[v];
				double ev2 = coordInfoVecs.edges2Nodes_2[v];
				if (ev1 == INT_MAX || ev2 == INT_MAX){
					continue;
				}
				double ed = sqrt((coordInfoVecs.nodeLocX[ev2] - coordInfoVecs.nodeLocX[ev1])*(coordInfoVecs.nodeLocX[ev2] - coordInfoVecs.nodeLocX[ev1]) +
							(coordInfoVecs.nodeLocY[ev2] - coordInfoVecs.nodeLocY[ev1])*(coordInfoVecs.nodeLocY[ev2] - coordInfoVecs.nodeLocY[ev1]) +
							(coordInfoVecs.nodeLocZ[ev2] - coordInfoVecs.nodeLocZ[ev1])*(coordInfoVecs.nodeLocZ[ev2] - coordInfoVecs.nodeLocZ[ev1]));
				if (ed >= 2.0){
					std::cout<<"Edge over extension, possibly some instability occuring. Aborting the simulation."<<std::endl;
					runSim = false;
					initial_kT = -1;
					break;
				}
			}
			// generalParams.angle_per_edge.clear();
			generalParams.true_num_edges = 0;
			for (int i = 0; i < coordInfoVecs.num_edges; i++){
				if (coordInfoVecs.edges2Nodes_1[i] != INT_MAX && coordInfoVecs.edges2Nodes_2[i] != INT_MAX){
					generalParams.true_num_edges += 1;
				}
				}
				storage->print_VTK_File();
				storage->storeVariables();
				
				double current_bud_area = 0.0;
				for (int k = 0; k < coordInfoVecs.num_triangles; k++){
				if (coordInfoVecs.triangles2Nodes_1[k] >= (INT_MAX - 1000.0) || coordInfoVecs.triangles2Nodes_1[k] <= (-INT_MAX + 1000.0) ||
					coordInfoVecs.triangles2Nodes_2[k] >= (INT_MAX - 1000.0) || coordInfoVecs.triangles2Nodes_2[k] <= (-INT_MAX + 1000.0) ||
					coordInfoVecs.triangles2Nodes_3[k] >= (INT_MAX - 1000.0) || coordInfoVecs.triangles2Nodes_3[k] <= (-INT_MAX + 1000.0)){
							continue;
						}
				else{
					if (generalParams.triangles_in_upperhem[k] == 1){
						double r1x = coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_1[k]];
						double r1y = coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_1[k]];
						double r1z = coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_1[k]];
						double r2x = coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_2[k]];
						double r2y = coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_2[k]];
						double r2z = coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_2[k]];
						double r3x = coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_3[k]];
						double r3y = coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_3[k]];
						double r3z = coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_3[k]];
						double norm_r1r2 = sqrt((r2x-r1x)*(r2x-r1x) + (r2y-r1y)*(r2y-r1y) + (r2z-r1z)*(r2z-r1z));
						double norm_r2r3 = sqrt((r3x-r2x)*(r3x-r2x) + (r3y-r2y)*(r3y-r2y) + (r3z-r2z)*(r3z-r2z));
						double norm_r3r1 = sqrt((r3x-r1x)*(r3x-r1x) + (r3y-r1y)*(r3y-r1y) + (r3z-r1z)*(r3z-r1z));
						double s = (norm_r1r2 + norm_r2r3 + norm_r3r1)/2.0;
						double area = sqrt(s*(s-norm_r1r2)*(s-norm_r2r3)*(s-norm_r3r1));
						current_bud_area += area;
					}
				}
				}
				std::cout<<"Current bud surface area = "<<current_bud_area<<std::endl;
				std::cout<<"Current number of edgeswap performed = "<<edgeswap_iteration<<std::endl;
			//  std::cout<<"current Hill equation constant = "<<generalParams.hilleqnconst<<std::endl;
				//storage->storeVariables();
				std::cout<<"current total energy = "<< new_total_energy<<std::endl;
			//  std::cout<<"LINEAR ENERGY = "<<linearSpringInfoVecs.linear_spring_energy<<std::endl;
			// std::cout<<"BEND ENERGY = "<<bendingTriangleInfoVecs.bending_triangle_energy<<std::endl;
			// std::cout<<"AREA ENERGY = "<<areaTriangleInfoVecs.area_triangle_energy<<std::endl;
			//std::cout<<"REPULSION ENERGY = "<<energy_rep<<std::endl;
			// std::cout<<"VOLUME ENERGY = "<<generalParams.volume_energy<<std::endl;
				std::cout<<"energy_gradient = "<<energy_gradient<<std::endl;
				std::cout<<"true current total volume = "<<generalParams.true_current_total_volume<<std::endl;
			std::cout<<"equilibrium total volume = "<<generalParams.eq_total_volume<<std::endl;
			std::cout<<"u_scalingPower updated to : "<<u_scalingPower<<std::endl;
		}
		if (edgeswap_iteration == NKBT-1 ){
			true_num_edges_in_upperhem = 0;
			for (int i = 0; i < coordInfoVecs.num_edges; i++){
				if (generalParams.edges_in_upperhem_list[i] != INT_MAX && generalParams.edges_in_upperhem_list[i] >= 0){
					true_num_edges_in_upperhem += 1;
					//break;
				}
			}
			storage->print_VTK_File();
			storage->storeVariables();
			std::cout<<"Allowed number of simulation steps reached. Simulation terminated."<<std::endl;
			runSim = false;
			initial_kT = -1;
			break;
			}

//std::cout<<"ERROR BEFORE GROWTH"<<std::endl;
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////// GROWTH OF THE CELL (MEMBRANE) ////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
		// GROWTH_COUNTER = 0;
		// if (1 < 0){
		// if (edgeswap_iteration % GROWTH_FREQUENCY == 0){
		if (edgeswap_iteration % GROWTH_FREQUENCY == 0 && TOTAL_GROWTH_ATTEMPT < NUMBER_OF_GROWTH_EVENT){
			std::cout<<"Entering growth algorithm"<<std::endl;
			GROWTH_COUNTER += 1;
			// VectorShuffleForGrowthLoop.clear();
			// int VectorShuffleForGrowthLoop_COUNT = 0;
			// max_height = -10000.0;
			// double current_center_x = 0.0;
			// double current_center_y = 0.0;

			// for (int k = 0; k < generalParams.maxNodeCount; k++){
			// 	if (generalParams.nodes_in_upperhem[k] == 1){
			// 		current_center_x += coordInfoVecs.nodeLocX[k];
			// 		current_center_y += coordInfoVecs.nodeLocX[k];
			// 	}
				
			// 	if (coordInfoVecs. nodeLocZ[k] >= max_height){
			// 		max_height = coordInfoVecs.nodeLocZ[k];
			// 		max_height_index = k;
			// 	}

			// }
			// current_center_x = current_center_x/generalParams.maxNodeCount;
			// current_center_y = current_center_y/generalParams.maxNodeCount;
			// // double bdry_to_tip = 0.0;
			// double bdry_to_tip_height = 0.0;
			// for (int y = 0; y < boundary_edge_list.size(); y++){
			// 	// double edge_mdpt_x = (coordInfoVecs.nodeLocX[coordInfoVecs.edges2Nodes_1[boundary_edge_list[y]]] +
			// 	// 						coordInfoVecs.nodeLocX[coordInfoVecs.edges2Nodes_2[boundary_edge_list[y]]])/2.0;
			// 	// double edge_mdpt_y = (coordInfoVecs.nodeLocY[coordInfoVecs.edges2Nodes_1[boundary_edge_list[y]]] +
			// 	// 						coordInfoVecs.nodeLocY[coordInfoVecs.edges2Nodes_2[boundary_edge_list[y]]])/2.0;
			// 	double edge_mdpt_z = (coordInfoVecs.nodeLocZ[coordInfoVecs.edges2Nodes_1[boundary_edge_list[y]]] +
			// 							coordInfoVecs.nodeLocZ[coordInfoVecs.edges2Nodes_2[boundary_edge_list[y]]])/2.0;
			// 	// bdry_to_tip += sqrt(pow(current_center_x - edge_mdpt_x,2.0)+pow(current_center_y - edge_mdpt_y,2.0)+pow(coordInfoVecs.nodeLocZ[max_height_index] - edge_mdpt_z,2.0));
			// 	bdry_to_tip_height += sqrt(pow(coordInfoVecs.nodeLocZ[max_height_index] - edge_mdpt_z,2.0));
			// }
			// // bdry_to_tip = bdry_to_tip/boundary_edge_list.size();
			// bdry_to_tip_height = bdry_to_tip_height/boundary_edge_list.size();
			double max_conc = -INT_MAX;
			double min_conc = INT_MAX;
			double temp_max_conc, temp_min_conc;
			int max_conc_idx;
			for (int i = 0; i < coordInfoVecs.num_triangles; i++){
				if (coordInfoVecs.triangles2Nodes_1[i] >= (INT_MAX - 1000) || coordInfoVecs.triangles2Nodes_1[i] <= (-INT_MAX + 1000)){
					continue;
				}
				else if (coordInfoVecs.triangles2Nodes_2[i] >= (INT_MAX - 1000) || coordInfoVecs.triangles2Nodes_2[i] <= (-INT_MAX + 1000)){
					continue;
				}
				else if (coordInfoVecs.triangles2Nodes_3[i] >= (INT_MAX - 1000) || coordInfoVecs.triangles2Nodes_3[i] <= (-INT_MAX + 1000)){
					continue;
				}

				temp_max_conc = coordInfoVecs.soln_per_triangle[i];
				if (temp_max_conc > max_conc){
					max_conc = temp_max_conc;
					max_conc_idx = i;
				}
				if (temp_min_conc < min_conc){
					min_conc = temp_min_conc;
				}
			}
			std::cout<<"max_conc = "<<max_conc<<std::endl;
			std::cout<<"min_conc = "<<min_conc<<std::endl;
			int max_conc_triNode1 = coordInfoVecs.triangles2Nodes_1[max_conc_idx];
			int max_conc_triNode2 = coordInfoVecs.triangles2Nodes_2[max_conc_idx];
			int max_conc_triNode3 = coordInfoVecs.triangles2Nodes_3[max_conc_idx];
			std::cout<<"location of max conc triangle: "<<coordInfoVecs.nodeLocX[max_conc_triNode1]<<", "<<coordInfoVecs.nodeLocY[max_conc_triNode1]<<", "<<coordInfoVecs.nodeLocZ[max_conc_triNode1]<<std::endl;
			std::cout<<"location of max conc triangle: "<<coordInfoVecs.nodeLocX[max_conc_triNode2]<<", "<<coordInfoVecs.nodeLocY[max_conc_triNode2]<<", "<<coordInfoVecs.nodeLocZ[max_conc_triNode2]<<std::endl;
			std::cout<<"location of max conc triangle: "<<coordInfoVecs.nodeLocX[max_conc_triNode3]<<", "<<coordInfoVecs.nodeLocY[max_conc_triNode3]<<", "<<coordInfoVecs.nodeLocZ[max_conc_triNode3]<<std::endl;

			VectorShuffleForGrowthLoop.clear();
			int VectorShuffleForGrowthLoop_COUNT = 0;

			for (int y = 0; y < coordInfoVecs.num_edges; y++){
				// std::cout<<y<<std::endl;
				if (generalParams.edges_in_upperhem_list[y] >= 0 &&
					generalParams.edges_in_upperhem_list[y] != INT_MAX &&
					generalParams.edges_in_upperhem_list[y] <= (INT_MAX-1000) &&
					generalParams.edges_in_upperhem_list[y] >= (-INT_MAX+1000) &&
					generalParams.boundaries_in_upperhem[y] != 1){
						// std::cout<<"IF condition satisfied"<<std::endl;
						// std::cout<<"generalParams.edges_in_upperhem_list = "<<generalParams.edges_in_upperhem_list[y]<<std::endl;
						if (coordInfoVecs.edges2Nodes_1[y] < 0 || coordInfoVecs.edges2Nodes_1[y] >= (INT_MAX-1000)){
							continue;
						}
						else if (coordInfoVecs.edges2Nodes_2[y] < 0 || coordInfoVecs.edges2Nodes_2[y] >= (INT_MAX-1000)){
							continue;
						}
						// double edge_mdpt_x = (coordInfoVecs.nodeLocX[coordInfoVecs.edges2Nodes_1[y]] +
						// 					coordInfoVecs.nodeLocX[coordInfoVecs.edges2Nodes_2[y]])/2.0;
						// // std::cout<<edge_mdpt_x<<std::endl;
						// double edge_mdpt_y = (coordInfoVecs.nodeLocY[coordInfoVecs.edges2Nodes_1[y]] +
						// 						coordInfoVecs.nodeLocY[coordInfoVecs.edges2Nodes_2[y]])/2.0;
						// // std::cout<<edge_mdpt_y<<std::endl;
						// double edge_mdpt_z = (coordInfoVecs.nodeLocZ[coordInfoVecs.edges2Nodes_1[y]] +
						// 						coordInfoVecs.nodeLocZ[coordInfoVecs.edges2Nodes_2[y]])/2.0;
						// // std::cout<<edge_mdpt_z<<std::endl;
						// // double current_edge_to_tip = sqrt(pow(current_center_x - edge_mdpt_x,2.0)+pow(current_center_y - edge_mdpt_y,2.0)+pow(coordInfoVecs.nodeLocZ[max_height_index] - edge_mdpt_z,2.0));
						// double current_edge_to_tip_height = sqrt(pow(coordInfoVecs.nodeLocZ[max_height_index] - edge_mdpt_z,2.0));
						// double current_edge_to_tip_dist = sqrt(pow(coordInfoVecs.nodeLocX[max_height_index] - edge_mdpt_x,2.0) + pow(coordInfoVecs.nodeLocY[max_height_index] - edge_mdpt_y,2.0) + pow(coordInfoVecs.nodeLocZ[max_height_index] - edge_mdpt_z,2.0));
						// std::cout<<"current_edge_to_tip = "<<current_edge_to_tip<<std::endl;
					// if ((current_edge_to_tip/bdry_to_tip) <= 0.8 && bdry_to_tip >= (dtb*1.5)){
					// if ((current_edge_to_tip_height/bdry_to_tip_height) <= 0.33 && bdry_to_tip_height >= (dtb*1.5)){
					// if ((current_edge_to_tip_height) <= dtb*1.0 && bdry_to_tip_height >= (dtb*1.5)){
					// if ((current_edge_to_tip_height) <= generalParams.Rmin*current_edge_to_tip_height_scale && bdry_to_tip_height >= (generalParams.Rmin*bdry_to_tip_height_scale)){
					double tri1 = coordInfoVecs.edges2Triangles_1[y];
					double tri2 = coordInfoVecs.edges2Triangles_2[y];
					double avg_conc = (coordInfoVecs.soln_per_triangle[tri1] + coordInfoVecs.soln_per_triangle[tri2])/2.0;
					if (avg_conc > max_conc_scaler_for_material_insert*max_conc){
					// if (avg_conc > (max_conc - 0.1*(max_conc - min_conc))){
						VectorShuffleForGrowthLoop.push_back(y);
						VectorShuffleForGrowthLoop_COUNT += 1;
						// if (VectorShuffleForEdgeswapLoop.size() == 1){
						// 	std::cout<<"avg_conc = "<<avg_conc<<std::endl;
						// }
					}
					// if ((current_edge_to_tip_dist) <= generalParams.Rmin*current_edge_to_tip_dist_scale && bdry_to_tip_height >= (generalParams.Rmin*bdry_to_tip_height_scale)){
					// 	VectorShuffleForGrowthLoop.push_back(y);
					// 	VectorShuffleForGrowthLoop_COUNT += 1;
					// }
					// else if(bdry_to_tip < (dtb*1.5)){
					// else if(bdry_to_tip_height < (generalParams.Rmin*bdry_to_tip_height_scale)){
					// 	VectorShuffleForGrowthLoop.push_back(y);
					// 	VectorShuffleForGrowthLoop_COUNT += 1;
					// }
				}
				/*if (generalParams.edges_in_upperhem_list[y] >= 0 &&
					generalParams.edges_in_upperhem_list[y] != INT_MAX &&
					generalParams.boundaries_in_upperhem[y] != 1 &&
					edges_in_growth[y] == 1){
					VectorShuffleForGrowthLoop.push_back(y);
				}*/
				
				
			}
			
			// for (int y = 0; y < coordInfoVecs.num_edges; y++){
			// 	if (generalParams.edges_in_upperhem_list[y] >= 0 &&
			// 		generalParams.edges_in_upperhem_list[y] != INT_MAX &&
			// 		generalParams.boundaries_in_upperhem[y] != 1){
			// 		VectorShuffleForGrowthLoop.push_back(y);
			// 		VectorShuffleForGrowthLoop_COUNT += 1;
			// 	}
			// 	/*if (generalParams.edges_in_upperhem_list[y] >= 0 &&
			// 		generalParams.edges_in_upperhem_list[y] != INT_MAX &&
			// 		generalParams.boundaries_in_upperhem[y] != 1 &&
			// 		edges_in_growth[y] == 1){
			// 		VectorShuffleForGrowthLoop.push_back(y);
			// 	}*/
				
				
			// }
			std::cout<<"VectorShuffleForGrowthLoop_COUNT = "<<VectorShuffleForGrowthLoop_COUNT<<std::endl;

			std::random_device rand_dev;
			std::mt19937 generator3(rand_dev());
			std::shuffle(std::begin(VectorShuffleForGrowthLoop), std::end(VectorShuffleForGrowthLoop), generator3);
			int MAX_GROWTH_TEST = VectorShuffleForGrowthLoop.size();
			// bool triggered = false;
			int true_DELTA = 0;
			int MAX_GROWTH_PER_GROWTH_EVENT = 1;

			generalParams.triangle_undergoing_growth.clear();

			std::cout<<"BEGIN GROWTH ALGORITHM"<<std::endl;
			utilities_ptr->transferDtoH(generalParams, coordInfoVecs, build_ptr->hostSetInfoVecs);
			int GROWTH_COUNT = 0;
			for (int p = 0; p < MAX_GROWTH_TEST; p++){
				if (coordInfoVecs.edges2Nodes_1[VectorShuffleForGrowthLoop[p]] < 0 || coordInfoVecs.edges2Nodes_1[VectorShuffleForGrowthLoop[p]] == INT_MAX){
					continue;
				}
				else if (coordInfoVecs.edges2Nodes_2[VectorShuffleForGrowthLoop[p]] < 0 || coordInfoVecs.edges2Nodes_2[VectorShuffleForGrowthLoop[p]] == INT_MAX){
					continue;
				}
				//std::cout<<"begin growth test"<<std::endl;
				int DELTA = utilities_ptr->growth_host_vecs(
					VectorShuffleForGrowthLoop[p],
					generalParams,
					build_ptr->hostSetInfoVecs,
					coordInfoVecs,
					linearSpringInfoVecs,
					bendingTriangleInfoVecs,
					areaTriangleInfoVecs);
				if (DELTA >= 0){
					double tri1 = coordInfoVecs.edges2Triangles_1[DELTA];
					double tri2 = coordInfoVecs.edges2Triangles_2[DELTA];
					double avg_conc = (coordInfoVecs.soln_per_triangle[tri1] + coordInfoVecs.soln_per_triangle[tri2])/2.0;
					std::cout<<" edge chosen for growth has avg_conc of : "<<avg_conc<<std::endl;
					GROWTH_COUNT += 1;
					TOTAL_GROWTH_COUNTER += 1;
				}
				else{
					generalParams.triangle_undergoing_growth.clear();
				}

				if (GROWTH_COUNT >= MAX_GROWTH_PER_GROWTH_EVENT){
					break;
				}
			}
			TOTAL_GROWTH_ATTEMPT += 1;
			utilities_ptr->transferHtoD(generalParams, coordInfoVecs, build_ptr->hostSetInfoVecs);
			std::cout<<"END GROWTH ALGORITHM"<<std::endl;
			
			if (generalParams.triangle_undergoing_growth.size() > 1){
				triggered = true;
				std::cout<<"Begin rebuilding data structure post growth"<<std::endl;
					// // This means that growth actually happened. Time to rebuild the data structure and update
					// //the chemical concentration data structure as well.
					utilities_ptr->LDG_Surface_Diffusion_Structure_Rebuild_postGrowth(coordInfoVecs,
																			generalParams,
																			build_ptr->hostSetInfoVecs,
																			auxVecs,
																			rbc,
																			n_rbc);
				std::cout<<"End rebuilding data structure post growth"<<std::endl;
			//	storage->print_VTK_File();
			//	storage->storeVariables();
			}
			
			std::cout<<"number of cell wall insertion = "<<GROWTH_COUNT<<std::endl;
			std::cout<<"Total growth event triggered = "<<TOTAL_GROWTH_COUNTER<<std::endl;
			std::cout<<"Total growth event attempt = "<<TOTAL_GROWTH_ATTEMPT<<std::endl;
				// if (triggered == true){	
				// 	true_num_edges_in_upperhem = 0;
				// 	for (int i = 0; i < coordInfoVecs.num_edges; i++){
				// 		if (generalParams.edges_in_upperhem_list[i] != INT_MAX && generalParams.edges_in_upperhem_list[i] >= 0){
				// 			true_num_edges_in_upperhem += 1;
				// 			//break;
				// 		}
				// 	}
				// 	//std::cout<<"WHERE iS THE PROBLEM 3"<<std::endl;
				// }
				std::cout<<"Exiting growth algorithm"<<std::endl;
			}	
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
	coordInfoVecs.scaling_per_edge.resize(mem_prealloc*coordInfoVecs.num_edges, 0.0);
	hostSetInfoVecs.scaling_per_edge.resize(coordInfoVecs.scaling_per_edge.size(), 0.0);

	coordInfoVecs.soln_per_triangle.resize(mem_prealloc*coordInfoVecs.num_triangles, INT_MAX);
	coordInfoVecs.b_per_triangle.resize(mem_prealloc*coordInfoVecs.num_triangles, INT_MAX);

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

	coordInfoVecs.triangles2Triangles_1.resize( mem_prealloc*coordInfoVecs.num_triangles, -INT_MAX );
	coordInfoVecs.triangles2Triangles_2.resize( mem_prealloc*coordInfoVecs.num_triangles, -INT_MAX );
	coordInfoVecs.triangles2Triangles_3.resize( mem_prealloc*coordInfoVecs.num_triangles, -INT_MAX );

	hostSetInfoVecs.triangles2Triangles_1.resize( mem_prealloc*coordInfoVecs.num_triangles, -INT_MAX );
	hostSetInfoVecs.triangles2Triangles_2.resize( mem_prealloc*coordInfoVecs.num_triangles, -INT_MAX );
	hostSetInfoVecs.triangles2Triangles_3.resize( mem_prealloc*coordInfoVecs.num_triangles, -INT_MAX );

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
	//coordInfoVecs.nndata10.resize( mem_prealloc*generalParams.maxNodeCount);
	//coordInfoVecs.nndata11.resize( mem_prealloc*generalParams.maxNodeCount);
	//coordInfoVecs.nndata12.resize( mem_prealloc*generalParams.maxNodeCount);

	coordInfoVecs.SurfaceNormalX.resize( mem_prealloc*generalParams.maxNodeCount);
	coordInfoVecs.SurfaceNormalY.resize( mem_prealloc*generalParams.maxNodeCount);
	coordInfoVecs.SurfaceNormalZ.resize( mem_prealloc*generalParams.maxNodeCount);

	generalParams.nodes_in_upperhem.resize(mem_prealloc*generalParams.maxNodeCount);
	generalParams.triangles_in_upperhem.resize(mem_prealloc*coordInfoVecs.num_triangles);
	generalParams.edges_in_upperhem.resize(mem_prealloc*coordInfoVecs.num_edges);
	generalParams.edges_in_upperhem_list.resize(mem_prealloc*coordInfoVecs.num_edges);
	generalParams.boundaries_in_upperhem.resize(mem_prealloc*coordInfoVecs.num_edges, -1);

	hostSetInfoVecs.nodes_in_upperhem.resize(generalParams.nodes_in_upperhem.size());
	hostSetInfoVecs.triangles_in_upperhem.resize(generalParams.triangles_in_upperhem.size());
	hostSetInfoVecs.edges_in_upperhem.resize(generalParams.edges_in_upperhem.size());
	hostSetInfoVecs.edges_in_upperhem_list.resize(mem_prealloc*coordInfoVecs.num_edges);
	hostSetInfoVecs.boundaries_in_upperhem.resize(mem_prealloc*coordInfoVecs.num_edges, -1);

	hostSetInfoVecs.nodes2Triangles_1.resize(mem_prealloc*generalParams.maxNodeCount,-INT_MAX);
	hostSetInfoVecs.nodes2Triangles_2.resize(mem_prealloc*generalParams.maxNodeCount,-INT_MAX);
	hostSetInfoVecs.nodes2Triangles_3.resize(mem_prealloc*generalParams.maxNodeCount,-INT_MAX);
	hostSetInfoVecs.nodes2Triangles_4.resize(mem_prealloc*generalParams.maxNodeCount,-INT_MAX);
	hostSetInfoVecs.nodes2Triangles_5.resize(mem_prealloc*generalParams.maxNodeCount,-INT_MAX);
	hostSetInfoVecs.nodes2Triangles_6.resize(mem_prealloc*generalParams.maxNodeCount,-INT_MAX);
	hostSetInfoVecs.nodes2Triangles_7.resize(mem_prealloc*generalParams.maxNodeCount,-INT_MAX);
	hostSetInfoVecs.nodes2Triangles_8.resize(mem_prealloc*generalParams.maxNodeCount,-INT_MAX);
	hostSetInfoVecs.nodes2Triangles_9.resize(mem_prealloc*generalParams.maxNodeCount,-INT_MAX);
	

	coordInfoVecs.nodes2Triangles_1.resize(mem_prealloc*generalParams.maxNodeCount,-INT_MAX);
	coordInfoVecs.nodes2Triangles_2.resize(mem_prealloc*generalParams.maxNodeCount,-INT_MAX);
	coordInfoVecs.nodes2Triangles_3.resize(mem_prealloc*generalParams.maxNodeCount,-INT_MAX);
	coordInfoVecs.nodes2Triangles_4.resize(mem_prealloc*generalParams.maxNodeCount,-INT_MAX);
	coordInfoVecs.nodes2Triangles_5.resize(mem_prealloc*generalParams.maxNodeCount,-INT_MAX);
	coordInfoVecs.nodes2Triangles_6.resize(mem_prealloc*generalParams.maxNodeCount,-INT_MAX);
	coordInfoVecs.nodes2Triangles_7.resize(mem_prealloc*generalParams.maxNodeCount,-INT_MAX);
	coordInfoVecs.nodes2Triangles_8.resize(mem_prealloc*generalParams.maxNodeCount,-INT_MAX);
	coordInfoVecs.nodes2Triangles_9.resize(mem_prealloc*generalParams.maxNodeCount,-INT_MAX);
	

	thrust::copy(coordInfoVecs.nodes2Triangles_1.begin(), coordInfoVecs.nodes2Triangles_1.end(), hostSetInfoVecs.nodes2Triangles_1.begin() );
	thrust::copy(coordInfoVecs.nodes2Triangles_2.begin(), coordInfoVecs.nodes2Triangles_2.end(), hostSetInfoVecs.nodes2Triangles_2.begin() );
	thrust::copy(coordInfoVecs.nodes2Triangles_3.begin(), coordInfoVecs.nodes2Triangles_3.end(), hostSetInfoVecs.nodes2Triangles_3.begin() );
	thrust::copy(coordInfoVecs.nodes2Triangles_4.begin(), coordInfoVecs.nodes2Triangles_4.end(), hostSetInfoVecs.nodes2Triangles_4.begin() );
	thrust::copy(coordInfoVecs.nodes2Triangles_5.begin(), coordInfoVecs.nodes2Triangles_5.end(), hostSetInfoVecs.nodes2Triangles_5.begin() );
	thrust::copy(coordInfoVecs.nodes2Triangles_6.begin(), coordInfoVecs.nodes2Triangles_6.end(), hostSetInfoVecs.nodes2Triangles_6.begin() );
	thrust::copy(coordInfoVecs.nodes2Triangles_7.begin(), coordInfoVecs.nodes2Triangles_7.end(), hostSetInfoVecs.nodes2Triangles_7.begin() );
	thrust::copy(coordInfoVecs.nodes2Triangles_8.begin(), coordInfoVecs.nodes2Triangles_8.end(), hostSetInfoVecs.nodes2Triangles_8.begin() );
	thrust::copy(coordInfoVecs.nodes2Triangles_9.begin(), coordInfoVecs.nodes2Triangles_9.end(), hostSetInfoVecs.nodes2Triangles_9.begin() );
	//thrust::copy(coordInfoVecs.nodes2Triangles_10.begin(), coordInfoVecs.nodes2Triangles_10.end(), hostInfoVecs.nodes2Triangles_10.begin() );
	//thrust::copy(coordInfoVecs.nodes2Triangles_11.begin(), coordInfoVecs.nodes2Triangles_11.end(), hostInfoVecs.nodes2Triangles_11.begin() );
	//thrust::copy(coordInfoVecs.nodes2Triangles_12.begin(), coordInfoVecs.nodes2Triangles_12.end(), hostInfoVecs.nodes2Triangles_12.begin() );

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
	//thrust::copy(hostSetInfoVecs.nndata10.begin(), hostSetInfoVecs.nndata10.end(), coordInfoVecs.nndata10.begin() );
	//thrust::copy(hostSetInfoVecs.nndata11.begin(), hostSetInfoVecs.nndata11.end(), coordInfoVecs.nndata11.begin() );
	//thrust::copy(hostSetInfoVecs.nndata12.begin(), hostSetInfoVecs.nndata12.end(), coordInfoVecs.nndata12.begin() );

	coordInfoVecs.u.resize(mem_prealloc*coordInfoVecs.num_triangles);


 
	//allocate memory for other data structures.   

	//area triangle info vec
	//number of area springs is the number of triangles
	std::cout<<"Mem"<<std::endl;
	areaTriangleInfoVecs.tempNodeIdUnreduced.resize(mem_prealloc*areaTriangleInfoVecs.factor * coordInfoVecs.num_triangles);
	areaTriangleInfoVecs.tempNodeForceXUnreduced.resize(mem_prealloc*areaTriangleInfoVecs.factor * coordInfoVecs.num_triangles);
	areaTriangleInfoVecs.tempNodeForceYUnreduced.resize(mem_prealloc*areaTriangleInfoVecs.factor * coordInfoVecs.num_triangles);
	areaTriangleInfoVecs.tempNodeForceZUnreduced.resize(mem_prealloc*areaTriangleInfoVecs.factor * coordInfoVecs.num_triangles);
	// std::cout<<"HERE 1?"<<std::endl;
	areaTriangleInfoVecs.tempNodeIdReduced.resize(mem_prealloc*areaTriangleInfoVecs.factor * coordInfoVecs.num_triangles);
	areaTriangleInfoVecs.tempNodeForceXReduced.resize(mem_prealloc*areaTriangleInfoVecs.factor * coordInfoVecs.num_triangles);
	areaTriangleInfoVecs.tempNodeForceYReduced.resize(mem_prealloc*areaTriangleInfoVecs.factor * coordInfoVecs.num_triangles);
	areaTriangleInfoVecs.tempNodeForceZReduced.resize(mem_prealloc*areaTriangleInfoVecs.factor * coordInfoVecs.num_triangles);
	// std::cout<<"HERE 2?"<<std::endl;
	//beinding triangle info vec
	//num bending springs is the number of times each edge is between two triangles. 
	bendingTriangleInfoVecs.numBendingSprings = coordInfoVecs.num_edges;//coordInfoVecs.edges2Triangles_1.size();
	// std::cout<<"HERE 2.5?"<<std::endl;
	bendingTriangleInfoVecs.tempNodeIdUnreduced.resize(mem_prealloc*bendingTriangleInfoVecs.factor * bendingTriangleInfoVecs.numBendingSprings);
	bendingTriangleInfoVecs.tempNodeForceXUnreduced.resize(mem_prealloc*bendingTriangleInfoVecs.factor * bendingTriangleInfoVecs.numBendingSprings);
	bendingTriangleInfoVecs.tempNodeForceYUnreduced.resize(mem_prealloc*bendingTriangleInfoVecs.factor * bendingTriangleInfoVecs.numBendingSprings);
	bendingTriangleInfoVecs.tempNodeForceZUnreduced.resize(mem_prealloc*bendingTriangleInfoVecs.factor * bendingTriangleInfoVecs.numBendingSprings);
	// std::cout<<"HERE 3?"<<std::endl;
	bendingTriangleInfoVecs.tempNodeIdReduced.resize(mem_prealloc*bendingTriangleInfoVecs.factor * bendingTriangleInfoVecs.numBendingSprings);
	bendingTriangleInfoVecs.tempNodeForceXReduced.resize(mem_prealloc*bendingTriangleInfoVecs.factor * bendingTriangleInfoVecs.numBendingSprings);
	bendingTriangleInfoVecs.tempNodeForceYReduced.resize(mem_prealloc*bendingTriangleInfoVecs.factor * bendingTriangleInfoVecs.numBendingSprings);
	bendingTriangleInfoVecs.tempNodeForceZReduced.resize(mem_prealloc*bendingTriangleInfoVecs.factor * bendingTriangleInfoVecs.numBendingSprings);
	// std::cout<<"HERE 4?"<<std::endl;
	//linear springs
	
	linearSpringInfoVecs.tempNodeIdUnreduced.resize(mem_prealloc*linearSpringInfoVecs.factor*coordInfoVecs.num_edges);
	linearSpringInfoVecs.tempNodeForceXUnreduced.resize(mem_prealloc*linearSpringInfoVecs.factor*coordInfoVecs.num_edges);
	linearSpringInfoVecs.tempNodeForceYUnreduced.resize(mem_prealloc*linearSpringInfoVecs.factor*coordInfoVecs.num_edges);
	linearSpringInfoVecs.tempNodeForceZUnreduced.resize(mem_prealloc*linearSpringInfoVecs.factor*coordInfoVecs.num_edges);
	// std::cout<<"HERE 5?"<<std::endl;
	linearSpringInfoVecs.tempNodeIdReduced.resize(mem_prealloc*linearSpringInfoVecs.factor * coordInfoVecs.num_edges);
	linearSpringInfoVecs.tempNodeForceXReduced.resize(mem_prealloc*linearSpringInfoVecs.factor * coordInfoVecs.num_edges);
	linearSpringInfoVecs.tempNodeForceYReduced.resize(mem_prealloc*linearSpringInfoVecs.factor * coordInfoVecs.num_edges);
	linearSpringInfoVecs.tempNodeForceZReduced.resize(mem_prealloc*linearSpringInfoVecs.factor * coordInfoVecs.num_edges);
	// std::cout<<"HERE 6?"<<std::endl;
	linearSpringInfoVecs.edge_initial_length.clear();
	//linearSpringInfoVecs.edge_initial_length.resize(mem_prealloc*coordInfoVecs.num_edges,1.0);
	
	//thrust::copy(hostSetInfoVecs.edge_initial_length.begin(), hostSetInfoVecs.edge_initial_length.end(), linearSpringInfoVecs.edge_initial_length.begin() );

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
	//hostSetInfoVecs.nndata10.resize( mem_prealloc*generalParams.maxNodeCount);
	//hostSetInfoVecs.nndata11.resize( mem_prealloc*generalParams.maxNodeCount);
	//hostSetInfoVecs.nndata12.resize( mem_prealloc*generalParams.maxNodeCount);

	//std::cout<<"initial lengths: "<< linearSpringInfoVecs.edge_initial_length.size()<<std::endl;

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


