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
//#include "LJEnergy.h"
#include "NodeAdvance.h"
#include "BucketScheme.h"
#include "Storage.h" 
#include "Edgeswap_test.h"
#include "SystemBuilder.h"
#include <vector>
#include "VolumeComp.h"
#include "VolumeSprings.h"

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
	//std::cout<<"Error here! before linear"<<std::endl;
	ComputeLinearSprings( 
		generalParams, 
		coordInfoVecs,
		linearSpringInfoVecs, 
		ljInfoVecs);
	//	std::cout<<"LINEAR"<<std::endl;
	
	//if (coordInfoVecs.nodeLocX.size() > 162){
	//	std::cout<<"force ="<<coordInfoVecs.nodeForceX[generalParams.maxNodeCount-1]<<" "<<coordInfoVecs.nodeForceY[generalParams.maxNodeCount-1]<<" "<<coordInfoVecs.nodeForceZ[generalParams.maxNodeCount-1]<<std::endl;
	//}
	//std::cout<<"Error here! before area"<<std::endl;
	ComputeAreaTriangleSprings(
		generalParams,
		coordInfoVecs,
		areaTriangleInfoVecs);
	//	std::cout<<"AREA"<<std::endl;
	
	//std::cout<<"Error here! before bending"<<std::endl;
	ComputeCosTriangleSprings(
		generalParams,
		coordInfoVecs,  
		bendingTriangleInfoVecs); 
	
	//std::cout<<"Error here! before memrepul"<<std::endl;
	/*ComputeMemRepulsionSprings(
		coordInfoVecs,
		linearSpringInfoVecs, 
		capsidInfoVecs,
		generalParams,
		auxVecs);*/

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

		
};


void System::solveSystem() {

	double displacementX, displacementY, displacementZ;
	double newcenterX, newcenterY, newcenterZ;

	//generalParams.num_of_triangles = coordInfoVecs.triangles2Nodes_1.size();
	//generalParams.maxNodeCount = coordInfoVecs.nodeLocX.size();
	//generalParams.num_of_edges = coordInfoVecs.edges2Nodes_1.size();
	//generalParams.eq_total_volume = 200.0;

	////IDENTIFY CENTER OF THE SPHERE////////////////////////////////////////
	////this is necessary to choose where to generate new lj points//////////
	////which will be located within a distance away from the center/////////
	for (int i = 0; i < generalParams.maxNodeCount; i++){
		generalParams.centerX += coordInfoVecs.nodeLocX[i];
		generalParams.centerY += coordInfoVecs.nodeLocY[i];
		generalParams.centerZ += coordInfoVecs.nodeLocZ[i];
	}
	generalParams.centerX = generalParams.centerX/generalParams.maxNodeCount;
	generalParams.centerY = generalParams.centerY/generalParams.maxNodeCount;
	generalParams.centerZ = generalParams.centerZ/generalParams.maxNodeCount;
	////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////

	std::vector<double> V1 = {0.0, 0.0};/*{ 0.2294,   -0.7103,    1.21788 ,  -1.3525 ,  -0.3047,
		0.8073,   -0.9284 ,  -0.9918 ,  -1.7507 ,  -2.1332,
		1.1312 ,  -0.2897 ,  -1.0996 ,  -0.3552 ,   0.9047};*/ //{-0.02283, -0.02283, -0.02283};

	std::vector<double> V2 = {0.0, 0.0};/*{-0.7403 ,   0.0554 ,  -0.21368  ,  1.0354 ,   2.0994,
		0.1924 ,  -1.7749 ,   0.4700 ,   0.5478 ,  -1.3549,
		0.4430 ,   0.6601 ,  -0.6469 ,  -1.0153 ,   0.9085};*/ //{2.384208, 2.384208, 1.68};

	std::vector<double> V3 = {-0.5, 0.5};/*{-0.9997 ,  -0.8749 ,  -1.82367 ,  -0.7812,    0.4490,
		2.0808 ,  -0.2730 ,   0.4642 ,   1.1056  , 0.6132,
		0.3975 ,  -0.9000 ,   2.1327  , -0.8614 ,  -0.6783};*/ //{1.5243396, 1.1043396, 1.5243396};
	
	
	for (int i = 0; i < V1.size(); i++){
		ljInfoVecs.LJ_PosX_all.push_back(V1[i]); 
		ljInfoVecs.LJ_PosY_all.push_back(V2[i]);
		ljInfoVecs.LJ_PosZ_all.push_back(V3[i]);
	}  
	//ljInfoVecs.LJ_PosX_all.push_back();
	//ljInfoVecs.LJ_PosY_all.push_back();
	//ljInfoVecs.LJ_PosZ_all.push_back();

	ljInfoVecs.forceX_all.reserve(ljInfoVecs.LJ_PosX_all.size());
	ljInfoVecs.forceY_all.reserve(ljInfoVecs.LJ_PosX_all.size());
	ljInfoVecs.forceZ_all.reserve(ljInfoVecs.LJ_PosX_all.size());

	generalParams.maxNodeCountLJ = ljInfoVecs.LJ_PosX_all.size();


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

	//std::vector<int> nodes_to_center;
	generalParams.nodes_in_upperhem.resize(generalParams.maxNodeCount);
	for (int i = 0; i < generalParams.maxNodeCount; i++){
		if (coordInfoVecs.nodeLocZ[i] > (generalParams.centerZ + 2.0)){
			generalParams.nodes_in_upperhem[i] = 1;
		}
		else{
			generalParams.nodes_in_upperhem[i] = -1;
		}
		
	}

	////////////////////////////////////////////////////////////////////////////////////////////



	////////////////////////////////////////////////////////////////////////////////////////////
	
	std::vector<int> edges_in_upperhem;
	generalParams.edges_in_upperhem.resize(coordInfoVecs.edges2Nodes_1.size());
	for (int i = 0; i < coordInfoVecs.edges2Nodes_1.size(); i++){
		int aaa = generalParams.nodes_in_upperhem[coordInfoVecs.edges2Nodes_1[i]];
		int bbb = generalParams.nodes_in_upperhem[coordInfoVecs.edges2Nodes_2[i]];
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
	generalParams.triangles_in_upperhem.resize(coordInfoVecs.triangles2Nodes_1.size());
	for (int i = 0; i < coordInfoVecs.triangles2Nodes_1.size(); i++){
		int aaa = generalParams.nodes_in_upperhem[coordInfoVecs.triangles2Nodes_1[i]];
		int bbb = generalParams.nodes_in_upperhem[coordInfoVecs.triangles2Nodes_2[i]];
		int ccc = generalParams.nodes_in_upperhem[coordInfoVecs.triangles2Nodes_3[i]];
		if ((aaa+bbb+ccc)==3){
			generalParams.triangles_in_upperhem[i] = 1;
			//triangles_in_upperhem.push_back(i);
		}
		if ((aaa+bbb+ccc)==2){
			generalParams.triangles_in_upperhem[i] = 0;
			//triangles_in_upperhem.push_back(i);
		}
		else{
			generalParams.triangles_in_upperhem[i] = -1;
		}
		
	}

	//std::vector<int> edge_to_ljparticle;
	//generalParams.edge_to_ljparticle.reserve(coordInfoVecs.num_edges);
	for (int i = 0; i < coordInfoVecs.num_edges; i++){
		generalParams.edge_to_ljparticle.push_back(-1);
	};
	
	int node1, node2;
	double R1, R2;
	int last_index;
	//double Influence_Range = ljInfoVecs.Rmin - 0.5;
	//std::cout<<"Influence_Range = "<<Influence_Range<<std::endl;
	int num_edge_loop;
	double LJ_PosX_backup, LJ_PosY_backup, LJ_PosZ_backup;
	
	double Max_Runtime = 0.05;
	double Max_RunStep = Max_Runtime/generalParams.dt;
	std::cout<<"Max runtime = "<<Max_Runtime<<std::endl;
	std::cout<<"Max runstep = "<<Max_RunStep<<std::endl;
	//linearSpringInfoVecs.spring_constant_att1 = 50.0;
	//linearSpringInfoVecs.spring_constant_att2 = 0.75;
	//linearSpringInfoVecs.spring_constant_rep1 = 0.01;
	//linearSpringInfoVecs.spring_constant_rep2 = 12.0;
	//std::cout<<"spring_constnat_att1 = "<<linearSpringInfoVecs.spring_constant_att1<<std::endl;
	//std::cout<<"spring_constnat_att2 = "<<linearSpringInfoVecs.spring_constant_att2<<std::endl;
	//std::cout<<"spring_constnat_rep1 = "<<linearSpringInfoVecs.spring_constant_rep1<<std::endl;
	//std::cout<<"spring_constnat_rep2 = "<<linearSpringInfoVecs.spring_constant_rep2<<std::endl;

	generalParams.volume_spring_constant = 10.0;//1.25;

	double scale = 20.0;
	std::cout<<"scaling of different region = "<<scale<<std::endl;
	linearSpringInfoVecs.spring_constant_weak = linearSpringInfoVecs.spring_constant/scale;
	bendingTriangleInfoVecs.spring_constant_weak = bendingTriangleInfoVecs.spring_constant/scale;
	areaTriangleInfoVecs.spring_constant_weak = areaTriangleInfoVecs.spring_constant/scale;
	
	//areaTriangleInfoVecs.initial_area = 0.6;

	bendingTriangleInfoVecs.initial_angle = 0.17549;//0.15;//0.0906;
	bendingTriangleInfoVecs.initial_angle_raft = 0.17549;//0.15;
	bendingTriangleInfoVecs.initial_angle_coat = 0.17549;//0.15;//0.167448079;
	
	bendingTriangleInfoVecs.spring_constant_raft = 0.0;//bendingTriangleInfoVecs.spring_constant;
	bendingTriangleInfoVecs.spring_constant_coat = 0.0;//bendingTriangleInfoVecs.spring_constant;
	bendingTriangleInfoVecs.spring_constant = bendingTriangleInfoVecs.spring_constant*(2/sqrt(3));
	bendingTriangleInfoVecs.spring_constant_raft = bendingTriangleInfoVecs.spring_constant_raft*(2/sqrt(3));
	bendingTriangleInfoVecs.spring_constant_coat = bendingTriangleInfoVecs.spring_constant_coat*(2/sqrt(3));
	std::cout<<"effective bending = "<<bendingTriangleInfoVecs.spring_constant<<std::endl;
	std::cout<<"effective bending raft = "<<bendingTriangleInfoVecs.spring_constant_raft<<std::endl;
	std::cout<<"effective bending coat = "<<bendingTriangleInfoVecs.spring_constant_coat<<std::endl;
	
	//std::cout<<"coat angle = "<<bendingTriangleInfoVecs.initial_angle_coat<<std::endl;
	generalParams.Rmin = 1.0;
	generalParams.abs_Rmin = 1.0;//0.586955;
	ljInfoVecs.Rmin_M = 0.75;
	ljInfoVecs.Rcutoff_M = 3.25;
	ljInfoVecs.Rmin_LJ = 1.0;//3.0//1.0;
	ljInfoVecs.Rcutoff_LJ = 1.1;//3.0;//1.0;
	//ljInfoVecs.epsilon_M = 1.0;
	ljInfoVecs.epsilon_M_rep1 = 0.5;//16.0;
	ljInfoVecs.epsilon_M_rep2 = 0.9;//1.0;
	//ljInfoVecs.epsilon_LJ = 0.25;
	ljInfoVecs.epsilon_LJ_rep1 = 0.125;//0.5;// 0.06;//7.5;
	ljInfoVecs.epsilon_LJ_rep2 = 1.0;//1.0;//1.0;//1.0;
	std::cout<<"Absolute minimum edge size = "<<generalParams.abs_Rmin<<std::endl;
	//std::cout<<"Morse D = "<<linearSpringInfoVecs.spring_constant_rep1<<std::endl;
	//std::cout<<"Morse a = "<<linearSpringInfoVecs.spring_constant_rep2<<std::endl;
	std::cout<<"spontaneous angle = "<<bendingTriangleInfoVecs.initial_angle<<std::endl;

	double initial_kT;
	initial_kT = generalParams.kT;//This is for the acceptance of change after looping through every edge within proximity.
		double SAMPLE_SIZE = 0.5;
	std::cout<<"Sample size: "<<SAMPLE_SIZE<<std::endl;
	auto edgeswap_ptr = std::make_shared<Edgeswap>(coordInfoVecs, generalParams);


	bool runSim = true;

	int GROWTH_TIME = 1;
	int RECORD_TIME = 20;//round(Max_RunStep/2);
	std::cout<<"Record frequency = "<<RECORD_TIME<<std::endl;
	std::cout<<"Growth frequency = "<<GROWTH_TIME<<std::endl;
	int NKBT = 400; //The max number of edge-swap attempt per kBT value
	std::cout<<"Number of edge-swap per kBT value = "<<NKBT<<std::endl;
	double min_kT = 0.21;
	std::cout<<"min kT for sim. termination = "<<min_kT<<std::endl;
	int WHEN = 0;
	double old_total_energy = 0.0;
	double new_total_energy = 0.0;
	double energy_gradient = 0.0;
	int Num_of_step_run = 0;
	auto build_ptr = weak_bld_ptr.lock();//upgrade weak builder to access host variables.
	std::cout<<"initial LJ-x : "<< ljInfoVecs.LJ_PosX <<std::endl;
	std::cout<<"initial LJ-y : "<< ljInfoVecs.LJ_PosY <<std::endl;
	std::cout<<"initial LJ-z : "<< ljInfoVecs.LJ_PosZ <<std::endl;
		

    
	double min_energy;
	//storage->print_VTK_File();
	////storage->storeVariables();

	////////////////////////////////////////
	
	ComputeVolume(
		generalParams,
		coordInfoVecs,
		linearSpringInfoVecs,
		ljInfoVecs
	);

	generalParams.eq_total_volume = generalParams.true_current_total_volume*1.50;//This is for setting different equilibrium volume to mimic growth or shirnkage.
	std::cout<<"true_current_total_volume = "<<generalParams.true_current_total_volume<<std::endl;
	std::cout<<"eq_total_volume = "<<generalParams.eq_total_volume<<std::endl;

	////////////////////////////////////////
	
	while (runSim == true){
		//WHEN += 1;
		double current_time = 0.0;

		//generalParams.kT = 1.0;//reset kT before simulations starts.
		//Max_Runtime = 0.0;//2.5;
		int translate_counter = 0;
			while (current_time < (2.0*Max_Runtime)){
					translate_counter += 1;
					Solve_Forces();
				
					/*double energy_rep =
					ComputeMemRepulsionEnergy(
						coordInfoVecs,
						linearSpringInfoVecs, 
						capsidInfoVecs,
						generalParams,
						auxVecs);*/

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
						
					}


					for (int i = 0; i < ljInfoVecs.LJ_PosX_all.size(); i++){
						ljInfoVecs.LJ_PosX_all[i] = ljInfoVecs.LJ_PosX_all[i] + generalParams.dt * ljInfoVecs.forceX_all[i];
						ljInfoVecs.LJ_PosY_all[i] = ljInfoVecs.LJ_PosY_all[i] + generalParams.dt * ljInfoVecs.forceY_all[i];
						ljInfoVecs.LJ_PosZ_all[i] = ljInfoVecs.LJ_PosZ_all[i] + generalParams.dt * ljInfoVecs.forceZ_all[i];
					
					}*/

				AdvancePositions(
					coordInfoVecs,
					generalParams,
					domainParams);
				if (translate_counter % 20 == 1){

					newcenterX = 0.0;
					newcenterY = 0.0;
					newcenterZ = 0.0;
					for (int i = 0; i < coordInfoVecs.nodeLocX.size(); i++){
						newcenterX += coordInfoVecs.nodeLocX[i];
						newcenterY += coordInfoVecs.nodeLocY[i];
						newcenterZ += coordInfoVecs.nodeLocZ[i];
					}
					newcenterX = newcenterX/coordInfoVecs.nodeLocX.size();
					newcenterY = newcenterY/coordInfoVecs.nodeLocX.size();
					newcenterZ = newcenterZ/coordInfoVecs.nodeLocX.size();
					displacementX = newcenterX - generalParams.centerX;
					displacementY = newcenterY - generalParams.centerY;
					displacementZ = newcenterZ - generalParams.centerZ;

					for (int i = 0; i < generalParams.maxNodeCount; i++){
						coordInfoVecs.nodeLocX[i] += -displacementX;
						coordInfoVecs.nodeLocY[i] += -displacementY;
						coordInfoVecs.nodeLocZ[i] += -displacementZ;
					}
				}
							
					new_total_energy = linearSpringInfoVecs.linear_spring_energy + 
						areaTriangleInfoVecs.area_triangle_energy + 
						bendingTriangleInfoVecs.bending_triangle_energy + 
						//0.5*energy_rep + 
						ljInfoVecs.lj_energy_M +
						ljInfoVecs.lj_energy_LJ;

				energy_gradient = sqrt((new_total_energy - old_total_energy)*(new_total_energy - old_total_energy));
				old_total_energy = new_total_energy;
				current_time+=generalParams.dt;
				

			}
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
		double preswap_energy = new_total_energy;
		double postswap_energy;
		double Ediff = 0.0;
		//initial_kT = generalParams.kT;
		num_edge_loop = 5;//round(edges_in_upperhem.size()*SAMPLE_SIZE);	
	
 		while (initial_kT > 0){
 					////////////////////NOW RELAX THE ATTEMPTED EDGESWAP//////////////////////
					 current_time = 0.0;
					 translate_counter = 0;
 					while (current_time < Max_Runtime){
						 translate_counter += 1;
						 Solve_Forces();

 						/*double energy_rep =
 						ComputeMemRepulsionEnergy(
 							coordInfoVecs,
 							linearSpringInfoVecs, 
 							capsidInfoVecs,
 							generalParams,
							 auxVecs);*/
					
				
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
							 
 						}
					
 						//now forces are computed, move nodes.
 						for (int i = 0; i < ljInfoVecs.LJ_PosX_all.size(); i++){
 							ljInfoVecs.LJ_PosX_all[i] = ljInfoVecs.LJ_PosX_all[i] + generalParams.dt * ljInfoVecs.forceX_all[i];
     						ljInfoVecs.LJ_PosY_all[i] = ljInfoVecs.LJ_PosY_all[i] + generalParams.dt * ljInfoVecs.forceY_all[i];
 							ljInfoVecs.LJ_PosZ_all[i] = ljInfoVecs.LJ_PosZ_all[i] + generalParams.dt * ljInfoVecs.forceZ_all[i];
						 }*/
						 
 						AdvancePositions(
 							coordInfoVecs,
 							generalParams,
							 domainParams);
						
						if (translate_counter % 20 == 1){
							newcenterX = 0.0;
							newcenterY = 0.0;
							newcenterZ = 0.0;
							for (int i = 0; i < coordInfoVecs.nodeLocX.size(); i++){
								newcenterX += coordInfoVecs.nodeLocX[i];
								newcenterY += coordInfoVecs.nodeLocY[i];
								newcenterZ += coordInfoVecs.nodeLocZ[i];
							}
							newcenterX = newcenterX/coordInfoVecs.nodeLocX.size();
							newcenterY = newcenterY/coordInfoVecs.nodeLocX.size();
							newcenterZ = newcenterZ/coordInfoVecs.nodeLocX.size();
							displacementX = newcenterX - generalParams.centerX;
							displacementY = newcenterY - generalParams.centerY;
							displacementZ = newcenterZ - generalParams.centerZ;
			
							for (int i = 0; i < generalParams.maxNodeCount; i++){
							coordInfoVecs.nodeLocX[i] += -displacementX;
							coordInfoVecs.nodeLocY[i] += -displacementY;
							coordInfoVecs.nodeLocZ[i] += -displacementZ;
							}

							//std::cout<<"ERROR 1"<<std::endl;
							edgeswap_ptr->transferDtoH(coordInfoVecs, build_ptr->hostSetInfoVecs);
							//std::cout<<"ERROR 1.5"<<std::endl;
							for (int edge_loop = 0; edge_loop < num_edge_loop; edge_loop++) {
				
								//std::random_device generator; //This is the older approach to generate random number, which may give bias 20190218
								std::random_device rand_dev;
								std::mt19937 generator(rand_dev());
							   //std::uniform_int_distribution<int> distribution(1,edgeIndices.size());
							   //std::uniform_int_distribution<int> distribution(1,edges_in_upperhem.size());
							   std::uniform_int_distribution<int> distribution(1,coordInfoVecs.num_edges);	
							   int dice_roll = distribution(generator);
							   //int edge = edgeIndices[dice_roll - 1];
							   int edge = dice_roll - 1;
							   //int edge = edges_in_upperhem[dice_roll - 1];
								//int ALPHA = edgeswap_ptr->edge_swap_device_vecs(
								int ALPHA = edgeswap_ptr->edge_swap_host_vecs(
									edge,
									generalParams,
									//coordInfoVecs,
									build_ptr->hostSetInfoVecs,
									linearSpringInfoVecs,
									bendingTriangleInfoVecs,
									areaTriangleInfoVecs);
								//num_edge_loop = num_edge_loop + ALPHA;
							}
							//std::cout<<"ERROR 2"<<std::endl;
							edgeswap_ptr->transferHtoD(coordInfoVecs, build_ptr->hostSetInfoVecs);//Currently this is treated as a backup of coordInfoVecs
							//std::cout<<"ERROR 2.5"<<std::endl;
							
						}
						
 						new_total_energy = linearSpringInfoVecs.linear_spring_energy + 
 							areaTriangleInfoVecs.area_triangle_energy + 
 							bendingTriangleInfoVecs.bending_triangle_energy +
 							//0.5*energy_rep +
 							ljInfoVecs.lj_energy_M +  
 							ljInfoVecs.lj_energy_LJ;
 						//std::cout<<"new_total_energy = "<<new_total_energy<<std::endl;

 						energy_gradient = sqrt((new_total_energy - old_total_energy)*(new_total_energy - old_total_energy));
 						
 					old_total_energy = new_total_energy;
 					current_time+=generalParams.dt;
					 }				
											
			
 								
 					if (edgeswap_iteration % RECORD_TIME == 0){
 						storage->print_VTK_File();
						 std::cout<<"current total energy = "<< new_total_energy<<std::endl;
						 std::cout<<"true current total volume = "<<generalParams.true_current_total_volume<<std::endl;
 					}
 					if (edgeswap_iteration % (1*RECORD_TIME) == 0){
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

					 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////// GROWTH OF THE CELL (MEMBRANE) ////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool triggered = false;
int triggered_counter = 0;
/*for (int p = 0; p < edges_in_upperhem_for_loop.size(); p++){//(int p = 0; p < generalParams.num_of_edges; p++){
	
		int k = generalParams.edges_in_upperhem_index[p];
		if (generalParams.edges_in_upperhem[k] == INT_MAX){
			continue;
		}
		//int k = p;
		//k -= triggered_counter;
		int iedge = k;
		//std::cout<<"node1 of iedge = "<<coordInfoVecs.edges2Nodes_1[iedge]<<std::endl;
		//std::cout<<"node2 of iedge = "<<coordInfoVecs.edges2Nodes_2[iedge]<<std::endl;
		int elem1 = coordInfoVecs.edges2Triangles_1[iedge];
		int elem2 = coordInfoVecs.edges2Triangles_2[iedge];
		//std::cout<<"elem1 of iedge = "<<elem1<<std::endl;
		//std::cout<<"elem2 of iedge = "<<elem2<<std::endl;
		
        int first_v = coordInfoVecs.triangles2Nodes_1[elem1];
        int second_v = coordInfoVecs.triangles2Nodes_2[elem1];
        int third_v = coordInfoVecs.triangles2Nodes_3[elem1];
        double v1x = coordInfoVecs.nodeLocX[second_v] - coordInfoVecs.nodeLocX[first_v];
        double v1y = coordInfoVecs.nodeLocY[second_v] - coordInfoVecs.nodeLocY[first_v];
        double v1z = coordInfoVecs.nodeLocZ[second_v] - coordInfoVecs.nodeLocZ[first_v];
        double v2x = coordInfoVecs.nodeLocX[third_v] - coordInfoVecs.nodeLocX[first_v];
        double v2y = coordInfoVecs.nodeLocY[third_v] - coordInfoVecs.nodeLocY[first_v];
        double v2z = coordInfoVecs.nodeLocZ[third_v] - coordInfoVecs.nodeLocZ[first_v];
        double This_area_v = sqrt((v1y*v2z - v2y*v1z)*(v1y*v2z - v2y*v1z) + 
                                ((-v1x*v2z) + (v2x*v1z))*((-v1x*v2z) + (v2x*v1z)) +
								(v1x*v2y - v2x*v1y)*(v1x*v2y - v2x*v1y))/2.0;
		int first_w = coordInfoVecs.triangles2Nodes_1[elem2];
        int second_w = coordInfoVecs.triangles2Nodes_2[elem2];
        int third_w = coordInfoVecs.triangles2Nodes_3[elem2];
        double w1x = coordInfoVecs.nodeLocX[second_w] - coordInfoVecs.nodeLocX[first_w];
        double w1y = coordInfoVecs.nodeLocY[second_w] - coordInfoVecs.nodeLocY[first_w];
        double w1z = coordInfoVecs.nodeLocZ[second_w] - coordInfoVecs.nodeLocZ[first_w];
        double w2x = coordInfoVecs.nodeLocX[third_w] - coordInfoVecs.nodeLocX[first_w];
        double w2y = coordInfoVecs.nodeLocY[third_w] - coordInfoVecs.nodeLocY[first_w];
        double w2z = coordInfoVecs.nodeLocZ[third_w] - coordInfoVecs.nodeLocZ[first_w];
        double This_area_w = sqrt((w1y*v2z - w2y*v1z)*(w1y*w2z - w2y*w1z) + 
                                ((-w1x*w2z) + (w2x*w1z))*((-w1x*w2z) + (w2x*w1z)) +
                                (w1x*w2y - w2x*w1y)*(w1x*w2y - w2x*w1y))/2.0;
        if ((This_area_v/ areaTriangleInfoVecs.initial_area >= 2.0) && (This_area_w/ areaTriangleInfoVecs.initial_area >= 2.0)){
            //triggered = true;
            triggered_counter += 1;
        }
        else{continue;}
		
		int t1e1, t1e2, t1e3, t2e1, t2e2, t2e3;

		if (coordInfoVecs.triangles2Edges_1[elem1] == iedge){
			t1e1 = coordInfoVecs.triangles2Edges_2[elem1];
			t1e2 = coordInfoVecs.triangles2Edges_3[elem1];
			t1e3 = coordInfoVecs.triangles2Edges_1[elem1];
		}
		else if (coordInfoVecs.triangles2Edges_2[elem1] == iedge){
			t1e1 = coordInfoVecs.triangles2Edges_3[elem1];
			t1e2 = coordInfoVecs.triangles2Edges_1[elem1];
			t1e3 = coordInfoVecs.triangles2Edges_2[elem1];
		} 
		else if (coordInfoVecs.triangles2Edges_3[elem1] == iedge){
			t1e1 = coordInfoVecs.triangles2Edges_1[elem1];
			t1e2 = coordInfoVecs.triangles2Edges_2[elem1];
			t1e3 = coordInfoVecs.triangles2Edges_3[elem1];
		}

		if (coordInfoVecs.triangles2Edges_1[elem2] == iedge){
			t2e1 = coordInfoVecs.triangles2Edges_2[elem2];
			t2e2 = coordInfoVecs.triangles2Edges_3[elem2];
			t2e3 = coordInfoVecs.triangles2Edges_1[elem2];
		}
		else if (coordInfoVecs.triangles2Edges_2[elem2] == iedge){
			t2e1 = coordInfoVecs.triangles2Edges_3[elem2];
			t2e2 = coordInfoVecs.triangles2Edges_1[elem2];
			t2e3 = coordInfoVecs.triangles2Edges_2[elem2];
		} 
		else if (coordInfoVecs.triangles2Edges_3[elem2] == iedge){
			t2e1 = coordInfoVecs.triangles2Edges_1[elem2];
			t2e2 = coordInfoVecs.triangles2Edges_2[elem2];
			t2e3 = coordInfoVecs.triangles2Edges_3[elem2];
		}

		
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

		if (coordInfoVecs.edges2Nodes_1[t2e1] == coordInfoVecs.edges2Nodes_1[iedge] || coordInfoVecs.edges2Nodes_1[t2e1] == coordInfoVecs.edges2Nodes_2[iedge]){
			n4 = coordInfoVecs.edges2Nodes_2[t2e1];
		}
		else if (coordInfoVecs.edges2Nodes_2[t2e1] == coordInfoVecs.edges2Nodes_1[iedge] || coordInfoVecs.edges2Nodes_2[t2e1] == coordInfoVecs.edges2Nodes_2[iedge]){
			n4 = coordInfoVecs.edges2Nodes_1[t2e1];
		}

		//std::cout<<"n1 = "<<n1<<std::endl;
		//std::cout<<"n2 = "<<n2<<std::endl;
		//std::cout<<"n3 = "<<n3<<std::endl;
		//std::cout<<"n4 = "<<n4<<std::endl;
		//These extract the indices of vertices of the selected triangles "elem1" and "elem2". Now we have n1, n2, n3, n4 in the correct orientation (supposedly).


		int edgeindex, a, a1, a2, a3, temp1, temp2;
		double newx = (coordInfoVecs.nodeLocX[coordInfoVecs.edges2Nodes_1[iedge]] + coordInfoVecs.nodeLocX[coordInfoVecs.edges2Nodes_2[iedge]])/2.0;
		coordInfoVecs.nodeLocX[generalParams. num_of_nodes] = newx;
		double newy = (coordInfoVecs.nodeLocY[coordInfoVecs.edges2Nodes_1[iedge]] + coordInfoVecs.nodeLocY[coordInfoVecs.edges2Nodes_2[iedge]])/2.0;
		coordInfoVecs.nodeLocY[generalParams. num_of_nodes] = newy;
		double newz = (coordInfoVecs.nodeLocZ[coordInfoVecs.edges2Nodes_1[iedge]] + coordInfoVecs.nodeLocZ[coordInfoVecs.edges2Nodes_2[iedge]])/2.0;
		coordInfoVecs.nodeLocZ[generalParams. num_of_nodes] = newz;
		//These are the coordinate of the new vertex. Its index is "coordInfoVecs.nodeLocX.size()-1"
		//int NODESIZE= generalParams.maxNodeCount;//coordInfoVecs.nodeLocX.size();
				
		coordInfoVecs.triangles2Nodes_1[generalParams. num_of_triangles] = n1;
		coordInfoVecs.triangles2Nodes_2[generalParams. num_of_triangles] = n2;
		coordInfoVecs.triangles2Nodes_3[generalParams. num_of_triangles] = generalParams.num_of_nodes;
		generalParams.num_of_triangles += 1;
		//NOTE: What this +1 actually does is that it specifies the location to write
		//any new data. Here it points to the location to write new triangles information.
		//This is a new triangle associated with (tn1, tn2, newnode). Its index is "coordInfoVecs.triangles2Nodes_1.size()-4".
		coordInfoVecs.triangles2Nodes_1[generalParams. num_of_triangles] =(n2);
		coordInfoVecs.triangles2Nodes_2[generalParams. num_of_triangles] =(n3);
		coordInfoVecs.triangles2Nodes_3[generalParams. num_of_triangles] = generalParams.num_of_nodes;
		generalParams.num_of_triangles += 1;
		//This is a new triangle associated with (tn2, tn3, newnode). Its index is "coordInfoVecs.triangles2Nodes_1.size()-3".
		coordInfoVecs.triangles2Nodes_1[generalParams. num_of_triangles] =(n3);
		coordInfoVecs.triangles2Nodes_2[generalParams. num_of_triangles] =(n4);
		coordInfoVecs.triangles2Nodes_3[generalParams. num_of_triangles] =generalParams.num_of_nodes;
		generalParams.num_of_triangles += 1;
		//This is a new triangle associated with (tn3, tn1, newnode). Its index is "coordInfoVecs.triangles2Nodes_1.size()-2".
		coordInfoVecs.triangles2Nodes_1[generalParams. num_of_triangles] =(n4);
		coordInfoVecs.triangles2Nodes_2[generalParams. num_of_triangles] =(n1);
		coordInfoVecs.triangles2Nodes_3[generalParams. num_of_triangles] =generalParams.num_of_nodes;
		generalParams.num_of_triangles += 1;
		//This is a new triangle associated with (tn3, tn1, newnode). Its index is "coordInfoVecs.triangles2Nodes_1.size()-1".
		
		//Now we add new edges formed by the addition of the new node.
		coordInfoVecs.edges2Nodes_1[generalParams.num_of_edges] = (generalParams.num_of_nodes);
		coordInfoVecs.edges2Nodes_2[generalParams.num_of_edges] = (n1);
		coordInfoVecs.edges2Triangles_1[generalParams.num_of_edges] = generalParams.num_of_triangles - 4;
		coordInfoVecs.edges2Triangles_2[generalParams.num_of_edges] = generalParams.num_of_triangles - 1;
		generalParams.num_of_edges += 1;
		//This is a new edge associated with (newnode, tn1). Its index is "edges2Nodes_1.size()-4".
		coordInfoVecs.edges2Nodes_1[generalParams.num_of_edges] = (generalParams.num_of_nodes);
		coordInfoVecs.edges2Nodes_2[generalParams.num_of_edges] = (n2);
		coordInfoVecs.edges2Triangles_1[generalParams.num_of_edges] = generalParams.num_of_triangles - 3;
		coordInfoVecs.edges2Triangles_2[generalParams.num_of_edges] = generalParams.num_of_triangles - 4;
		generalParams.num_of_edges += 1;
		//This is a new edge associated with (newnode, tn2). Its index is "edges2Nodes_1.size()-3".
		coordInfoVecs.edges2Nodes_1[generalParams.num_of_edges] = (generalParams.num_of_nodes);
		coordInfoVecs.edges2Nodes_2[generalParams.num_of_edges] = (n3);
		coordInfoVecs.edges2Triangles_1[generalParams.num_of_edges] = generalParams.num_of_triangles - 2;
		coordInfoVecs.edges2Triangles_2[generalParams.num_of_edges] = generalParams.num_of_triangles - 3;
		generalParams.num_of_edges += 1;
		//This is a new edge associated with (newnode, tn3). Its index is "edges2Nodes_1.size()-2".
		coordInfoVecs.edges2Nodes_1[generalParams.num_of_edges] = (generalParams.num_of_nodes);
		coordInfoVecs.edges2Nodes_2[generalParams.num_of_edges] = (n4);
		coordInfoVecs.edges2Triangles_1[generalParams.num_of_edges] = generalParams.num_of_triangles - 1;
		coordInfoVecs.edges2Triangles_2[generalParams.num_of_edges] = generalParams.num_of_triangles - 2;
		generalParams.num_of_edges += 1;
		for (int j = 0; j < 4; j++){			
			//Now we check to see if the order of update is correct, i.e. are edges2Triangles data in correct orientation.
			//This is crucial in the bendingspring computation.
			edgeindex = (generalParams.num_of_edges - (4-j));
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
		//generalParams.maxNodeCount += 1;
		generalParams.num_of_nodes += 1;

		


		
		//Now we update the edges2Triangles data structure with new edges.
		int TRIANGLESIZE = generalParams.num_of_triangles;//coordInfoVecs.triangles2Nodes_1.size();
		if (coordInfoVecs.edges2Triangles_1[t1e1] == elem1){
			coordInfoVecs.edges2Triangles_1[t1e1] = TRIANGLESIZE-4;
		}
		else if (coordInfoVecs.edges2Triangles_2[t1e1] == elem1){
			coordInfoVecs.edges2Triangles_2[t1e1] = TRIANGLESIZE-4;
		}
		else{}
		if (coordInfoVecs.edges2Triangles_1[t1e2] == elem1){
			coordInfoVecs.edges2Triangles_1[t1e2] = TRIANGLESIZE-3;
		}
		else if (coordInfoVecs.edges2Triangles_2[t1e2] == elem1){
			coordInfoVecs.edges2Triangles_2[t1e2] = TRIANGLESIZE-3;
		}
		else{}
		if (coordInfoVecs.edges2Triangles_1[t2e1] == elem2){
			coordInfoVecs.edges2Triangles_1[t2e1] = TRIANGLESIZE-2;
		}
		else if (coordInfoVecs.edges2Triangles_2[t2e1] == elem2){
			coordInfoVecs.edges2Triangles_2[t2e1] = TRIANGLESIZE-2;
		}
		if (coordInfoVecs.edges2Triangles_1[t2e2] == elem2){
			coordInfoVecs.edges2Triangles_1[t2e2] = TRIANGLESIZE-1;
		}
		else if (coordInfoVecs.edges2Triangles_2[t2e2] == elem2){
			coordInfoVecs.edges2Triangles_2[t2e2] = TRIANGLESIZE-1;
		}
		else{}
		//The above change the existing edges2Triangles data structure to accomodate new triangles added.

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		
		//Now we will take care of the last unedited data structure "triangles2Edges".
		//int aa, bb;
		int EDGESIZE = generalParams.num_of_edges;//coordInfoVecs.edges2Nodes_1.size();
		for (int j = 0; j < 4; j++){
			if (j == 0){
				coordInfoVecs.triangles2Edges_1[generalParams.num_of_triangles - 4] = (EDGESIZE-4);
				coordInfoVecs.triangles2Edges_2[generalParams.num_of_triangles - 4] = (t1e1);
				coordInfoVecs.triangles2Edges_3[generalParams.num_of_triangles - 4] = (EDGESIZE-3);   
			}
			else if (j == 1){
				coordInfoVecs.triangles2Edges_1[generalParams.num_of_triangles - 3] = (EDGESIZE-3);
				coordInfoVecs.triangles2Edges_2[generalParams.num_of_triangles - 3] = (t1e2);
				coordInfoVecs.triangles2Edges_3[generalParams.num_of_triangles - 3] = (EDGESIZE-2);   
			}
			else if (j ==2){
				coordInfoVecs.triangles2Edges_1[generalParams.num_of_triangles - 2] = (EDGESIZE-2);
				coordInfoVecs.triangles2Edges_2[generalParams.num_of_triangles - 2] = (t2e1);
				coordInfoVecs.triangles2Edges_3[generalParams.num_of_triangles - 2] = (EDGESIZE-1);   
			}
			else if (j ==3){
				coordInfoVecs.triangles2Edges_1[generalParams.num_of_triangles - 1] = (EDGESIZE-1);
				coordInfoVecs.triangles2Edges_2[generalParams.num_of_triangles - 1] = (t2e2);
				coordInfoVecs.triangles2Edges_3[generalParams.num_of_triangles - 1] = (EDGESIZE-4);   
			}
			
		}
	

		//Finally, we will fill the edge data chosen for growth (expansion) with INT_MAX so its data is no longer relevant to the computation
		
		coordInfoVecs.edges2Nodes_1[iedge] = INT_MAX;
		coordInfoVecs.edges2Nodes_2[iedge] = INT_MAX;
		for (int i = 0; i < generalParams.num_of_triangles; i++){
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
		

		
			coordInfoVecs.triangles2Nodes_1[elem1] = INT_MAX;
			coordInfoVecs.triangles2Nodes_2[elem1] = INT_MAX;
			coordInfoVecs.triangles2Nodes_3[elem1] = INT_MAX;
			coordInfoVecs.triangles2Nodes_1[elem2] = INT_MAX;
			coordInfoVecs.triangles2Nodes_2[elem2] = INT_MAX;
			coordInfoVecs.triangles2Nodes_3[elem2] = INT_MAX;
			//Delete the associated vertices information of the selected triangle.
			//Since we delete the chosen triangles, any triangle indexed lower than the deleted one will have its index reduced (or moved up) by 1.
			//Hence, we need to sweep through all data structures using the triangle index to change the index accordingly.
			for (int i = 0; i < generalParams.num_of_edges; i++){
				if (coordInfoVecs.edges2Triangles_1[i] == elem1 || coordInfoVecs.edges2Triangles_1[i] == elem2){
					coordInfoVecs.edges2Triangles_1[i] = INT_MAX;
				}
				if (coordInfoVecs.edges2Triangles_2[i] == elem1 || coordInfoVecs.edges2Triangles_2[i] == elem2){
					coordInfoVecs.edges2Triangles_2[i] = INT_MAX;
				}
			}
			//This completes the sweep. After this, the indices of triangle used in edges2Triangles data structure should be the correct one.
		
			coordInfoVecs.triangles2Edges_1[elem1] = INT_MAX;
			coordInfoVecs.triangles2Edges_2[elem1] = INT_MAX;
			coordInfoVecs.triangles2Edges_3[elem1] = INT_MAX;
			coordInfoVecs.triangles2Edges_1[elem2] = INT_MAX;
			coordInfoVecs.triangles2Edges_2[elem2] = INT_MAX;
			coordInfoVecs.triangles2Edges_3[elem2] = INT_MAX;
			for (int i = 0; i < generalParams.num_of_triangles; i++){
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
		if (generalParams.nodes_in_upperhem[coordInfoVecs.edges2Nodes_1[iedge]] == 1 && generalParams.nodes_in_upperhem[coordInfoVecs.edges2Nodes_2[iedge]] == 1){
			generalParams.nodes_in_upperhem.push_back(1);
		}
		else{
			generalParams.nodes_in_upperhem.push_back(-1);
		}

		
						//Erase the edge infomation related to the deleted triangle.

						//Now we update the nodes_in_upperhem and edges_in_upperhem data structures.
						//This ensures that the newly created edges will have the correct associated spring constant.
//std::cout<<"ERROR HERE?"<<std::endl;
		generalParams.edges_in_upperhem[iedge] = INT_MAX;
		for (int i = 0; i < generalParams.edges_in_upperhem_index.size(); i++){
			if (generalParams.edges_in_upperhem_index[i] == iedge){
				generalParams.edges_in_upperhem_index[i] == INT_MAX;
				//break;
			}
		}

						
		
		generalParams.triangles_in_upperhem[elem1] = INT_MAX;
		generalParams.triangles_in_upperhem[elem2] = INT_MAX;
		
						for (int q = 0; q < 4; q++){
							int nodeP = coordInfoVecs.triangles2Nodes_1[generalParams.num_of_triangles - (4-q)]; 
							int nodeQ = coordInfoVecs.triangles2Nodes_2[generalParams.num_of_triangles - (4-q)];
							int nodeR = coordInfoVecs.triangles2Nodes_3[generalParams.num_of_triangles - (4-q)];

							if (generalParams.nodes_in_upperhem[nodeP]==1 && generalParams.nodes_in_upperhem[nodeQ] ==1 && generalParams.nodes_in_upperhem[nodeR] ==1){
								generalParams.triangles_in_upperhem.push_back(1);
							}
							else{
								generalParams.triangles_in_upperhem.push_back(INT_MAX);
							}
						}
						for (int q = 0; q < 4; q++){
							int elem_1 = coordInfoVecs.edges2Triangles_1[generalParams.num_of_edges-(4 - q)];
							//std::cout<<nodeP<<std::endl;
							//std::cout<<generalParams.nodes_in_upperhem[nodeP]<<std::endl;
							int elem_2 = coordInfoVecs.edges2Triangles_2[generalParams.num_of_edges-(4 - q)];
							//std::cout<<nodeQ<<std::endl;
							//std::cout<<generalParams.nodes_in_upperhem[nodeQ]<<std::endl;
							
							if (generalParams.triangles_in_upperhem[elem_1] == 1 && generalParams.triangles_in_upperhem[elem_2] == 1){
								
								generalParams.edges_in_upperhem.push_back(1);
								generalParams.edges_in_upperhem_index.push_back(generalParams.num_of_edges - (4 - q));
							}
							else{
								
								generalParams.edges_in_upperhem.push_back(INT_MAX);
							}
							
							
						}

						//Here we regenerate the edges that we will loop through in both edgeswap and growth, if any edge deletion actually happened.
						
						if (triggered_counter > 0){
							edges_in_upperhem_for_loop.clear();
							for (int i = 0; i < generalParams.edges_in_upperhem_index.size(); i++){
								if (generalParams.edges_in_upperhem_index[i] != INT_MAX){
									edges_in_upperhem_for_loop.push_back(generalParams.edges_in_upperhem_index[i]);
								}
							}
						}
						
						//This should completes the dreadful data structure update associated with cell (membrane) growth.
						//Have fun modifying it if you need more function!
						//if (triggered == true){
						//	break;
						//}
					}*/
// if (triggered_counter > 0){
// generalParams.maxNodeCount = coordInfoVecs.nodeLocX.size();
// //std::cout<<"maxnodecount = " <<generalParams.maxNodeCount<<std::endl;
// generalParams.num_of_edges = coordInfoVecs.edges2Nodes_1.size();
// //std::cout<<"num_of_edges = " <<generalParams.num_of_edges<<std::endl;
// generalParams.num_of_triangles = coordInfoVecs.triangles2Nodes_1.size();
// //std::cout<<"num_of_triangles = " <<generalParams.num_of_triangles<<std::endl;

// }

//std::cout<<"GROWTH DONE!"<<std::endl;
 ////storage->print_VTK_File();
////storage->storeVariables();

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////// END OF GROWTH SECTION //////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

 					
					
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
	coordInfoVecs.isNodeFixed.resize(2*hostSetInfoVecs.nodeLocX.size());
	coordInfoVecs.prevNodeLocX.resize(2*hostSetInfoVecs.nodeLocX.size());
	coordInfoVecs.prevNodeLocY.resize(2*hostSetInfoVecs.nodeLocX.size());
	coordInfoVecs.prevNodeLocZ.resize(2*hostSetInfoVecs.nodeLocX.size());

	coordInfoVecs.prevNodeForceX.resize(2*hostSetInfoVecs.nodeLocX.size());
	coordInfoVecs.prevNodeForceY.resize(2*hostSetInfoVecs.nodeLocX.size());
	coordInfoVecs.prevNodeForceZ.resize(2*hostSetInfoVecs.nodeLocX.size());
	
	coordInfoVecs.nodeLocX.resize(2*hostSetInfoVecs.nodeLocX.size());
	coordInfoVecs.nodeLocY.resize(2*hostSetInfoVecs.nodeLocX.size());
	coordInfoVecs.nodeLocZ.resize(2*hostSetInfoVecs.nodeLocX.size());

	coordInfoVecs.nodeForceX.resize(2*hostSetInfoVecs.nodeLocX.size());
	coordInfoVecs.nodeForceY.resize(2*hostSetInfoVecs.nodeLocX.size());
	coordInfoVecs.nodeForceZ.resize(2*hostSetInfoVecs.nodeLocX.size());

	coordInfoVecs.triangles2Nodes_1.resize( 2*coordInfoVecs.num_triangles );
	coordInfoVecs.triangles2Nodes_2.resize( 2*coordInfoVecs.num_triangles );
	coordInfoVecs.triangles2Nodes_3.resize( 2*coordInfoVecs.num_triangles );
	
	coordInfoVecs.triangles2Edges_1.resize( 2*coordInfoVecs.num_triangles );
	coordInfoVecs.triangles2Edges_2.resize( 2*coordInfoVecs.num_triangles );
	coordInfoVecs.triangles2Edges_3.resize( 2*coordInfoVecs.num_triangles );

	coordInfoVecs.edges2Nodes_1.resize( 2*coordInfoVecs.num_edges );
	coordInfoVecs.edges2Nodes_2.resize( 2*coordInfoVecs.num_edges );
	
	coordInfoVecs.edges2Triangles_1.resize( 2*coordInfoVecs.num_edges );
	coordInfoVecs.edges2Triangles_2.resize( 2*coordInfoVecs.num_edges );

	coordInfoVecs.nndata1.resize( 2*generalParams.maxNodeCount);
	coordInfoVecs.nndata2.resize( 2*generalParams.maxNodeCount);
	coordInfoVecs.nndata3.resize( 2*generalParams.maxNodeCount);
	coordInfoVecs.nndata4.resize( 2*generalParams.maxNodeCount);
	coordInfoVecs.nndata5.resize( 2*generalParams.maxNodeCount);
	coordInfoVecs.nndata6.resize( 2*generalParams.maxNodeCount);
	coordInfoVecs.nndata7.resize( 2*generalParams.maxNodeCount);
	coordInfoVecs.nndata8.resize( 2*generalParams.maxNodeCount);
	coordInfoVecs.nndata9.resize( 2*generalParams.maxNodeCount);
	coordInfoVecs.nndata10.resize( 2*generalParams.maxNodeCount);
	coordInfoVecs.nndata11.resize( 2*generalParams.maxNodeCount);
	coordInfoVecs.nndata12.resize( 2*generalParams.maxNodeCount);



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
	areaTriangleInfoVecs.tempNodeIdUnreduced.resize(2*areaTriangleInfoVecs.factor * coordInfoVecs.num_triangles);
	areaTriangleInfoVecs.tempNodeForceXUnreduced.resize(2*areaTriangleInfoVecs.factor * coordInfoVecs.num_triangles);
	areaTriangleInfoVecs.tempNodeForceYUnreduced.resize(2*areaTriangleInfoVecs.factor * coordInfoVecs.num_triangles);
	areaTriangleInfoVecs.tempNodeForceZUnreduced.resize(2*areaTriangleInfoVecs.factor * coordInfoVecs.num_triangles);
	
	areaTriangleInfoVecs.tempNodeIdReduced.resize(2*areaTriangleInfoVecs.factor * coordInfoVecs.num_triangles);
	areaTriangleInfoVecs.tempNodeForceXReduced.resize(2*areaTriangleInfoVecs.factor * coordInfoVecs.num_triangles);
	areaTriangleInfoVecs.tempNodeForceYReduced.resize(2*areaTriangleInfoVecs.factor * coordInfoVecs.num_triangles);
	areaTriangleInfoVecs.tempNodeForceZReduced.resize(2*areaTriangleInfoVecs.factor * coordInfoVecs.num_triangles);

	//beinding triangle info vec
	//num bending springs is the number of times each edge is between two triangles. 
	bendingTriangleInfoVecs.numBendingSprings = coordInfoVecs.num_edges;//coordInfoVecs.edges2Triangles_1.size();

	bendingTriangleInfoVecs.tempNodeIdUnreduced.resize(2*bendingTriangleInfoVecs.factor * bendingTriangleInfoVecs.numBendingSprings);
	bendingTriangleInfoVecs.tempNodeForceXUnreduced.resize(2*bendingTriangleInfoVecs.factor * bendingTriangleInfoVecs.numBendingSprings);
	bendingTriangleInfoVecs.tempNodeForceYUnreduced.resize(2*bendingTriangleInfoVecs.factor * bendingTriangleInfoVecs.numBendingSprings);
	bendingTriangleInfoVecs.tempNodeForceZUnreduced.resize(2*bendingTriangleInfoVecs.factor * bendingTriangleInfoVecs.numBendingSprings);
	
	bendingTriangleInfoVecs.tempNodeIdReduced.resize(2*bendingTriangleInfoVecs.factor * bendingTriangleInfoVecs.numBendingSprings);
	bendingTriangleInfoVecs.tempNodeForceXReduced.resize(2*bendingTriangleInfoVecs.factor * bendingTriangleInfoVecs.numBendingSprings);
	bendingTriangleInfoVecs.tempNodeForceYReduced.resize(2*bendingTriangleInfoVecs.factor * bendingTriangleInfoVecs.numBendingSprings);
	bendingTriangleInfoVecs.tempNodeForceZReduced.resize(2*bendingTriangleInfoVecs.factor * bendingTriangleInfoVecs.numBendingSprings);

	//linear springs
	
	linearSpringInfoVecs.tempNodeIdUnreduced.resize(2*linearSpringInfoVecs.factor*coordInfoVecs.num_edges);
	linearSpringInfoVecs.tempNodeForceXUnreduced.resize(2*linearSpringInfoVecs.factor*coordInfoVecs.num_edges);
	linearSpringInfoVecs.tempNodeForceYUnreduced.resize(2*linearSpringInfoVecs.factor*coordInfoVecs.num_edges);
	linearSpringInfoVecs.tempNodeForceZUnreduced.resize(2*linearSpringInfoVecs.factor*coordInfoVecs.num_edges);
	
	linearSpringInfoVecs.tempNodeIdReduced.resize(2*linearSpringInfoVecs.factor * coordInfoVecs.num_edges);
	linearSpringInfoVecs.tempNodeForceXReduced.resize(2*linearSpringInfoVecs.factor * coordInfoVecs.num_edges);
	linearSpringInfoVecs.tempNodeForceYReduced.resize(2*linearSpringInfoVecs.factor * coordInfoVecs.num_edges);
	linearSpringInfoVecs.tempNodeForceZReduced.resize(2*linearSpringInfoVecs.factor * coordInfoVecs.num_edges);
	
	linearSpringInfoVecs.edge_initial_length.resize(2*coordInfoVecs.num_edges,1.0);
	
	thrust::copy(hostSetInfoVecs.edge_initial_length.begin(), hostSetInfoVecs.edge_initial_length.end(), linearSpringInfoVecs.edge_initial_length.begin() );

	//Resize the hostSetInfoVecs so that we can copy data back and forth between hostSetinfoVecs and coordInfoVecs without problem.
	hostSetInfoVecs.isNodeFixed.resize(2*hostSetInfoVecs.nodeLocX.size());
	//hostSetInfoVecs.prevNodeLocX.resize(2*hostSetInfoVecs.nodeLocX.size());
	//hostSetInfoVecs.prevNodeLocY.resize(2*hostSetInfoVecs.nodeLocX.size());
	//hostSetInfoVecs.prevNodeLocZ.resize(2*hostSetInfoVecs.nodeLocX.size());

	//hostSetInfoVecs.prevNodeForceX.resize(2*hostSetInfoVecs.nodeLocX.size());
	//hostSetInfoVecs.prevNodeForceY.resize(2*hostSetInfoVecs.nodeLocX.size());
	//hostSetInfoVecs.prevNodeForceZ.resize(2*hostSetInfoVecs.nodeLocX.size());
	
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

	/*hostSetInfoVecs.nndata1.resize( 2*generalParams.maxNodeCount);
	hostSetInfoVecs.nndata2.resize( 2*generalParams.maxNodeCount);
	hostSetInfoVecs.nndata3.resize( 2*generalParams.maxNodeCount);
	hostSetInfoVecs.nndata4.resize( 2*generalParams.maxNodeCount);
	hostSetInfoVecs.nndata5.resize( 2*generalParams.maxNodeCount);
	hostSetInfoVecs.nndata6.resize( 2*generalParams.maxNodeCount);
	hostSetInfoVecs.nndata7.resize( 2*generalParams.maxNodeCount);
	hostSetInfoVecs.nndata8.resize( 2*generalParams.maxNodeCount);
	hostSetInfoVecs.nndata9.resize( 2*generalParams.maxNodeCount);
	hostSetInfoVecs.nndata10.resize( 2*generalParams.maxNodeCount);
	hostSetInfoVecs.nndata11.resize( 2*generalParams.maxNodeCount);
	hostSetInfoVecs.nndata12.resize( 2*generalParams.maxNodeCount);*/

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


