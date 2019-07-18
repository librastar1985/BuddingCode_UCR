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
#include "LJEnergy.h"
#include "NodeAdvance.h"
#include "BucketScheme.h"
#include "Storage.h" 
#include "Edgeswap_test.h"
#include "SystemBuilder.h"
#include <vector>

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
	
	//std::cout<<"Error here! before area"<<std::endl;
	ComputeAreaTriangleSprings(
		generalParams,
		coordInfoVecs,
		areaTriangleInfoVecs);
	
	//std::cout<<"Error here! before bending"<<std::endl;
	ComputeCosTriangleSprings(
		generalParams,
		coordInfoVecs,  
		bendingTriangleInfoVecs); 
	
	//std::cout<<"Error here! before memrepul"<<std::endl;
	ComputeMemRepulsionSprings(
		coordInfoVecs,
		linearSpringInfoVecs, 
		capsidInfoVecs,
		generalParams,
		auxVecs);
	
	ComputeLJSprings(
		coordInfoVecs,
		ljInfoVecs,
		generalParams);
		
};


void System::solveSystem() {
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
		//std::cout<<edgeIndices[i]<<std::endl;
	}
	

	auto it = remove_if(edgeIndices.begin(), edgeIndices.end(),  [](const int i) {return i < 0; });
	edgeIndices.erase(it, edgeIndices.end());

	//std::vector<int> edge_to_ljparticle;
	//generalParams.edge_to_ljparticle.reserve(coordInfoVecs.num_edges);
	for (int i = 0; i < coordInfoVecs.num_edges; i++){
		generalParams.edge_to_ljparticle.push_back(-1);
	};
	
	int node1, node2;
	double R1, R2;
	int last_index;
	double Influence_Range = ljInfoVecs.Rmin - 0.5;
	std::cout<<"Influence_Range = "<<Influence_Range<<std::endl;
	int num_edge_loop;
	double LJ_PosX_backup, LJ_PosY_backup, LJ_PosZ_backup;
	
	double Max_Runtime = 1.0;
	double Max_RunStep = Max_Runtime/generalParams.dt;
	std::cout<<"Max runtime = "<<Max_Runtime<<std::endl;
	std::cout<<"Max runstep = "<<Max_RunStep<<std::endl;
	linearSpringInfoVecs.spring_constant_att1 = 50.0;
	linearSpringInfoVecs.spring_constant_att2 = 0.75;
	linearSpringInfoVecs.spring_constant_rep1 = 0.01;
	linearSpringInfoVecs.spring_constant_rep2 = 12.0;
	std::cout<<"spring_constnat_att1 = "<<linearSpringInfoVecs.spring_constant_att1<<std::endl;
	std::cout<<"spring_constnat_att2 = "<<linearSpringInfoVecs.spring_constant_att2<<std::endl;
	std::cout<<"spring_constnat_rep1 = "<<linearSpringInfoVecs.spring_constant_rep1<<std::endl;
	std::cout<<"spring_constnat_rep2 = "<<linearSpringInfoVecs.spring_constant_rep2<<std::endl;
	
	bendingTriangleInfoVecs.initial_angle = 0.0;
	bendingTriangleInfoVecs.initial_angle_raft = 0.0;
	bendingTriangleInfoVecs.initial_angle_coat = 0.0;//0.167448079;
	
	bendingTriangleInfoVecs.spring_constant_raft = 0.0;//bendingTriangleInfoVecs.spring_constant;
	bendingTriangleInfoVecs.spring_constant_coat = 0.0;//bendingTriangleInfoVecs.spring_constant;
	
	std::cout<<"coat angle = "<<bendingTriangleInfoVecs.initial_angle_coat<<std::endl;
	generalParams.Rmin = 1.0;
	generalParams.abs_Rmin = 1.0;//0.586955;
	std::cout<<"Absolute minimum edge size = "<<generalParams.abs_Rmin<<std::endl;
	//std::cout<<"Morse D = "<<linearSpringInfoVecs.spring_constant_rep1<<std::endl;
	//std::cout<<"Morse a = "<<linearSpringInfoVecs.spring_constant_rep2<<std::endl;
	std::cout<<"spontaneous angle = "<<bendingTriangleInfoVecs.initial_angle<<std::endl;

	double initial_kT;
	initial_kT = generalParams.kT;//This is for the acceptance of change after looping through every edge within proximity.
		double SAMPLE_SIZE = 0.2;
	std::cout<<"Sample size: "<<SAMPLE_SIZE<<std::endl;
	auto edgeswap_ptr = std::make_shared<Edgeswap>(coordInfoVecs);

	int RECORD_TIME = 200;//round(Max_RunStep/2);
	std::cout<<"Record frequency = "<<RECORD_TIME<<std::endl;
	

	bool runSim = true;
	int NKBT = 4000; //The max number of edge-swap attempt per kBT value
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
	
	while (runSim == true){
		//WHEN += 1;
		double current_time = 0.0;

		//generalParams.kT = 1.0;//reset kT before simulations starts.
		//Max_Runtime = 0.0;//2.5;
			while (current_time < Max_Runtime){
					Solve_Forces();
					double energy_rep =
					ComputeMemRepulsionEnergy(
						coordInfoVecs,
						linearSpringInfoVecs, 
						capsidInfoVecs,
						generalParams,
						auxVecs);
						
					//now forces are computed, move nodes.
					AdvancePositions(
						coordInfoVecs,
						generalParams,
						domainParams);
					
					AdvanceLJParticle(
						generalParams,
						coordInfoVecs,
						ljInfoVecs);
							
					new_total_energy = linearSpringInfoVecs.linear_spring_energy + 
						areaTriangleInfoVecs.area_triangle_energy + 
						bendingTriangleInfoVecs.bending_triangle_energy + 
						0.5*energy_rep + 
						ljInfoVecs.lj_energy;

				energy_gradient = sqrt((new_total_energy - old_total_energy)*(new_total_energy - old_total_energy));
				//std::cout<<"energy_gradient "<<energy_gradient<<std::endl;
				if (sqrt((energy_gradient/old_total_energy)*(energy_gradient/old_total_energy)) <= 0.0001){
					old_total_energy = new_total_energy;
					//std::cout<<"STOPPED DUE TO CONVERGENCE"<<std::endl;
					break;
				}
				else{
					old_total_energy = new_total_energy;
				}
				current_time+=generalParams.dt;
				
			}
		std::cout<<"current time (1st iter before edgeswap): "<< current_time << std::endl;
		//std::cout<<"energy_difference "<<energy_gradient<<std::endl;
		std::cout<<"current total energy (1st iter before edgeswap) = "<<new_total_energy<<std::endl;
		std::cout<<"current KBT = "<<generalParams.kT<<std::endl;
		if (isnan(new_total_energy)==1){
			std::cout<<"Nan or Inf position update !!!!"<<std::endl;
			runSim = false;
			break;
		}
	
		//edgeswap_ptr->transferDtoH(coordInfoVecs, build_ptr->hostSetInfoVecs);//Currently this is treated as a backup of coordInfoVecs
		//storage->print_VTK_File();
		//storage->storeVariables();


		int edgeswap_iteration = 0;
		double preswap_energy = new_total_energy;
		double postswap_energy;
		double Ediff = 0.0;
		//initial_kT = generalParams.kT;
	
		while (initial_kT > 0){
					////////////////////NOW RELAX THE ATTEMPTED EDGESWAP//////////////////////
					current_time = 0.0;
					while (current_time < Max_Runtime){
						Solve_Forces();
						double energy_rep =
						ComputeMemRepulsionEnergy(
							coordInfoVecs,
							linearSpringInfoVecs, 
							capsidInfoVecs,
							generalParams,
							auxVecs);
							
						//now forces are computed, move nodes.
						AdvancePositions(
							coordInfoVecs,
							generalParams,
							domainParams);
						AdvanceLJParticle(
							generalParams,
							coordInfoVecs,
							ljInfoVecs);
							
								
						new_total_energy = linearSpringInfoVecs.linear_spring_energy + 
							areaTriangleInfoVecs.area_triangle_energy + 
							bendingTriangleInfoVecs.bending_triangle_energy +
							0.5*energy_rep +
							ljInfoVecs.lj_energy;
						//std::cout<<"new_total_energy = "<<new_total_energy<<std::endl;
		
						energy_gradient = sqrt((new_total_energy - old_total_energy)*(new_total_energy - old_total_energy));

						//std::cout<<"energy_gradient "<<energy_gradient<<std::endl;
						if (isnan(new_total_energy)==1){
							//runSim = false;
							std::cout<<"Infinite energy!!!!!!!!!!!!!!!!!"<<std::endl;
							edgeswap_ptr->transferHtoD(coordInfoVecs,build_ptr->hostSetInfoVecs);
							ljInfoVecs.LJ_PosX = LJ_PosX_backup;
							ljInfoVecs.LJ_PosY = LJ_PosY_backup;
							ljInfoVecs.LJ_PosZ = LJ_PosZ_backup;
							break;
							
						}
					//	else{
							if (sqrt((energy_gradient/old_total_energy)*(energy_gradient/old_total_energy)) <= 0.0001){
								old_total_energy = new_total_energy;
								std::cout<<"STOPPED DUE TO CONVERGENCE"<<std::endl;
								//std::cout<<"new_total_energy = "<<new_total_energy<<std::endl;
								WHEN +=1;
								//if (WHEN % RECORD_TIME == 0){
									//storage->print_VTK_File();
									std::cout<<"energy_gradient snapshot = "<<energy_gradient<<std::endl;
									std::cout<<"new_total_energy = "<<new_total_energy<<std::endl;
									//std::cout<<"current KBT = "<<generalParams.kT<<std::endl;
									//std::cout<<"WHEN = "<<WHEN<<std::endl;
									//storage->storeVariables();
								//}
								break;
							}
							else {
								old_total_energy = new_total_energy;
								if (current_time == Max_Runtime){
									std::cout<<"Max runtime reached"<<std::endl;
									//storage->print_VTK_File();
									std::cout<<"new_total_energy = "<<new_total_energy<<std::endl;
									}
								else{
									//if (WHEN % RECORD_TIME == 0){
									//	storage->print_VTK_File();
									//}									
								}
								WHEN +=1;
								//if (WHEN % RECORD_TIME == 0){
								//	storage->print_VTK_File();
								//	std::cout<<"energy_gradient snapshot = "<<energy_gradient<<std::endl;
								//	std::cout<<"new_total_energy = "<<new_total_energy<<std::endl;
									//std::cout<<"current KBT = "<<generalParams.kT<<std::endl;
									//std::cout<<"WHEN = "<<WHEN<<std::endl;
									//storage->storeVariables();
								}
								//break;
					//}
						
					//}
					current_time+=generalParams.dt;
					//std::cout<<current_time<<std::endl;
					}
					
					edgeswap_ptr->transferDtoH(coordInfoVecs, build_ptr->hostSetInfoVecs);//Currently this is treated as a backup of coordInfoVecs
					LJ_PosX_backup = ljInfoVecs.LJ_PosX;
					LJ_PosY_backup = ljInfoVecs.LJ_PosY;
					LJ_PosZ_backup = ljInfoVecs.LJ_PosZ;
					
					
					/*generalParams.edge_to_ljparticle.clear();
					for (int w = 0; w < coordInfoVecs.num_edges; w++){
						node1 = coordInfoVecs.edges2Nodes_1[w];
						node2 = coordInfoVecs.edges2Nodes_2[w];
						R1 = sqrt((coordInfoVecs.nodeLocX[node1] - ljInfoVecs.LJ_PosX)*(coordInfoVecs.nodeLocX[node1] - ljInfoVecs.LJ_PosX) +
							(coordInfoVecs.nodeLocY[node1] - ljInfoVecs.LJ_PosY)*(coordInfoVecs.nodeLocY[node1] - ljInfoVecs.LJ_PosY));

						R2 = sqrt((coordInfoVecs.nodeLocX[node2] - ljInfoVecs.LJ_PosX)*(coordInfoVecs.nodeLocX[node2] - ljInfoVecs.LJ_PosX) +
							(coordInfoVecs.nodeLocY[node2] - ljInfoVecs.LJ_PosY)*(coordInfoVecs.nodeLocY[node2] - ljInfoVecs.LJ_PosY));

						if (R1 < Influence_Range || R2 <Influence_Range){
							generalParams.edge_to_ljparticle.push_back(w);//if within range, record the index of the edge.
						}
						else{
							generalParams.edge_to_ljparticle.push_back(-1);//otherwise insert an -1.
						}
					}*/

					
					
					auto it = remove_if(edgeIndices.begin(), edgeIndices.end(),  [](const int i) {return i < 0; });

 					edgeIndices.erase(it, edgeIndices.end());
					
					
					num_edge_loop = round(edgeIndices.size()*SAMPLE_SIZE);
					
					if (edgeswap_iteration % RECORD_TIME == 0){
						storage->print_VTK_File();
						//storage->storeVariables();
					}
					edgeswap_iteration += 1;
					//std::cout<<"edgeswap_iteration = "<<edgeswap_iteration<<std::endl;
					if (edgeswap_iteration == NKBT){
						generalParams.kT = generalParams.kT - 0.072;
						std::cout<<"Current kBT = "<<generalParams.kT<<std::endl;
						edgeswap_iteration = 0;
					}

					if (generalParams.kT < min_kT){
						initial_kT = -1.0;
						runSim = false;
						break;
					}
					
					for (int edge_loop = 0; edge_loop < num_edge_loop; edge_loop++) {
						
						std::random_device rand_dev;
						std::mt19937 generator(rand_dev());
						std::uniform_int_distribution<int> distribution(1,edgeIndices.size());
						int dice_roll = distribution(generator);
						int edge = edgeIndices[dice_roll - 1];
						//int edge = edge_loop;
						int ALPHA = edgeswap_ptr->edge_swap_device_vecs(
							edge,
							generalParams,
							coordInfoVecs,
							//build_ptr->hostSetInfoVecs,
							linearSpringInfoVecs,
							bendingTriangleInfoVecs,
							areaTriangleInfoVecs);
						//num_edge_loop = num_edge_loop + ALPHA;
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
	coordInfoVecs.isNodeFixed.resize(hostSetInfoVecs.nodeLocX.size());
	coordInfoVecs.prevNodeLocX.resize(hostSetInfoVecs.nodeLocX.size());
	coordInfoVecs.prevNodeLocY.resize(hostSetInfoVecs.nodeLocX.size());
	coordInfoVecs.prevNodeLocZ.resize(hostSetInfoVecs.nodeLocX.size());

	coordInfoVecs.prevNodeForceX.resize(hostSetInfoVecs.nodeLocX.size());
	coordInfoVecs.prevNodeForceY.resize(hostSetInfoVecs.nodeLocX.size());
	coordInfoVecs.prevNodeForceZ.resize(hostSetInfoVecs.nodeLocX.size());
	
	coordInfoVecs.nodeLocX.resize(hostSetInfoVecs.nodeLocX.size());
	coordInfoVecs.nodeLocY.resize(hostSetInfoVecs.nodeLocX.size());
	coordInfoVecs.nodeLocZ.resize(hostSetInfoVecs.nodeLocX.size());

	coordInfoVecs.nodeForceX.resize(hostSetInfoVecs.nodeLocX.size());
	coordInfoVecs.nodeForceY.resize(hostSetInfoVecs.nodeLocX.size());
	coordInfoVecs.nodeForceZ.resize(hostSetInfoVecs.nodeLocX.size());

	coordInfoVecs.triangles2Nodes_1.resize( coordInfoVecs.num_triangles );
	coordInfoVecs.triangles2Nodes_2.resize( coordInfoVecs.num_triangles );
	coordInfoVecs.triangles2Nodes_3.resize( coordInfoVecs.num_triangles );
	
	coordInfoVecs.triangles2Edges_1.resize( coordInfoVecs.num_triangles );
	coordInfoVecs.triangles2Edges_2.resize( coordInfoVecs.num_triangles );
	coordInfoVecs.triangles2Edges_3.resize( coordInfoVecs.num_triangles );

	coordInfoVecs.edges2Nodes_1.resize( coordInfoVecs.num_edges );
	coordInfoVecs.edges2Nodes_2.resize( coordInfoVecs.num_edges );
	
	coordInfoVecs.edges2Triangles_1.resize( coordInfoVecs.num_edges );
	coordInfoVecs.edges2Triangles_2.resize( coordInfoVecs.num_edges );

	coordInfoVecs.nndata1.resize( generalParams.maxNodeCount);
	coordInfoVecs.nndata2.resize( generalParams.maxNodeCount);
	coordInfoVecs.nndata3.resize( generalParams.maxNodeCount);
	coordInfoVecs.nndata4.resize( generalParams.maxNodeCount);
	coordInfoVecs.nndata5.resize( generalParams.maxNodeCount);
	coordInfoVecs.nndata6.resize( generalParams.maxNodeCount);
	coordInfoVecs.nndata7.resize( generalParams.maxNodeCount);
	coordInfoVecs.nndata8.resize( generalParams.maxNodeCount);
	coordInfoVecs.nndata9.resize( generalParams.maxNodeCount);
	coordInfoVecs.nndata10.resize( generalParams.maxNodeCount);
	coordInfoVecs.nndata11.resize( generalParams.maxNodeCount);
	coordInfoVecs.nndata12.resize( generalParams.maxNodeCount);



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

	for (int k = 0; k < coordInfoVecs.isNodeFixed.size(); k++){
		bool isFixedHost = hostSetInfoVecs.isNodeFixed[k];
		bool isFixedDevice = coordInfoVecs.isNodeFixed[k];
		if (isFixedDevice != isFixedHost){

			std::cout<<"pos "<< k << " dev val = " << coordInfoVecs.isNodeFixed[k]
				<< " host val = " <<  hostSetInfoVecs.isNodeFixed[k] <<std::endl;
		}
	}
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
	areaTriangleInfoVecs.tempNodeIdUnreduced.resize(areaTriangleInfoVecs.factor * coordInfoVecs.num_triangles);
	areaTriangleInfoVecs.tempNodeForceXUnreduced.resize(areaTriangleInfoVecs.factor * coordInfoVecs.num_triangles);
	areaTriangleInfoVecs.tempNodeForceYUnreduced.resize(areaTriangleInfoVecs.factor * coordInfoVecs.num_triangles);
	areaTriangleInfoVecs.tempNodeForceZUnreduced.resize(areaTriangleInfoVecs.factor * coordInfoVecs.num_triangles);
	
	areaTriangleInfoVecs.tempNodeIdReduced.resize(areaTriangleInfoVecs.factor * coordInfoVecs.num_triangles);
	areaTriangleInfoVecs.tempNodeForceXReduced.resize(areaTriangleInfoVecs.factor * coordInfoVecs.num_triangles);
	areaTriangleInfoVecs.tempNodeForceYReduced.resize(areaTriangleInfoVecs.factor * coordInfoVecs.num_triangles);
	areaTriangleInfoVecs.tempNodeForceZReduced.resize(areaTriangleInfoVecs.factor * coordInfoVecs.num_triangles);

	//beinding triangle info vec
	//num bending springs is the number of times each edge is between two triangles. 
	bendingTriangleInfoVecs.numBendingSprings = coordInfoVecs.edges2Triangles_1.size();

	bendingTriangleInfoVecs.tempNodeIdUnreduced.resize(bendingTriangleInfoVecs.factor * bendingTriangleInfoVecs.numBendingSprings);
	bendingTriangleInfoVecs.tempNodeForceXUnreduced.resize(bendingTriangleInfoVecs.factor * bendingTriangleInfoVecs.numBendingSprings);
	bendingTriangleInfoVecs.tempNodeForceYUnreduced.resize(bendingTriangleInfoVecs.factor * bendingTriangleInfoVecs.numBendingSprings);
	bendingTriangleInfoVecs.tempNodeForceZUnreduced.resize(bendingTriangleInfoVecs.factor * bendingTriangleInfoVecs.numBendingSprings);
	
	bendingTriangleInfoVecs.tempNodeIdReduced.resize(bendingTriangleInfoVecs.factor * bendingTriangleInfoVecs.numBendingSprings);
	bendingTriangleInfoVecs.tempNodeForceXReduced.resize(bendingTriangleInfoVecs.factor * bendingTriangleInfoVecs.numBendingSprings);
	bendingTriangleInfoVecs.tempNodeForceYReduced.resize(bendingTriangleInfoVecs.factor * bendingTriangleInfoVecs.numBendingSprings);
	bendingTriangleInfoVecs.tempNodeForceZReduced.resize(bendingTriangleInfoVecs.factor * bendingTriangleInfoVecs.numBendingSprings);

	//linear springs
	linearSpringInfoVecs.tempNodeIdUnreduced.resize(linearSpringInfoVecs.factor * coordInfoVecs.num_edges);
	linearSpringInfoVecs.tempNodeForceXUnreduced.resize(linearSpringInfoVecs.factor * coordInfoVecs.num_edges);
	linearSpringInfoVecs.tempNodeForceYUnreduced.resize(linearSpringInfoVecs.factor * coordInfoVecs.num_edges);
	linearSpringInfoVecs.tempNodeForceZUnreduced.resize(linearSpringInfoVecs.factor * coordInfoVecs.num_edges);
	
	linearSpringInfoVecs.tempNodeIdReduced.resize(linearSpringInfoVecs.factor * coordInfoVecs.num_edges);
	linearSpringInfoVecs.tempNodeForceXReduced.resize(linearSpringInfoVecs.factor * coordInfoVecs.num_edges);
	linearSpringInfoVecs.tempNodeForceYReduced.resize(linearSpringInfoVecs.factor * coordInfoVecs.num_edges);
	linearSpringInfoVecs.tempNodeForceZReduced.resize(linearSpringInfoVecs.factor * coordInfoVecs.num_edges);
	
	linearSpringInfoVecs.edge_initial_length.resize(coordInfoVecs.num_edges);
	
	thrust::copy(hostSetInfoVecs.edge_initial_length.begin(), hostSetInfoVecs.edge_initial_length.end(), linearSpringInfoVecs.edge_initial_length.begin() );
	std::cout<<"initial lengths: "<< linearSpringInfoVecs.edge_initial_length.size()<<std::endl;

	std::cout<<"System Ready"<<std::endl;

	//Generate LJ particle list. and set LJ particle midpoint.
	double maxX_lj = *(thrust::max_element(coordInfoVecs.nodeLocX.begin(),coordInfoVecs.nodeLocX.end()));
	double minX_lj = *(thrust::min_element(coordInfoVecs.nodeLocX.begin(),coordInfoVecs.nodeLocX.end()));
	double maxY_lj = *(thrust::max_element(coordInfoVecs.nodeLocY.begin(),coordInfoVecs.nodeLocY.end()));
	double minY_lj = *(thrust::min_element(coordInfoVecs.nodeLocY.begin(),coordInfoVecs.nodeLocY.end()));
	
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


