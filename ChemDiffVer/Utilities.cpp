#include "System.h"
#include <random>
#include "Utilities.h"
#include <math.h>
#include <cmath>
#include <stdio.h>
#include "cdecs.h"
#include "vmalloc.h"
#include "cleanup.c"
#include "vectormalloc.c"
#include "screen.c"
#include "ppsub.c"
#include "error.c"
#include "debug.c"
#include "output.c"
#include "simpleio.c"
#include "times.c"
#include "other.c"
#include "machine.c"
#include "fgetstrin.c"
//#include "General.h"
//#include <gintrp>
// #include "gas/gdecs/gdecs.h"
// #include "driver/ddecs.h"
#include <omp.h>


int       Gauss_N = 2;
double    q[9], qw[9];

typedef  void	    *POINTER;
#define Dot2d(A,B)							\
	((A)[0]*(B)[0] + (A)[1]*(B)[1])

#define Dot3d(A,B)							\
	((A)[0]*(B)[0] + (A)[1]*(B)[1] + (A)[2]*(B)[2])

#define	Mag3d(A) sqrt(Dot3d(A,A))

#define Cross3d(B,C,ans)						\
	{								\
		(ans)[0] = ((B)[1])*((C)[2]) - ((B)[2])*((C)[1]);	\
		(ans)[1] = ((B)[2])*((C)[0]) - ((B)[0])*((C)[2]);	\
		(ans)[2] = ((B)[0])*((C)[1]) - ((B)[1])*((C)[0]);	\
	}

#define QCross3d(B,C,ans)						\
		ans ## 0 = B ## 1*C ## 2 - B ## 2*C ## 1;		\
		ans ## 1 = B ## 2*C ## 0 - B ## 0*C ## 2;		\
		ans ## 2 = B ## 0*C ## 1 - B ## 1*C ## 0
#define     SWAP(a,b)    {temp=(a); (a) = (b); (b) = temp;}

struct loc_EleNumber: public std::binary_function<ELEMENT, int, bool>
{
  bool operator () ( const ELEMENT &ele, const int &number ) const {
    return ele.Id == number;
    }
};

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//REMEMBER TO CHANGE THE NEXT LINE IF YOU CHANGE THE ACCEPTANCE RULE!
//CURRENT SETUP: swap is always accepted for boltzmann.
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Utilities::Utilities(CoordInfoVecs& coordInfoVecs,
GeneralParams& generalParams) {
            
    /*int nnode = generalParams.maxNodeCount;
    std::vector<bool> boundary_node_temp(nnode,false);
    for (int i = 0; i < nnode; i++){
        if (coordInfoVecs.edges2Triangles_1[i] == coordInfoVecs.edges2Triangles_2[i]){
            boundary_node_temp[coordInfoVecs.edges2Nodes_1[i]] = true;
            boundary_node_temp[coordInfoVecs.edges2Nodes_2[i]] = true;
        }
    }
    
    //This creates a int vector whose length equals to number of nodes.
    //The initial mesh has every node paired with 6 neighboring nodes. 
    //During the simulation, the number will be updated accordingly. Therefore this has to be moved
    //to another location to avoid re-initialization every time Edgeswap is called.

    std::vector<int> nndata_temp(nnode, 6);
    
    boundary_node = boundary_node_temp;
    nndata = nndata_temp;*/
};

//This function aims to have a fixed range or neighborhood around the tip be weakened uniformly
void Utilities::gradient_weakening_update_host_vecs_tip(double sigma,
    //double max_height_index,
    double max_height_x,
    double max_height_y,
    double max_height_z,
    double distance_to_boundary,
    double distance_uniform_weak,
    GeneralParams& generalParams,
    CoordInfoVecs& coordInfoVecs,
    HostSetInfoVecs& hostSetInfoVecs){

        double pi = 3.1415927;
    /* Scaling by gaussian distribution in the form of 
    (1/sqrt(2*pi*sigma^2))*Exp(-(d/|d|)^2/(sigma^2))/(1/sqrt(2*pi*sigma^2))
    */
   double scale;
   //double tip_threshold = 2.05*generalParams.Rmin;
   bool scale_need = false;
    if (pi < 0){
    //if (sigma == INT_MAX){
        //  for (int i = 0; i < coordInfoVecs.num_edges; i++){
        //     if (hostSetInfoVecs.edges_in_upperhem[i] != 1){
        //         hostSetInfoVecs.scaling_per_edge[i] = 0.0;
        //         continue;
        //     }
        //     else{
        //         hostSetInfoVecs.scaling_per_edge[i] = 1.0;
        //     }
        //  }
    }
    else{
        for (int i = 0; i < coordInfoVecs.num_edges; i++){
            int v1 = hostSetInfoVecs.edges2Nodes_1[i];
            int v2 = hostSetInfoVecs.edges2Nodes_2[i];
            if (hostSetInfoVecs.edges2Nodes_1[i] > (INT_MAX-100) || hostSetInfoVecs.edges2Nodes_1[i] < 0){
                hostSetInfoVecs.scaling_per_edge[i] = -INT_MAX;
                continue;
            }
            else if (hostSetInfoVecs.edges2Nodes_2[i] > (INT_MAX-100) || hostSetInfoVecs.edges2Nodes_2[i] < 0){
                hostSetInfoVecs.scaling_per_edge[i] = -INT_MAX;
                continue;
            }
            double avg_z = (hostSetInfoVecs.nodeLocZ[v1] + hostSetInfoVecs.nodeLocZ[v2])/2.0;
            if (hostSetInfoVecs.edges_in_upperhem[i] != 1 && hostSetInfoVecs.edges_in_upperhem[i] != 0){
                if (avg_z < generalParams.boundary_z){
                    hostSetInfoVecs.scaling_per_edge[i] = 1.0;
                    continue;
                }
                else{
                    double avg_x = (hostSetInfoVecs.nodeLocX[v1] + hostSetInfoVecs.nodeLocX[v2])/2.0;
                    double avg_y = (hostSetInfoVecs.nodeLocY[v1] + hostSetInfoVecs.nodeLocY[v2])/2.0;
                    double avg_z = (hostSetInfoVecs.nodeLocZ[v1] + hostSetInfoVecs.nodeLocZ[v2])/2.0;
                    double dtt = sqrt((max_height_x - avg_x)*(max_height_x - avg_x) +
                                                (max_height_y - avg_y)*(max_height_y - avg_y) +
                                                (max_height_z - avg_z)*(max_height_z - avg_z));
                    //double dtt = sqrt((hostSetInfoVecs.nodeLocX[max_height_index] - avg_x)*(hostSetInfoVecs.nodeLocX[max_height_index] - avg_x) +
                      //                          (hostSetInfoVecs.nodeLocY[max_height_index] - avg_y)*(hostSetInfoVecs.nodeLocY[max_height_index] - avg_y) +
                        //                        (hostSetInfoVecs.nodeLocZ[max_height_index] - avg_z)*(hostSetInfoVecs.nodeLocZ[max_height_index] - avg_z)); //dtt := distance to tip
                    //scale = (1.0/sqrt(2.0*pi*sigma*sigma))*exp(-(dtt/distance_to_boundary)*(dtt/distance_to_boundary)/(sigma*sigma));///(1.0/sqrt(2.0*pi*sigma*sigma));
                    //double dtt = sqrt((hostSetInfoVecs.nodeLocZ[max_height_index] - avg_z)*(hostSetInfoVecs.nodeLocZ[max_height_index] - avg_z));
                    scale = ((dtt-distance_uniform_weak)/(distance_to_boundary-distance_uniform_weak));
                    
                    if ((distance_to_boundary-distance_uniform_weak)<0.0){
                        scale = 0.0;
                    }
                    else if (scale > 1.0){
                        scale = 1.0;
                    }
                    else if (scale < 0.0){
                        scale = 0.0;
                    }
                    else{}
                    scale_need = true;
                }
            }
           
            
           

            if (generalParams.edges_in_upperhem[i] == 1 || generalParams.edges_in_upperhem[i] == 0){//(generalParams.nodes_in_upperhem[v1] == 1 && generalParams.nodes_in_upperhem[v2] == 1){
                double avg_x = (hostSetInfoVecs.nodeLocX[v1] + hostSetInfoVecs.nodeLocX[v2])/2.0;
                double avg_y = (hostSetInfoVecs.nodeLocY[v1] + hostSetInfoVecs.nodeLocY[v2])/2.0;
                double avg_z = (hostSetInfoVecs.nodeLocZ[v1] + hostSetInfoVecs.nodeLocZ[v2])/2.0;
                double dtt = sqrt((max_height_x - avg_x)*(max_height_x - avg_x) +
                                                (max_height_y - avg_y)*(max_height_y - avg_y) +
                                                (max_height_z - avg_z)*(max_height_z - avg_z));
                //double dtt = sqrt((hostSetInfoVecs.nodeLocX[max_height_index] - avg_x)*(hostSetInfoVecs.nodeLocX[max_height_index] - avg_x) +
                  //                          (hostSetInfoVecs.nodeLocY[max_height_index] - avg_y)*(hostSetInfoVecs.nodeLocY[max_height_index] - avg_y) +
                    //                        (hostSetInfoVecs.nodeLocZ[max_height_index] - avg_z)*(hostSetInfoVecs.nodeLocZ[max_height_index] - avg_z)); //dtt := distance to tip
                //scale = (1.0/sqrt(2.0*pi*sigma*sigma))*exp(-(dtt/distance_to_boundary)*(dtt/distance_to_boundary)/(sigma*sigma));///(1.0/sqrt(2.0*pi*sigma*sigma));
                //double dtt = sqrt((hostSetInfoVecs.nodeLocZ[max_height_index] - avg_z)*(hostSetInfoVecs.nodeLocZ[max_height_index] - avg_z));
                scale = ((dtt-distance_uniform_weak)/(distance_to_boundary-distance_uniform_weak));
                
                
                //if (dtt < tip_threshold){
                //    scale = 1.0;
                //}
                //else {
                //    scale = ((dtt - tip_threshold)/(distance_to_boundary - tip_threshold));
                //}
                if ((distance_to_boundary-distance_uniform_weak) < 0.0){
                    scale = 0.0;
                }
                else if (scale > 1.0){
                    scale = 1.0;
                }
                else if (scale < 0.0){
                    scale = 0.0;
                }
                else{}
                //std::cout<<"dtt/distance_to_boundary = "<<dtt/distance_to_boundary<<std::endl;
            // coordInfoVecs.scaling_per_edge[i] = scale;
                scale_need = true;
            }
            /*else if (generalParams.nodes_in_upperhem[v1] == 1 && generalParams.nodes_in_upperhem[v2] == 0){
                double avg_x = (hostSetInfoVecs.nodeLocX[v1] + hostSetInfoVecs.nodeLocX[v2])/2.0;
                double avg_y = (hostSetInfoVecs.nodeLocY[v1] + hostSetInfoVecs.nodeLocY[v2])/2.0;
                double avg_z = (hostSetInfoVecs.nodeLocZ[v1] + hostSetInfoVecs.nodeLocZ[v2])/2.0;
                double dtt = sqrt((hostSetInfoVecs.nodeLocX[max_height_index] - avg_x)*(hostSetInfoVecs.nodeLocX[max_height_index] - avg_x) +
                                            (hostSetInfoVecs.nodeLocY[max_height_index] - avg_y)*(hostSetInfoVecs.nodeLocY[max_height_index] - avg_y) +
                                            (hostSetInfoVecs.nodeLocZ[max_height_index] - avg_z)*(hostSetInfoVecs.nodeLocZ[max_height_index] - avg_z)); //dtt := distance to tip
                //scale = (1.0/sqrt(2.0*pi*sigma*sigma))*exp(-(dtt/distance_to_boundary)*(dtt/distance_to_boundary)/(sigma*sigma));///(1.0/sqrt(2.0*pi*sigma*sigma));
                scale = (dtt/distance_to_boundary);
                if (scale > 1.0){
                    scale = 1.0;
                }
                else{}
                //std::cout<<"dtt/distance_to_boundary = "<<dtt/distance_to_boundary<<std::endl;
                //coordInfoVecs.scaling_per_edge[i] = scale;
                scale_need = true;
            }
            else if (generalParams.nodes_in_upperhem[v1] == 0 && generalParams.nodes_in_upperhem[v2] == 1){
                double avg_x = (hostSetInfoVecs.nodeLocX[v1] + hostSetInfoVecs.nodeLocX[v2])/2.0;
                double avg_y = (hostSetInfoVecs.nodeLocY[v1] + hostSetInfoVecs.nodeLocY[v2])/2.0;
                double avg_z = (hostSetInfoVecs.nodeLocZ[v1] + hostSetInfoVecs.nodeLocZ[v2])/2.0;
                double dtt = sqrt((hostSetInfoVecs.nodeLocX[max_height_index] - avg_x)*(hostSetInfoVecs.nodeLocX[max_height_index] - avg_x) +
                                            (hostSetInfoVecs.nodeLocY[max_height_index] - avg_y)*(hostSetInfoVecs.nodeLocY[max_height_index] - avg_y) +
                                            (hostSetInfoVecs.nodeLocZ[max_height_index] - avg_z)*(hostSetInfoVecs.nodeLocZ[max_height_index] - avg_z)); //dtt := distance to tip
                //scale = (1.0/sqrt(2.0*pi*sigma*sigma))*exp(-(dtt/distance_to_boundary)*(dtt/distance_to_boundary)/(sigma*sigma));///(1.0/sqrt(2.0*pi*sigma*sigma));
                scale = (dtt/distance_to_boundary);
                if (scale > 1.0){
                    scale = 1.0;
                }
                else{}
                //std::cout<<"dtt/distance_to_boundary = "<<dtt/distance_to_boundary<<std::endl;
                //coordInfoVecs.scaling_per_edge[i] = scale;
                scale_need = true;
            }
            else if (generalParams.nodes_in_upperhem[v1] == 0 && generalParams.nodes_in_upperhem[v2] == 0){
                double avg_x = (hostSetInfoVecs.nodeLocX[v1] + hostSetInfoVecs.nodeLocX[v2])/2.0;
                double avg_y = (hostSetInfoVecs.nodeLocY[v1] + hostSetInfoVecs.nodeLocY[v2])/2.0;
                double avg_z = (hostSetInfoVecs.nodeLocZ[v1] + hostSetInfoVecs.nodeLocZ[v2])/2.0;
                double dtt = sqrt((hostSetInfoVecs.nodeLocX[max_height_index] - avg_x)*(hostSetInfoVecs.nodeLocX[max_height_index] - avg_x) +
                                            (hostSetInfoVecs.nodeLocY[max_height_index] - avg_y)*(hostSetInfoVecs.nodeLocY[max_height_index] - avg_y) +
                                            (hostSetInfoVecs.nodeLocZ[max_height_index] - avg_z)*(hostSetInfoVecs.nodeLocZ[max_height_index] - avg_z)); //dtt := distance to tip
               // scale = (1.0/sqrt(2.0*pi*sigma*sigma))*exp(-(dtt/distance_to_boundary)*(dtt/distance_to_boundary)/(sigma*sigma));///(1.0/sqrt(2.0*pi*sigma*sigma));
                scale = (dtt/distance_to_boundary);
                if (scale > 1.0){
                    scale = 1.0;
                }
                else{}
                //std::cout<<"dtt/distance_to_boundary = "<<dtt/distance_to_boundary<<std::endl;
                //coordInfoVecs.scaling_per_edge[i] = 1.0;//scale;
                scale_need = true;
            }*/
            else{
                //coordInfoVecs.scaling_per_edge[i] = 0.0;
            }

            if (scale_need == true){
                hostSetInfoVecs.scaling_per_edge[i] = scale;
            }
            else{
                hostSetInfoVecs.scaling_per_edge[i] = 1.0;
            }
        }
    }
};

//This function is for hill function type weakening purpose
void Utilities::gradient_weakening_update_host_vecs(double sigma,
    //double max_height_index,
    double max_height_x,
    double max_height_y,
    double max_height_z,
    double distance_to_boundary,
    double distance_to_boundary_max,
    GeneralParams& generalParams,
    CoordInfoVecs& coordInfoVecs,
    HostSetInfoVecs& hostSetInfoVecs){

        double pi = 3.1415927;
    /* Scaling by gaussian distribution in the form of 
    (1/sqrt(2*pi*sigma^2))*Exp(-(d/|d|)^2/(sigma^2))/(1/sqrt(2*pi*sigma^2))
    */
   double scale;
   //double tip_threshold = 2.05*generalParams.Rmin;
   bool scale_need = false;
    if (pi < 0){
    //if (sigma == INT_MAX){
        //  for (int i = 0; i < coordInfoVecs.num_edges; i++){
        //     if (hostSetInfoVecs.edges_in_upperhem[i] != 1){
        //         hostSetInfoVecs.scaling_per_edge[i] = 0.0;
        //         continue;
        //     }
        //     else{
        //         hostSetInfoVecs.scaling_per_edge[i] = 1.0;
        //     }
        //  }
    }
    else{
        for (int i = 0; i < coordInfoVecs.num_edges; i++){
            int v1 = hostSetInfoVecs.edges2Nodes_1[i];
            int v2 = hostSetInfoVecs.edges2Nodes_2[i];
            if (hostSetInfoVecs.edges2Nodes_1[i] == INT_MAX || hostSetInfoVecs.edges2Nodes_1[i] < 0){
                hostSetInfoVecs.scaling_per_edge[i] = -INT_MAX;
                continue;
            }
            else if (hostSetInfoVecs.edges2Nodes_2[i] == INT_MAX || hostSetInfoVecs.edges2Nodes_2[i] < 0){
                hostSetInfoVecs.scaling_per_edge[i] = -INT_MAX;
                continue;
            }
            double avg_z = (hostSetInfoVecs.nodeLocZ[v1] + hostSetInfoVecs.nodeLocZ[v2])/2.0;
            if (hostSetInfoVecs.edges_in_upperhem[i] != 1 || hostSetInfoVecs.edges_in_upperhem[i] != 0){
                if (avg_z < generalParams.boundary_z){
                    hostSetInfoVecs.scaling_per_edge[i] = 1.0;
                    continue;
                }
                else{
                    double avg_x = (hostSetInfoVecs.nodeLocX[v1] + hostSetInfoVecs.nodeLocX[v2])/2.0;
                    double avg_y = (hostSetInfoVecs.nodeLocY[v1] + hostSetInfoVecs.nodeLocY[v2])/2.0;
                    double avg_z = (hostSetInfoVecs.nodeLocZ[v1] + hostSetInfoVecs.nodeLocZ[v2])/2.0;
                    double dtt = sqrt((max_height_x - avg_x)*(max_height_x - avg_x) +
                                                (max_height_y - avg_y)*(max_height_y - avg_y) +
                                                (max_height_z - avg_z)*(max_height_z - avg_z));
                    //double dtt = sqrt((hostSetInfoVecs.nodeLocX[max_height_index] - avg_x)*(hostSetInfoVecs.nodeLocX[max_height_index] - avg_x) +
                      //                          (hostSetInfoVecs.nodeLocY[max_height_index] - avg_y)*(hostSetInfoVecs.nodeLocY[max_height_index] - avg_y) +
                        //                        (hostSetInfoVecs.nodeLocZ[max_height_index] - avg_z)*(hostSetInfoVecs.nodeLocZ[max_height_index] - avg_z)); //dtt := distance to tip
                    //scale = (1.0/sqrt(2.0*pi*sigma*sigma))*exp(-(dtt/distance_to_boundary)*(dtt/distance_to_boundary)/(sigma*sigma));///(1.0/sqrt(2.0*pi*sigma*sigma));
                    //double dtt = sqrt((hostSetInfoVecs.nodeLocZ[max_height_index] - avg_z)*(hostSetInfoVecs.nodeLocZ[max_height_index] - avg_z));
                    scale = (dtt/distance_to_boundary_max);
                    
                    if (scale > 1.0){
                        scale = 1.0;
                    }
                    else{}
                    scale_need = true;
                }
            }
           
            
           

            if (generalParams.edges_in_upperhem[i] == 1 || generalParams.edges_in_upperhem[i] == 0){//(generalParams.nodes_in_upperhem[v1] == 1 && generalParams.nodes_in_upperhem[v2] == 1){
                double avg_x = (hostSetInfoVecs.nodeLocX[v1] + hostSetInfoVecs.nodeLocX[v2])/2.0;
                double avg_y = (hostSetInfoVecs.nodeLocY[v1] + hostSetInfoVecs.nodeLocY[v2])/2.0;
                double avg_z = (hostSetInfoVecs.nodeLocZ[v1] + hostSetInfoVecs.nodeLocZ[v2])/2.0;
                double dtt = sqrt((max_height_x - avg_x)*(max_height_x - avg_x) +
                                                (max_height_y - avg_y)*(max_height_y - avg_y) +
                                                (max_height_z - avg_z)*(max_height_z - avg_z));
                //double dtt = sqrt((hostSetInfoVecs.nodeLocX[max_height_index] - avg_x)*(hostSetInfoVecs.nodeLocX[max_height_index] - avg_x) +
                  //                          (hostSetInfoVecs.nodeLocY[max_height_index] - avg_y)*(hostSetInfoVecs.nodeLocY[max_height_index] - avg_y) +
                    //                        (hostSetInfoVecs.nodeLocZ[max_height_index] - avg_z)*(hostSetInfoVecs.nodeLocZ[max_height_index] - avg_z)); //dtt := distance to tip
                //scale = (1.0/sqrt(2.0*pi*sigma*sigma))*exp(-(dtt/distance_to_boundary)*(dtt/distance_to_boundary)/(sigma*sigma));///(1.0/sqrt(2.0*pi*sigma*sigma));
                //double dtt = sqrt((hostSetInfoVecs.nodeLocZ[max_height_index] - avg_z)*(hostSetInfoVecs.nodeLocZ[max_height_index] - avg_z));
                scale = (dtt/distance_to_boundary_max);
                
                
                //if (dtt < tip_threshold){
                //    scale = 1.0;
                //}
                //else {
                //    scale = ((dtt - tip_threshold)/(distance_to_boundary - tip_threshold));
                //}

                if (scale > 1.0){
                    scale = 1.0;
                }
                else{}
                //std::cout<<"dtt/distance_to_boundary = "<<dtt/distance_to_boundary<<std::endl;
            // coordInfoVecs.scaling_per_edge[i] = scale;
                scale_need = true;
            }
            /*else if (generalParams.nodes_in_upperhem[v1] == 1 && generalParams.nodes_in_upperhem[v2] == 0){
                double avg_x = (hostSetInfoVecs.nodeLocX[v1] + hostSetInfoVecs.nodeLocX[v2])/2.0;
                double avg_y = (hostSetInfoVecs.nodeLocY[v1] + hostSetInfoVecs.nodeLocY[v2])/2.0;
                double avg_z = (hostSetInfoVecs.nodeLocZ[v1] + hostSetInfoVecs.nodeLocZ[v2])/2.0;
                double dtt = sqrt((hostSetInfoVecs.nodeLocX[max_height_index] - avg_x)*(hostSetInfoVecs.nodeLocX[max_height_index] - avg_x) +
                                            (hostSetInfoVecs.nodeLocY[max_height_index] - avg_y)*(hostSetInfoVecs.nodeLocY[max_height_index] - avg_y) +
                                            (hostSetInfoVecs.nodeLocZ[max_height_index] - avg_z)*(hostSetInfoVecs.nodeLocZ[max_height_index] - avg_z)); //dtt := distance to tip
                //scale = (1.0/sqrt(2.0*pi*sigma*sigma))*exp(-(dtt/distance_to_boundary)*(dtt/distance_to_boundary)/(sigma*sigma));///(1.0/sqrt(2.0*pi*sigma*sigma));
                scale = (dtt/distance_to_boundary);
                if (scale > 1.0){
                    scale = 1.0;
                }
                else{}
                //std::cout<<"dtt/distance_to_boundary = "<<dtt/distance_to_boundary<<std::endl;
                //coordInfoVecs.scaling_per_edge[i] = scale;
                scale_need = true;
            }
            else if (generalParams.nodes_in_upperhem[v1] == 0 && generalParams.nodes_in_upperhem[v2] == 1){
                double avg_x = (hostSetInfoVecs.nodeLocX[v1] + hostSetInfoVecs.nodeLocX[v2])/2.0;
                double avg_y = (hostSetInfoVecs.nodeLocY[v1] + hostSetInfoVecs.nodeLocY[v2])/2.0;
                double avg_z = (hostSetInfoVecs.nodeLocZ[v1] + hostSetInfoVecs.nodeLocZ[v2])/2.0;
                double dtt = sqrt((hostSetInfoVecs.nodeLocX[max_height_index] - avg_x)*(hostSetInfoVecs.nodeLocX[max_height_index] - avg_x) +
                                            (hostSetInfoVecs.nodeLocY[max_height_index] - avg_y)*(hostSetInfoVecs.nodeLocY[max_height_index] - avg_y) +
                                            (hostSetInfoVecs.nodeLocZ[max_height_index] - avg_z)*(hostSetInfoVecs.nodeLocZ[max_height_index] - avg_z)); //dtt := distance to tip
                //scale = (1.0/sqrt(2.0*pi*sigma*sigma))*exp(-(dtt/distance_to_boundary)*(dtt/distance_to_boundary)/(sigma*sigma));///(1.0/sqrt(2.0*pi*sigma*sigma));
                scale = (dtt/distance_to_boundary);
                if (scale > 1.0){
                    scale = 1.0;
                }
                else{}
                //std::cout<<"dtt/distance_to_boundary = "<<dtt/distance_to_boundary<<std::endl;
                //coordInfoVecs.scaling_per_edge[i] = scale;
                scale_need = true;
            }
            else if (generalParams.nodes_in_upperhem[v1] == 0 && generalParams.nodes_in_upperhem[v2] == 0){
                double avg_x = (hostSetInfoVecs.nodeLocX[v1] + hostSetInfoVecs.nodeLocX[v2])/2.0;
                double avg_y = (hostSetInfoVecs.nodeLocY[v1] + hostSetInfoVecs.nodeLocY[v2])/2.0;
                double avg_z = (hostSetInfoVecs.nodeLocZ[v1] + hostSetInfoVecs.nodeLocZ[v2])/2.0;
                double dtt = sqrt((hostSetInfoVecs.nodeLocX[max_height_index] - avg_x)*(hostSetInfoVecs.nodeLocX[max_height_index] - avg_x) +
                                            (hostSetInfoVecs.nodeLocY[max_height_index] - avg_y)*(hostSetInfoVecs.nodeLocY[max_height_index] - avg_y) +
                                            (hostSetInfoVecs.nodeLocZ[max_height_index] - avg_z)*(hostSetInfoVecs.nodeLocZ[max_height_index] - avg_z)); //dtt := distance to tip
               // scale = (1.0/sqrt(2.0*pi*sigma*sigma))*exp(-(dtt/distance_to_boundary)*(dtt/distance_to_boundary)/(sigma*sigma));///(1.0/sqrt(2.0*pi*sigma*sigma));
                scale = (dtt/distance_to_boundary);
                if (scale > 1.0){
                    scale = 1.0;
                }
                else{}
                //std::cout<<"dtt/distance_to_boundary = "<<dtt/distance_to_boundary<<std::endl;
                //coordInfoVecs.scaling_per_edge[i] = 1.0;//scale;
                scale_need = true;
            }*/
            else{
                //coordInfoVecs.scaling_per_edge[i] = 0.0;
            }

            if (scale_need == true){
                hostSetInfoVecs.scaling_per_edge[i] = scale;
            }
            else{
                hostSetInfoVecs.scaling_per_edge[i] = 1.0;
            }
        }
    }
};

//This function is for chemical signaling driven weakening purpose
void Utilities::chemSignal_weakening_update_host_vecs(double max_conc,
    GeneralParams& generalParams,
    CoordInfoVecs& coordInfoVecs,
    HostSetInfoVecs& hostSetInfoVecs){

    double pi = 3.1415927;
    double scale;
    bool scale_need = false;
    if (pi < 0){
    }
    else{
        for (int i = 0; i < coordInfoVecs.num_edges; i++){
            int v1 = hostSetInfoVecs.edges2Nodes_1[i];
            int v2 = hostSetInfoVecs.edges2Nodes_2[i];
            if (hostSetInfoVecs.edges2Nodes_1[i] == INT_MAX || hostSetInfoVecs.edges2Nodes_1[i] < 0){
                hostSetInfoVecs.scaling_per_edge[i] = -INT_MAX;
                continue;
            }
            else if (hostSetInfoVecs.edges2Nodes_2[i] == INT_MAX || hostSetInfoVecs.edges2Nodes_2[i] < 0){
                hostSetInfoVecs.scaling_per_edge[i] = -INT_MAX;
                continue;
            }
            double avg_z = (hostSetInfoVecs.nodeLocZ[v1] + hostSetInfoVecs.nodeLocZ[v2])/2.0;
            if (hostSetInfoVecs.edges_in_upperhem[i] != 1 || hostSetInfoVecs.edges_in_upperhem[i] != 0){
                if (avg_z < generalParams.boundary_z){
                    hostSetInfoVecs.scaling_per_edge[i] = 1.0;
                    continue;
                }
                else{
                    double tri1 = coordInfoVecs.edges2Triangles_1[i];
					double tri2 = coordInfoVecs.edges2Triangles_2[i];
					double avg_conc = (coordInfoVecs.soln_per_triangle[tri1] + coordInfoVecs.soln_per_triangle[tri2])/2.0;
					
                    scale = (avg_conc/max_conc);
                    
                    if (scale > 1.0){
                        scale = 1.0;
                    }
                    scale_need = true;
                }
            }

            if (generalParams.edges_in_upperhem[i] == 1 || generalParams.edges_in_upperhem[i] == 0){//(generalParams.nodes_in_upperhem[v1] == 1 && generalParams.nodes_in_upperhem[v2] == 1){
                double tri1 = coordInfoVecs.edges2Triangles_1[i];
                double tri2 = coordInfoVecs.edges2Triangles_2[i];
                double avg_conc = (coordInfoVecs.soln_per_triangle[tri1] + coordInfoVecs.soln_per_triangle[tri2])/2.0;
                
                scale = (avg_conc/max_conc);
                
                
                if (scale > 1.0){
                    scale = 1.0;
                }
                else{}
                scale_need = true;
            }
            else{
                //coordInfoVecs.scaling_per_edge[i] = 0.0;
            }

            if (scale_need == true){
                hostSetInfoVecs.scaling_per_edge[i] = scale;
            }
            else{
                hostSetInfoVecs.scaling_per_edge[i] = 1.0;
            }
        }
    }
};

int Utilities::surfaceNormal_device_vecs(
    int inode,
    CoordInfoVecs& coordInfoVecs,
    GeneralParams& generalParams
){
    double current_normalX = 0.0;
    double current_normalY = 0.0;
    double current_normalZ = 0.0;
    int neighbor;
    double AREA;
    double scaled_forceX, scaled_forceY, scaled_forceZ;
    for (int j = 0; j < 9; j ++){
        if (j == 0 && (coordInfoVecs.nodes2Triangles_1[inode] >= 0)){neighbor = coordInfoVecs.nodes2Triangles_1[inode];}
        //else{continue;}
        else if (j == 1 && (coordInfoVecs.nodes2Triangles_2[inode] >= 0)){neighbor = coordInfoVecs.nodes2Triangles_2[inode];}
        //else{continue;}
        else if (j == 2 && (coordInfoVecs.nodes2Triangles_3[inode] >= 0)){neighbor = coordInfoVecs.nodes2Triangles_3[inode];}
        //else{continue;}
        else if (j == 3 && (coordInfoVecs.nodes2Triangles_4[inode] >= 0)){neighbor = coordInfoVecs.nodes2Triangles_4[inode];}
        //else{continue;}
        else if (j == 4 && (coordInfoVecs.nodes2Triangles_5[inode] >= 0)){neighbor = coordInfoVecs.nodes2Triangles_5[inode];}
        //else{continue;}
        else if (j == 5 && (coordInfoVecs.nodes2Triangles_6[inode] >= 0)){neighbor = coordInfoVecs.nodes2Triangles_6[inode];}
        //else{continue;}
        else if (j == 6 && (coordInfoVecs.nodes2Triangles_7[inode] >= 0)){neighbor = coordInfoVecs.nodes2Triangles_7[inode];}
        //else{continue;}
        else if (j == 7 && (coordInfoVecs.nodes2Triangles_8[inode] >= 0)){neighbor = coordInfoVecs.nodes2Triangles_8[inode];}
        //else{continue;}
        else if (j == 8 && (coordInfoVecs.nodes2Triangles_9[inode] >= 0)){neighbor = coordInfoVecs.nodes2Triangles_9[inode];}
        //else{continue;}
       // else if (j == 9 && (coordInfoVecs.nodes2Triangles_10[inode] >= 0)){neighbor = coordInfoVecs.nodes2Triangles_10[inode];}
        //else{continue;}
       // else if (j == 10 && (coordInfoVecs.nodes2Triangles_11[inode] >= 0)){neighbor = coordInfoVecs.nodes2Triangles_11[inode];}
        //else{continue;}
        //else if (j == 11 && (coordInfoVecs.nodes2Triangles_12[inode] >= 0)){neighbor = coordInfoVecs.nodes2Triangles_12[inode];}
        else{continue;}
       //for (int i = 0; i < 12; i++){
        
        int node1 = coordInfoVecs. triangles2Nodes_1[neighbor];
        int node2 = coordInfoVecs. triangles2Nodes_2[neighbor];
        int node3 = coordInfoVecs. triangles2Nodes_3[neighbor];

        double vec1x = coordInfoVecs.nodeLocX[node2] - coordInfoVecs.nodeLocX[node1];
        double vec1y = coordInfoVecs.nodeLocY[node2] - coordInfoVecs.nodeLocY[node1];
        double vec1z = coordInfoVecs.nodeLocZ[node2] - coordInfoVecs.nodeLocZ[node1];
        double vec2x = coordInfoVecs.nodeLocX[node3] - coordInfoVecs.nodeLocX[node1];
        double vec2y = coordInfoVecs.nodeLocY[node3] - coordInfoVecs.nodeLocY[node1];
        double vec2z = coordInfoVecs.nodeLocZ[node3] - coordInfoVecs.nodeLocZ[node1];

        double normalX = vec1y*vec2z - vec2y*vec1z;
        double normalY = -(vec1x*vec2z - vec2x*vec1z);
        double normalZ = vec1x*vec2y - vec2x*vec1y;
        double normalize = sqrt((normalX*normalX)+(normalY*normalY)+(normalZ*normalZ));
        normalX = normalX/normalize;
        normalY = normalY/normalize;
        normalZ = normalZ/normalize;
        AREA = normalize/2.0;
        scaled_forceX = ((1.0/3.0)*AREA*normalX);
        scaled_forceY = ((1.0/3.0)*AREA*normalY);
        scaled_forceZ = ((1.0/3.0)*AREA*normalZ);
        //current_normalX += normalX;
        //current_normalY += normalY;
        //current_normalZ += normalZ;
        coordInfoVecs.nodeForceX[inode] += generalParams.volume_spring_constant*scaled_forceX;
        coordInfoVecs.nodeForceY[inode] += generalParams.volume_spring_constant*scaled_forceY;
        coordInfoVecs.nodeForceZ[inode] += generalParams.volume_spring_constant*scaled_forceZ;
    }

    /*double UN = sqrt((current_normalX*current_normalX) +
                     (current_normalY*current_normalY) + 
                     (current_normalZ*current_normalZ));
    current_normalX = current_normalX/UN;
    current_normalY = current_normalY/UN;
    current_normalZ = current_normalZ/UN;*/

    /*for (int j = 0; j < AREA.size(); j++){
        coordInfoVecs.nodeForceX[inode] += generalParams.volume_spring_constant*current_normalX;
        coordInfoVecs.nodeForceY[inode] += generalParams.volume_spring_constant*current_normalY;
        coordInfoVecs.nodeForceZ[inode] += generalParams.volume_spring_constant*current_normalZ;
    }*/
};

// int Utilities::growth_host_vecs(
int Utilities::growth_host_vecs(
    int iedge, 
    GeneralParams& generalParams,
    HostSetInfoVecs& hostSetInfoVecs,
    CoordInfoVecs& coordInfoVecs,
    LinearSpringInfoVecs& linearSpringInfoVecs,
    BendingTriangleInfoVecs& bendingTriangleInfoVecs,
    AreaTriangleInfoVecs& areaTriangleInfoVecs
){
    int alpha;
    // std::vector<int> alpha;

    int HEAD,TAIL;
    int H0, T0,H1,H2,T1,T2;
    int edge_start, edge_end;
    int a1, b1, c1, a2, b2, c2;
    double temp_bend = 0.0;
    double linear_spring_constant;
    double bend_spring_constant;
    double vol_0, vol_1;
    double P0x_vol1, P0y_vol1, P0z_vol1, P0x_vol2, P0y_vol2, P0z_vol2;
    double N1x_vol, N1y_vol, N1z_vol, N2x_vol, N2y_vol, N2z_vol;
    bool GROWTH_ACCEPTED = false;
      
    if ( hostSetInfoVecs.edges2Triangles_1[iedge] != hostSetInfoVecs.edges2Triangles_2[iedge]){
        H0 = hostSetInfoVecs.edges2Triangles_1[iedge];//index of the 1st triangle to i-th edge
        generalParams.triangle_undergoing_growth.push_back(H0);
        //std::cout<<"H0 = "<<H0<<std::endl;
        T0 = hostSetInfoVecs.edges2Triangles_2[iedge];//index of the 2nd triangle to i-th edge
        generalParams.triangle_undergoing_growth.push_back(T0);
        //std::cout<<"T0 = "<<T0<<std::endl;
        edge_start = hostSetInfoVecs.edges2Nodes_1[iedge];//index of the 1st node of i-th edge
        edge_end = hostSetInfoVecs.edges2Nodes_2[iedge];//index of the 2nd node of i-th edge

        a1 = hostSetInfoVecs.triangles2Edges_1[H0];//index of the 1st node of triangle H0
        //std::cout<<"a1 = "<<a1<<std::endl;
        b1 = hostSetInfoVecs.triangles2Edges_2[H0];//index of the 2nd node of triangle H0
        //std::cout<<"b1 = "<<b1<<std::endl;
        c1 = hostSetInfoVecs.triangles2Edges_3[H0];//index of the 3rd node of triangle H0
        //std::cout<<"c1 = "<<c1<<std::endl;
        
        a2 = hostSetInfoVecs.triangles2Edges_1[T0];//index of the 1st node of triangle T0
        //std::cout<<"a2 = "<<a2<<std::endl;
        b2 = hostSetInfoVecs.triangles2Edges_2[T0];
        //std::cout<<"b2 = "<<b2<<std::endl;
        c2 = hostSetInfoVecs.triangles2Edges_3[T0];
        //std::cout<<"c2 = "<<c2<<std::endl;
        
        //Now we identify the edge indices associated with the small subsystem.
        //This gives us the indices for H1, H2, T1, T2 (see the figure below).
        if (a1 != iedge && hostSetInfoVecs.edges2Nodes_1[a1] == edge_start){H1 = a1;}
        else if (a1 != iedge && hostSetInfoVecs.edges2Nodes_2[a1] == edge_start){H1 = a1;}
        else if (a1 != iedge && hostSetInfoVecs.edges2Nodes_1[a1] == edge_end){H2 = a1;}
        else if (a1 != iedge && hostSetInfoVecs.edges2Nodes_2[a1] == edge_end){H2 = a1;}

        if (b1 != iedge && hostSetInfoVecs.edges2Nodes_1[b1] == edge_start){H1 = b1;}
        else if (b1 != iedge && hostSetInfoVecs.edges2Nodes_2[b1] == edge_start){H1 = b1;}
        else if (b1 != iedge && hostSetInfoVecs.edges2Nodes_1[b1] == edge_end){H2 = b1;}
        else if (b1 != iedge && hostSetInfoVecs.edges2Nodes_2[b1] == edge_end){H2 = b1;}

        if (c1 != iedge && hostSetInfoVecs.edges2Nodes_1[c1] == edge_start){H1 = c1;}
        else if (c1 != iedge && hostSetInfoVecs.edges2Nodes_2[c1] == edge_start){H1 = c1;}
        else if (c1 != iedge && hostSetInfoVecs.edges2Nodes_1[c1] == edge_end){H2 = c1;}
        else if (c1 != iedge && hostSetInfoVecs.edges2Nodes_2[c1] == edge_end){H2 = c1;}
        
        if (a2 != iedge && hostSetInfoVecs.edges2Nodes_1[a2] == edge_start){T1 = a2;}
        else if (a2 != iedge && hostSetInfoVecs.edges2Nodes_2[a2] == edge_start){T1 = a2;}
        else if (a2 != iedge && hostSetInfoVecs.edges2Nodes_1[a2] == edge_end){T2 = a2;}
        else if (a2 != iedge && hostSetInfoVecs.edges2Nodes_2[a2] == edge_end){T2 = a2;}

        if (b2 != iedge && hostSetInfoVecs.edges2Nodes_1[b2] == edge_start){T1 = b2;}
        else if (b2 != iedge && hostSetInfoVecs.edges2Nodes_2[b2] == edge_start){T1 = b2;}
        else if (b2 != iedge && hostSetInfoVecs.edges2Nodes_1[b2] == edge_end){T2 = b2;}
        else if (b2 != iedge && hostSetInfoVecs.edges2Nodes_2[b2] == edge_end){T2 = b2;}

        if (c2 != iedge && hostSetInfoVecs.edges2Nodes_1[c2] == edge_start){T1 = c2;}
        else if (c2 != iedge && hostSetInfoVecs.edges2Nodes_2[c2] == edge_start){T1 = c2;}
        else if (c2 != iedge && hostSetInfoVecs.edges2Nodes_1[c2] == edge_end){T2 = c2;}
        else if (c2 != iedge && hostSetInfoVecs.edges2Nodes_2[c2] == edge_end){T2 = c2;}

        //std::cout<<"H1 = "<<H1<<std::endl;
        //std::cout<<"H2 = "<<H2<<std::endl;
        //std::cout<<"T1 = "<<T1<<std::endl;
        //std::cout<<"T2 = "<<T2<<std::endl;

        //Now search for the associated 

        int CANDIDATE1_1 = hostSetInfoVecs.triangles2Nodes_1[H0];
        int CANDIDATE1_2 = hostSetInfoVecs.triangles2Nodes_2[H0];
        int CANDIDATE1_3 = hostSetInfoVecs.triangles2Nodes_3[H0];
        int CANDIDATE2_1 = hostSetInfoVecs.triangles2Nodes_1[T0];
        int CANDIDATE2_2 = hostSetInfoVecs.triangles2Nodes_2[T0];
        int CANDIDATE2_3 = hostSetInfoVecs.triangles2Nodes_3[T0];
        
        if ((CANDIDATE1_1 != edge_start) 
            && (CANDIDATE1_1 != edge_end)) {
            HEAD = CANDIDATE1_1;
        }
        else if ((CANDIDATE1_2 != edge_start) && (CANDIDATE1_2 != edge_end)){HEAD = CANDIDATE1_2;}
        else if (CANDIDATE1_3 != edge_start && CANDIDATE1_3 != edge_end){HEAD = CANDIDATE1_3;}
        else {std::cout<<"head not set" <<std::endl;}

        if (CANDIDATE2_1 != edge_start && CANDIDATE2_1 != edge_end){TAIL = CANDIDATE2_1;}
        else if (CANDIDATE2_2 != edge_start && CANDIDATE2_2 != edge_end){TAIL = CANDIDATE2_2;}
        else if (CANDIDATE2_3 != edge_start && CANDIDATE2_3 != edge_end){TAIL = CANDIDATE2_3;}
        else {std::cout<<"tail not set" <<std::endl;}

        bool BAD_CHOICE = false;
        for (int q = 0; q < 4; q++){
            double qq;
            if (q == 0){
                qq = edge_start;
            }
            else if (q == 1){
                qq = edge_end;
            }
            else if (q == 2){
                qq = HEAD;
            }
            else if (q == 3){
                qq = TAIL;
            }
            
            int safe_flip1 = 0;
            if (hostSetInfoVecs.nndata1[qq] >= 0){safe_flip1 += 1;        }
            if (hostSetInfoVecs.nndata2[qq] >= 0){safe_flip1 += 1;        }
            if (hostSetInfoVecs.nndata3[qq] >= 0){safe_flip1 += 1;        }
            if (hostSetInfoVecs.nndata4[qq] >= 0){safe_flip1 += 1;        }
            if (hostSetInfoVecs.nndata5[qq] >= 0){safe_flip1 += 1;        }
            if (hostSetInfoVecs.nndata6[qq] >= 0){safe_flip1 += 1;        }
            if (hostSetInfoVecs.nndata7[qq] >= 0){safe_flip1 += 1;        }
            if (hostSetInfoVecs.nndata8[qq] >= 0){safe_flip1 += 1;        }
            if (hostSetInfoVecs.nndata9[qq] >= 0){safe_flip1 += 1;        }
            //if (hostSetInfoVecs.nndata10[qq] >= 0){safe_flip1 += 1;        }
            //if (hostSetInfoVecs.nndata11[qq] >= 0){safe_flip1 += 1;        }
            //if (hostSetInfoVecs.nndata12[qq] >= 0){safe_flip1 += 1;        }

            if (q == 0){// && safe_flip1 == 4){
                //BAD_CHOICE = true;
               // std::cout<<"BAD_CHOICE = "<<BAD_CHOICE<<std::endl;
                //std::cout<<"SAFE_FLIP_start = "<<safe_flip1<<std::endl;
                //break;
            }
            else if (q == 1){// && safe_flip1 == 4){
                //BAD_CHOICE = true;
               // std::cout<<"BAD_CHOICE = "<<BAD_CHOICE<<std::endl;
        //std::cout<<"SAFE_FLIP_end = "<<safe_flip1<<std::endl;
                //break;
            }
            else if (q == 2 && safe_flip1 == generalParams.safeguardthreshold){
                BAD_CHOICE = true;
               // std::cout<<"BAD_CHOICE = "<<BAD_CHOICE<<std::endl;
        //std::cout<<"SAFE_FLIP_H = "<<safe_flip1<<std::endl;
                break;
            }
            else if (q == 3 && safe_flip1 == generalParams.safeguardthreshold){
                BAD_CHOICE = true;
               // std::cout<<"BAD_CHOICE = "<<BAD_CHOICE<<std::endl;
        //std::cout<<"SAFE_FLIP_T = "<<safe_flip1<<std::endl;
                break;
            }
        }
        

        if (BAD_CHOICE == false){//(safe_flip1 < generalParams.safeguardthreshold && safe_flip2 < generalParams.safeguardthreshold){

            int H0n1 = edge_start;//hostSetInfoVecs.triangles2Nodes_1[H0];
            int H0n2 = edge_end;//hostSetInfoVecs.triangles2Nodes_2[H0];
            int H0n3 = HEAD;//hostSetInfoVecs.triangles2Nodes_3[H0];
            int T0n1 = edge_start;//hostSetInfoVecs.triangles2Nodes_1[T0];
            int T0n2 = TAIL;//hostSetInfoVecs.triangles2Nodes_2[T0];
            int T0n3 = edge_end;//hostSetInfoVecs.triangles2Nodes_3[T0];
            double a = sqrt(pow((hostSetInfoVecs.nodeLocX[H0n2] - hostSetInfoVecs.nodeLocX[H0n1]),2.0) + 
                        pow((hostSetInfoVecs.nodeLocY[H0n2] - hostSetInfoVecs.nodeLocY[H0n1]),2.0) +
                        pow((hostSetInfoVecs.nodeLocZ[H0n2] - hostSetInfoVecs.nodeLocZ[H0n1]),2.0)
                        );
            double b = sqrt(pow((hostSetInfoVecs.nodeLocX[H0n3] - hostSetInfoVecs.nodeLocX[H0n1]),2.0) + 
                        pow((hostSetInfoVecs.nodeLocY[H0n3] - hostSetInfoVecs.nodeLocY[H0n1]),2.0) +
                        pow((hostSetInfoVecs.nodeLocZ[H0n3] - hostSetInfoVecs.nodeLocZ[H0n1]),2.0)
                        );
            double c = sqrt(pow((hostSetInfoVecs.nodeLocX[H0n3] - hostSetInfoVecs.nodeLocX[H0n2]),2.0) + 
                        pow((hostSetInfoVecs.nodeLocY[H0n3] - hostSetInfoVecs.nodeLocY[H0n2]),2.0) +
                        pow((hostSetInfoVecs.nodeLocZ[H0n3] - hostSetInfoVecs.nodeLocZ[H0n2]),2.0)
                        );
            double mean_abc = (a + b + c)/2;
            double d = sqrt(pow((hostSetInfoVecs.nodeLocX[T0n2] - hostSetInfoVecs.nodeLocX[T0n1]),2.0) + 
                        pow((hostSetInfoVecs.nodeLocY[T0n2] - hostSetInfoVecs.nodeLocY[T0n1]),2.0) +
                        pow((hostSetInfoVecs.nodeLocZ[T0n2] - hostSetInfoVecs.nodeLocZ[T0n1]),2.0)
                        );
            double e = sqrt(pow((hostSetInfoVecs.nodeLocX[T0n3] - hostSetInfoVecs.nodeLocX[T0n1]),2.0) + 
                        pow((hostSetInfoVecs.nodeLocY[T0n3] - hostSetInfoVecs.nodeLocY[T0n1]),2.0) +
                        pow((hostSetInfoVecs.nodeLocZ[T0n3] - hostSetInfoVecs.nodeLocZ[T0n1]),2.0)
                        );
            double f = sqrt(pow((hostSetInfoVecs.nodeLocX[T0n3] - hostSetInfoVecs.nodeLocX[T0n2]),2.0) + 
                        pow((hostSetInfoVecs.nodeLocY[T0n3] - hostSetInfoVecs.nodeLocY[T0n2]),2.0) +
                        pow((hostSetInfoVecs.nodeLocZ[T0n3] - hostSetInfoVecs.nodeLocZ[T0n2]),2.0)
                        );
            double mean_def = (d + e + f)/2.0;
            // double area_spring_constant_1, area_spring_constant_2;
            // if (hostSetInfoVecs.triangles_in_upperhem[H0] == 1){
            //     area_spring_constant_1 = areaTriangleInfoVecs.spring_constant_weak;
            // }
            // else if (hostSetInfoVecs.triangles_in_upperhem[H0] == 0){
            //     area_spring_constant_1 = (areaTriangleInfoVecs.spring_constant_weak + areaTriangleInfoVecs.spring_constant)/2.0;
            // }
            // else{
            //     area_spring_constant_1 = areaTriangleInfoVecs.spring_constant;
            // }
            // if (hostSetInfoVecs.triangles_in_upperhem[T0] == 1){
            //     area_spring_constant_2 = areaTriangleInfoVecs.spring_constant_weak;
            // }
            // else if (hostSetInfoVecs.triangles_in_upperhem[T0] == 0){
            //     area_spring_constant_2 = (areaTriangleInfoVecs.spring_constant_weak + areaTriangleInfoVecs.spring_constant)/2.0;
            // }
            // else{
            //     area_spring_constant_2 = areaTriangleInfoVecs.spring_constant;
            // }
            double area_H0 = sqrt(mean_abc*(mean_abc - a)*(mean_abc - b)*(mean_abc - c));
            double area_T0 = sqrt(mean_def*(mean_def - d)*(mean_def - e)*(mean_def - f));
            double avg_area = (area_H0 + area_T0)/2.0;
            // std::random_device rd;  //Will be used to obtain a seed for the random number engine
            // std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
            // std::uniform_real_distribution<> dis(0.0, 1.0);
            // double random_number = dis(gen);
            // //double random_number = 0;
            // double Edif = generalParams.insertion_energy_cost*(1.0 - (1.0/generalParams.strain_threshold)*((avg_area/areaTriangleInfoVecs.initial_area) - 1));
            // //std::cout<<Edif<<std::endl;
            // double prob = generalParams.tau*exp(-(Edif*generalParams.growth_energy_scaling)/generalParams.kT);
            // if (prob >= 1){
            //     prob = 1;
            // }
            // //std::cout<<"prob = "<<prob<<std::endl;
            // if (random_number < prob){
            //     GROWTH_ACCEPTED = true;
            // }
            if ((avg_area - areaTriangleInfoVecs.initial_area)/areaTriangleInfoVecs.initial_area >= generalParams.strain_threshold){
            GROWTH_ACCEPTED = true;
        }
      }
    //if (H0 > 20000.0){
      if (GROWTH_ACCEPTED == true){
        /*  int k = iedge;
		if (hostSetInfoVecs.edges2Nodes_1[k] == INT_MAX || hostSetInfoVecs.edges2Nodes_2[k] == INT_MAX){
			continue;
		}
		if (generalParams.boundaries_in_upperhem[k] == 1){
			continue;
		}
		if (generalParams.edges_in_upperhem[k] < 0 || generalParams.edges_in_upperhem[k] == INT_MAX){
			continue;
		}
		int iedge = k;*/
		
		
		////std::cout<<"GROWTH ERROR 2"<<std::endl;	
		int t1e1, t1e2, t1e3, t2e1, t2e2, t2e3;
        int elem1 = H0;
        int elem2 = T0;
		if (hostSetInfoVecs.triangles2Edges_1[elem1] == iedge){
			t1e1 = hostSetInfoVecs.triangles2Edges_2[elem1];
			t1e2 = hostSetInfoVecs.triangles2Edges_3[elem1];
			//t1e3 = hostSetInfoVecs.triangles2Edges_1[elem1];
		}
		else if (hostSetInfoVecs.triangles2Edges_2[elem1] == iedge){
			t1e1 = hostSetInfoVecs.triangles2Edges_3[elem1];
			t1e2 = hostSetInfoVecs.triangles2Edges_1[elem1];
			//t1e3 = hostSetInfoVecs.triangles2Edges_2[elem1];
		} 
		else if (hostSetInfoVecs.triangles2Edges_3[elem1] == iedge){
			t1e1 = hostSetInfoVecs.triangles2Edges_1[elem1];
			t1e2 = hostSetInfoVecs.triangles2Edges_2[elem1];
			//t1e3 = hostSetInfoVecs.triangles2Edges_3[elem1];
		}
		////std::cout<<"GROWTH ERROR 3"<<std::endl;	

		if (hostSetInfoVecs.triangles2Edges_1[elem2] == iedge){
			t2e1 = hostSetInfoVecs.triangles2Edges_2[elem2];
			t2e2 = hostSetInfoVecs.triangles2Edges_3[elem2];
			//t2e3 = hostSetInfoVecs.triangles2Edges_1[elem2];
		}
		else if (hostSetInfoVecs.triangles2Edges_2[elem2] == iedge){
			t2e1 = hostSetInfoVecs.triangles2Edges_3[elem2];
			t2e2 = hostSetInfoVecs.triangles2Edges_1[elem2];
			//t2e3 = hostSetInfoVecs.triangles2Edges_2[elem2];
		} 
		else if (hostSetInfoVecs.triangles2Edges_3[elem2] == iedge){
			t2e1 = hostSetInfoVecs.triangles2Edges_1[elem2];
			t2e2 = hostSetInfoVecs.triangles2Edges_2[elem2];
			//t2e3 = hostSetInfoVecs.triangles2Edges_3[elem2];
		}
		//Note that in the above assignment, t1e3 and t2e3 are the edges shared by the neighboring triangles T1 and T2.
		////std::cout<<"GROWTH ERROR 4"<<std::endl;	

		//VectorShuffleForEdgeswapLoop.push_back(t1e1);
		//VectorShuffleForEdgeswapLoop.push_back(t1e2);
		//VectorShuffleForEdgeswapLoop.push_back(t2e1);
		//VectorShuffleForEdgeswapLoop.push_back(t2e2);
		int n1, n2, n3, n4;
		
		if ((hostSetInfoVecs.edges2Nodes_1[t1e1] == hostSetInfoVecs.edges2Nodes_1[iedge]) || (hostSetInfoVecs.edges2Nodes_1[t1e1] == hostSetInfoVecs.edges2Nodes_2[iedge]) ){
			n1 = hostSetInfoVecs.edges2Nodes_1[t1e1];
			n2 = hostSetInfoVecs.edges2Nodes_2[t1e1];
			if (hostSetInfoVecs.edges2Nodes_1[t1e1] == hostSetInfoVecs.edges2Nodes_1[iedge]){
				n3 = hostSetInfoVecs.edges2Nodes_2[iedge];
			}
			else if (hostSetInfoVecs.edges2Nodes_1[t1e1] == hostSetInfoVecs.edges2Nodes_2[iedge]){
				n3 = hostSetInfoVecs.edges2Nodes_1[iedge];
			}
		}
		else if ((hostSetInfoVecs.edges2Nodes_2[t1e1] == hostSetInfoVecs.edges2Nodes_1[iedge]) || (hostSetInfoVecs.edges2Nodes_2[t1e1] == hostSetInfoVecs.edges2Nodes_2[iedge]) ){
			n1 = hostSetInfoVecs.edges2Nodes_2[t1e1];
			n2 = hostSetInfoVecs.edges2Nodes_1[t1e1];
			if (hostSetInfoVecs.edges2Nodes_2[t1e1] == hostSetInfoVecs.edges2Nodes_1[iedge]){
				n3 = hostSetInfoVecs.edges2Nodes_2[iedge];
			}
			else if (hostSetInfoVecs.edges2Nodes_2[t1e1] == hostSetInfoVecs.edges2Nodes_2[iedge]){
				n3 = hostSetInfoVecs.edges2Nodes_1[iedge];
			}
		}
		////std::cout<<"GROWTH ERROR 5"<<std::endl;	

		if (hostSetInfoVecs.edges2Nodes_1[t2e1] == hostSetInfoVecs.edges2Nodes_1[iedge] || hostSetInfoVecs.edges2Nodes_1[t2e1] == hostSetInfoVecs.edges2Nodes_2[iedge]){
			n4 = hostSetInfoVecs.edges2Nodes_2[t2e1];
		}
		else if (hostSetInfoVecs.edges2Nodes_2[t2e1] == hostSetInfoVecs.edges2Nodes_1[iedge] || hostSetInfoVecs.edges2Nodes_2[t2e1] == hostSetInfoVecs.edges2Nodes_2[iedge]){
			n4 = hostSetInfoVecs.edges2Nodes_1[t2e1];
		}
		

		//std::cout<<"n1 = "<<n1<<std::endl;
		//std::cout<<"n2 = "<<n2<<std::endl;
		//std::cout<<"n3 = "<<n3<<std::endl;
		//std::cout<<"n4 = "<<n4<<std::endl;
		//These extract the indices of vertices of the selected triangles "elem1" and "elem2". Now we have n1, n2, n3, n4 in the correct orientation (supposedly).

		////std::cout<<"GROWTH ERROR 6"<<std::endl;	
		int edgeindex, a, a1, a2, a3, temp1, temp2;
		//std::cout<<"maxNodeCount = "<< generalParams.maxNodeCount<<std::endl;
		double newx = (hostSetInfoVecs.nodeLocX[hostSetInfoVecs.edges2Nodes_1[iedge]] + hostSetInfoVecs.nodeLocX[hostSetInfoVecs.edges2Nodes_2[iedge]])/2.0;
		////std::cout<<"GROWTH ERROR 6.1"<<std::endl;	
		hostSetInfoVecs.nodeLocX[generalParams. maxNodeCount] = newx;
		////std::cout<<"GROWTH ERROR 6.2"<<std::endl;	
		double newy = (hostSetInfoVecs.nodeLocY[hostSetInfoVecs.edges2Nodes_1[iedge]] + hostSetInfoVecs.nodeLocY[hostSetInfoVecs.edges2Nodes_2[iedge]])/2.0;
		////std::cout<<"GROWTH ERROR 6.3"<<std::endl;	
		hostSetInfoVecs.nodeLocY[generalParams. maxNodeCount] = newy;
		////std::cout<<"GROWTH ERROR 6.4"<<std::endl;	
		double newz = (hostSetInfoVecs.nodeLocZ[hostSetInfoVecs.edges2Nodes_1[iedge]] + hostSetInfoVecs.nodeLocZ[hostSetInfoVecs.edges2Nodes_2[iedge]])/2.0;
		////std::cout<<"GROWTH ERROR 6.5"<<std::endl;	
		hostSetInfoVecs.nodeLocZ[generalParams. maxNodeCount] = newz;
		//These are the coordinate of the new vertex. Its index is "hostSetInfoVecs.nodeLocX.size()-1"

		//Before editing major data structures, we will update the nndata here since it is only affected by the addition of new nodes.

		//int NODESIZE= generalParams.maxNodeCount;//hostSetInfoVecs.nodeLocX.size();
		////std::cout<<"GROWTH ERROR 7"<<std::endl;			
		hostSetInfoVecs.triangles2Nodes_1[coordInfoVecs.num_triangles] = n1;
		hostSetInfoVecs.triangles2Nodes_2[coordInfoVecs.num_triangles] = n2;
		hostSetInfoVecs.triangles2Nodes_3[coordInfoVecs.num_triangles] = generalParams.maxNodeCount;
		coordInfoVecs.num_triangles += 1;
		//NOTE: What this +1 actually does is that it specifies the location to write
		//any new data. Here it points to the location to write new triangles information.
		//This is a new triangle associated with (tn1, tn2, newnode). Its index is "hostSetInfoVecs.triangles2Nodes_1.size()-4".
		////std::cout<<"GROWTH ERROR 8"<<std::endl;	
		hostSetInfoVecs.triangles2Nodes_1[coordInfoVecs.num_triangles] =(n2);
		hostSetInfoVecs.triangles2Nodes_2[coordInfoVecs.num_triangles] =(n3);
		hostSetInfoVecs.triangles2Nodes_3[coordInfoVecs.num_triangles] = generalParams.maxNodeCount;
		coordInfoVecs.num_triangles += 1;
		//This is a new triangle associated with (tn2, tn3, newnode). Its index is "hostSetInfoVecs.triangles2Nodes_1.size()-3".
		////std::cout<<"GROWTH ERROR 9"<<std::endl;	
		hostSetInfoVecs.triangles2Nodes_1[coordInfoVecs.num_triangles] =(n3);
		hostSetInfoVecs.triangles2Nodes_2[coordInfoVecs.num_triangles] =(n4);
		hostSetInfoVecs.triangles2Nodes_3[coordInfoVecs.num_triangles] = generalParams.maxNodeCount;
		coordInfoVecs.num_triangles += 1;
		//This is a new triangle associated with (tn3, tn1, newnode). Its index is "hostSetInfoVecs.triangles2Nodes_1.size()-2".
		////std::cout<<"GROWTH ERROR 10"<<std::endl;	
		hostSetInfoVecs.triangles2Nodes_1[coordInfoVecs.num_triangles] =(n4);
		hostSetInfoVecs.triangles2Nodes_2[coordInfoVecs.num_triangles] =(n1);
		hostSetInfoVecs.triangles2Nodes_3[coordInfoVecs.num_triangles] = generalParams.maxNodeCount;
		coordInfoVecs.num_triangles += 1;
		//This is a new triangle associated with (tn3, tn1, newnode). Its index is "hostSetInfoVecs.triangles2Nodes_1.size()-1".
		////std::cout<<"GROWTH ERROR 11"<<std::endl;	
		//Now we add new edges formed by the addition of the new node.
		hostSetInfoVecs.edges2Nodes_1[coordInfoVecs.num_edges] = (generalParams.maxNodeCount);
		hostSetInfoVecs.edges2Nodes_2[coordInfoVecs.num_edges] = (n1);
		hostSetInfoVecs.edges2Triangles_1[coordInfoVecs.num_edges] = coordInfoVecs.num_triangles - 4;
		hostSetInfoVecs.edges2Triangles_2[coordInfoVecs.num_edges] = coordInfoVecs.num_triangles - 1;
		coordInfoVecs.num_edges += 1;
		////std::cout<<"GROWTH ERROR 12"<<std::endl;	
		//This is a new edge associated with (newnode, tn1). Its index is "edges2Nodes_1.size()-4".
		hostSetInfoVecs.edges2Nodes_1[coordInfoVecs.num_edges] = (generalParams.maxNodeCount);
		hostSetInfoVecs.edges2Nodes_2[coordInfoVecs.num_edges] = (n2);
		hostSetInfoVecs.edges2Triangles_1[coordInfoVecs.num_edges] = coordInfoVecs.num_triangles - 3;
		hostSetInfoVecs.edges2Triangles_2[coordInfoVecs.num_edges] = coordInfoVecs.num_triangles - 4;
		coordInfoVecs.num_edges += 1;
		////std::cout<<"GROWTH ERROR 13"<<std::endl;	
		//This is a new edge associated with (newnode, tn2). Its index is "edges2Nodes_1.size()-3".
		hostSetInfoVecs.edges2Nodes_1[coordInfoVecs.num_edges] = (generalParams.maxNodeCount);
		hostSetInfoVecs.edges2Nodes_2[coordInfoVecs.num_edges] = (n3);
		hostSetInfoVecs.edges2Triangles_1[coordInfoVecs.num_edges] = coordInfoVecs.num_triangles - 2;
		hostSetInfoVecs.edges2Triangles_2[coordInfoVecs.num_edges] = coordInfoVecs.num_triangles - 3;
		coordInfoVecs.num_edges += 1;
		////std::cout<<"GROWTH ERROR 14"<<std::endl;	
		//This is a new edge associated with (newnode, tn3). Its index is "edges2Nodes_1.size()-2".
		hostSetInfoVecs.edges2Nodes_1[coordInfoVecs.num_edges] = (generalParams.maxNodeCount);
		hostSetInfoVecs.edges2Nodes_2[coordInfoVecs.num_edges] = (n4);
		hostSetInfoVecs.edges2Triangles_1[coordInfoVecs.num_edges] = coordInfoVecs.num_triangles - 1;
		hostSetInfoVecs.edges2Triangles_2[coordInfoVecs.num_edges] = coordInfoVecs.num_triangles - 2;
		coordInfoVecs.num_edges += 1;
		////std::cout<<"GROWTH ERROR 15"<<std::endl;	
		for (int j = 0; j < 4; j++){
		//	//std::cout<<"GROWTH ERROR 16"<<std::endl;				
			//Now we check to see if the order of update is correct, i.e. are edges2Triangles data in correct orientation.
			//This is crucial in the bendingspring computation.
			edgeindex = (coordInfoVecs.num_edges - (4-j));
			a = hostSetInfoVecs.edges2Triangles_1[edgeindex];
			if ((hostSetInfoVecs.triangles2Nodes_1[a] == hostSetInfoVecs.edges2Nodes_1[edgeindex]) && (hostSetInfoVecs.triangles2Nodes_2[a] == hostSetInfoVecs.edges2Nodes_2[edgeindex])){
				a1 = 1;
			}
			else{
				a1 = 0;
			}
			if ((hostSetInfoVecs.triangles2Nodes_2[a] == hostSetInfoVecs.edges2Nodes_1[edgeindex]) && (hostSetInfoVecs.triangles2Nodes_3[a] == hostSetInfoVecs.edges2Nodes_2[edgeindex])){
				a2 = 1;
			}
			else{
				a2 = 0;
			}
			if ((hostSetInfoVecs.triangles2Nodes_3[a] == hostSetInfoVecs.edges2Nodes_1[edgeindex]) && (hostSetInfoVecs.triangles2Nodes_1[a] == hostSetInfoVecs.edges2Nodes_2[edgeindex])){
				a3 = 1;
			}
			else{
				a3 = 0;
			}

			if ((a1+a2+a3) == 0){
				temp1 = hostSetInfoVecs.edges2Triangles_1[edgeindex];
				temp2 = hostSetInfoVecs.edges2Triangles_2[edgeindex];
				hostSetInfoVecs.edges2Triangles_1[edgeindex] = temp2;
				hostSetInfoVecs.edges2Triangles_2[edgeindex] = temp1;
			}
			else{}
			//This checks if the orientation is correct or not, if not, flip the ordering.
		}
		//This is a new edge associated with (newnode, tn3). Its index is "edges2Nodes_1.size()-1".
		generalParams.maxNodeCount += 1;

		hostSetInfoVecs.nndata1[generalParams.maxNodeCount-1] =  (n1);
		hostSetInfoVecs.nndata2[generalParams.maxNodeCount-1] =  (n2);
		hostSetInfoVecs.nndata3[generalParams.maxNodeCount-1] =  (n3);
		hostSetInfoVecs.nndata4[generalParams.maxNodeCount-1] =  (n4);
		hostSetInfoVecs.nndata5[generalParams.maxNodeCount-1] =  (-2);
		hostSetInfoVecs.nndata6[generalParams.maxNodeCount-1] =  (-2);
		hostSetInfoVecs.nndata7[generalParams.maxNodeCount-1] =  (-2);
		hostSetInfoVecs.nndata8[generalParams.maxNodeCount-1] =  (-2);
		hostSetInfoVecs.nndata9[generalParams.maxNodeCount-1] =  (-2);
		//hostSetInfoVecs.nndata10[generalParams.maxNodeCount-1] = (-2);
		//hostSetInfoVecs.nndata11[generalParams.maxNodeCount-1] = (-2);
		//hostSetInfoVecs.nndata12[generalParams.maxNodeCount-1] = (-2);
        hostSetInfoVecs.nodes2Triangles_1[generalParams.maxNodeCount-1] =  (coordInfoVecs.num_triangles - 4);
        hostSetInfoVecs.nodes2Triangles_2[generalParams.maxNodeCount-1] =  (coordInfoVecs.num_triangles - 3);
        hostSetInfoVecs.nodes2Triangles_3[generalParams.maxNodeCount-1] =  (coordInfoVecs.num_triangles - 2);
        hostSetInfoVecs.nodes2Triangles_4[generalParams.maxNodeCount-1] =  (coordInfoVecs.num_triangles - 1);
        hostSetInfoVecs.nodes2Triangles_5[generalParams.maxNodeCount-1] =  (-INT_MAX);
        hostSetInfoVecs.nodes2Triangles_6[generalParams.maxNodeCount-1] =  (-INT_MAX);
        hostSetInfoVecs.nodes2Triangles_7[generalParams.maxNodeCount-1] =  (-INT_MAX);
        hostSetInfoVecs.nodes2Triangles_8[generalParams.maxNodeCount-1] =  (-INT_MAX);
        hostSetInfoVecs.nodes2Triangles_9[generalParams.maxNodeCount-1] =  (-INT_MAX);
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
			if (hostSetInfoVecs.nndata1[nn] == nnn){
				hostSetInfoVecs.nndata1[nn] = nnnn;
			}
			else if (hostSetInfoVecs.nndata2[nn] == nnn){
				hostSetInfoVecs.nndata2[nn] = nnnn;
			}
			else if (hostSetInfoVecs.nndata3[nn] == nnn){
				hostSetInfoVecs.nndata3[nn] = nnnn;
			}
			else if (hostSetInfoVecs.nndata4[nn] == nnn){
				hostSetInfoVecs.nndata4[nn] = nnnn;
			}
			else if (hostSetInfoVecs.nndata5[nn] == nnn){
				hostSetInfoVecs.nndata5[nn] = nnnn;
			}
			else if (hostSetInfoVecs.nndata6[nn] == nnn){
				hostSetInfoVecs.nndata6[nn] = nnnn;
			}
			else if (hostSetInfoVecs.nndata7[nn] == nnn){
				hostSetInfoVecs.nndata7[nn] = nnnn;
			}
			else if (hostSetInfoVecs.nndata8[nn] == nnn){
				hostSetInfoVecs.nndata8[nn] = nnnn;
			}
			else if (hostSetInfoVecs.nndata9[nn] == nnn){
				hostSetInfoVecs.nndata9[nn] = nnnn;
			}
			/*else if (hostSetInfoVecs.nndata10[nn] == nnn){
				hostSetInfoVecs.nndata10[nn] = nnnn;
			}
			 else if (hostSetInfoVecs.nndata11[nn] == nnn){
				hostSetInfoVecs.nndata11[nn] = nnnn;
			}
			else if (hostSetInfoVecs.nndata12[nn] == nnn){
				hostSetInfoVecs.nndata12[nn] = nnnn;
			} */
		}

        for (int j = 0; j < 2; j++){
			int nn, nnn, nnnn;
			if (j == 0){
				nn = n1;
			}
			else if (j == 1){
				nn = n3;
			}

			nnn = elem1;
			nnnn = elem2;

			if (hostSetInfoVecs.nodes2Triangles_1[nn] == nnn || hostSetInfoVecs.nodes2Triangles_1[nn] == nnnn){
				hostSetInfoVecs.nodes2Triangles_1[nn] = -INT_MAX;
			}
			if (hostSetInfoVecs.nodes2Triangles_2[nn] == nnn || hostSetInfoVecs.nodes2Triangles_2[nn] == nnnn){
				hostSetInfoVecs.nodes2Triangles_2[nn] = -INT_MAX;
			}
			if (hostSetInfoVecs.nodes2Triangles_3[nn] == nnn || hostSetInfoVecs.nodes2Triangles_3[nn] == nnnn){
				hostSetInfoVecs.nodes2Triangles_3[nn] = -INT_MAX;
			}
			if (hostSetInfoVecs.nodes2Triangles_4[nn] == nnn || hostSetInfoVecs.nodes2Triangles_4[nn] == nnnn){
				hostSetInfoVecs.nodes2Triangles_4[nn] = -INT_MAX;
			}
			if (hostSetInfoVecs.nodes2Triangles_5[nn] == nnn || hostSetInfoVecs.nodes2Triangles_5[nn] == nnnn){
				hostSetInfoVecs.nodes2Triangles_5[nn] = -INT_MAX;
			}
			if (hostSetInfoVecs.nodes2Triangles_6[nn] == nnn || hostSetInfoVecs.nodes2Triangles_6[nn] == nnnn){
				hostSetInfoVecs.nodes2Triangles_6[nn] = -INT_MAX;
			}
			if (hostSetInfoVecs.nodes2Triangles_7[nn] == nnn || hostSetInfoVecs.nodes2Triangles_7[nn] == nnnn){
				hostSetInfoVecs.nodes2Triangles_7[nn] = -INT_MAX;
			}
			if (hostSetInfoVecs.nodes2Triangles_8[nn] == nnn || hostSetInfoVecs.nodes2Triangles_8[nn] == nnnn){
				hostSetInfoVecs.nodes2Triangles_8[nn] = -INT_MAX;
			}
			if (hostSetInfoVecs.nodes2Triangles_9[nn] == nnn || hostSetInfoVecs.nodes2Triangles_9[nn] == nnnn){
				hostSetInfoVecs.nodes2Triangles_9[nn] = -INT_MAX;
			}

			for (int k = 0; k < 2; k++){
				if (j == 0 && k == 0){
					nnn = coordInfoVecs.num_triangles - 4;
				}
				else if (j == 0 && k == 1){
					nnn = coordInfoVecs.num_triangles - 1;
				}
				else if (j == 1 && k == 0){
					nnn = coordInfoVecs.num_triangles - 3;
				}
				else if (j == 1 && k == 1){
					nnn = coordInfoVecs.num_triangles - 2;
				}
				if (hostSetInfoVecs.nodes2Triangles_1[nn] < 0 ){
					hostSetInfoVecs.nodes2Triangles_1[nn] = nnn;
				}
				else if (hostSetInfoVecs.nodes2Triangles_2[nn] < 0 ){
					hostSetInfoVecs.nodes2Triangles_2[nn] = nnn;
				}
				else if (hostSetInfoVecs.nodes2Triangles_3[nn] < 0 ){
					hostSetInfoVecs.nodes2Triangles_3[nn] = nnn;
				}
				else if (hostSetInfoVecs.nodes2Triangles_4[nn] < 0 ){
					hostSetInfoVecs.nodes2Triangles_4[nn] = nnn;
				}
				else if (hostSetInfoVecs.nodes2Triangles_5[nn] < 0 ){
					hostSetInfoVecs.nodes2Triangles_5[nn] = nnn;
				}
				else if (hostSetInfoVecs.nodes2Triangles_6[nn] < 0 ){
					hostSetInfoVecs.nodes2Triangles_6[nn] = nnn;
				}
				else if (hostSetInfoVecs.nodes2Triangles_7[nn] < 0 ){
					hostSetInfoVecs.nodes2Triangles_7[nn] = nnn;
				}
				else if (hostSetInfoVecs.nodes2Triangles_8[nn] < 0 ){
					hostSetInfoVecs.nodes2Triangles_8[nn] = nnn;
				}
				else if (hostSetInfoVecs.nodes2Triangles_9[nn] < 0 ){
					hostSetInfoVecs.nodes2Triangles_9[nn] = nnn;
				}
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
			if (hostSetInfoVecs.nndata1[nn] < 0){
				hostSetInfoVecs.nndata1[nn] = nnn;
			}
			else if (hostSetInfoVecs.nndata2[nn] < 0){
				hostSetInfoVecs.nndata2[nn] = nnn;
			}
			else if (hostSetInfoVecs.nndata3[nn] < 0){
				hostSetInfoVecs.nndata3[nn] = nnn;
			}
			else if (hostSetInfoVecs.nndata4[nn] < 0){
				hostSetInfoVecs.nndata4[nn] = nnn;
			}
			else if (hostSetInfoVecs.nndata5[nn] < 0){
				hostSetInfoVecs.nndata5[nn] = nnn;
			}
			else if (hostSetInfoVecs.nndata6[nn] < 0){
				hostSetInfoVecs.nndata6[nn] = nnn;
			}
			else if (hostSetInfoVecs.nndata7[nn] < 0){
				hostSetInfoVecs.nndata7[nn] = nnn;
			}
			else if (hostSetInfoVecs.nndata8[nn] < 0){
				hostSetInfoVecs.nndata8[nn] = nnn;
			}
			else if (hostSetInfoVecs.nndata9[nn] < 0){
				hostSetInfoVecs.nndata9[nn] = nnn;
			}
			/*else if (hostSetInfoVecs.nndata10[nn] < 0){
				hostSetInfoVecs.nndata10[nn] = nnn;
			}
			 else if (hostSetInfoVecs.nndata11[nn] < 0){
				hostSetInfoVecs.nndata11[nn] = nnn;
			}
			else if (hostSetInfoVecs.nndata12[nn] < 0){
				hostSetInfoVecs.nndata12[nn] = nnn;
			} */
		}

        for (int j = 0; j < 2; j++){
			int nn, nnn, nnnn;
			if (j == 0){
				nn = n2;
				nnn = elem1;
                //nnnn = elem2;
			}
			else if (j == 1){
				nn = n4;
				nnn = elem2;
                //nnnn = elem1;
			}
			if (hostSetInfoVecs.nodes2Triangles_1[nn]  == nnn){// || hostSetInfoVecs.nodes2Triangles_1[nn]  == nnnn){
				hostSetInfoVecs.nodes2Triangles_1[nn] = -INT_MAX;
			}
			else if (hostSetInfoVecs.nodes2Triangles_2[nn]  == nnn){// || hostSetInfoVecs.nodes2Triangles_2[nn]  == nnnn){
				hostSetInfoVecs.nodes2Triangles_2[nn] = -INT_MAX;
			}
			else if (hostSetInfoVecs.nodes2Triangles_3[nn]  == nnn){// || hostSetInfoVecs.nodes2Triangles_3[nn]  == nnnn){
				hostSetInfoVecs.nodes2Triangles_3[nn] = -INT_MAX;
			}
			else if (hostSetInfoVecs.nodes2Triangles_4[nn]  == nnn){// || hostSetInfoVecs.nodes2Triangles_4[nn]  == nnnn){
				hostSetInfoVecs.nodes2Triangles_4[nn] = -INT_MAX;
			}
			else if (hostSetInfoVecs.nodes2Triangles_5[nn]  == nnn){// || hostSetInfoVecs.nodes2Triangles_5[nn]  == nnnn){
				hostSetInfoVecs.nodes2Triangles_5[nn] = -INT_MAX;
			}
			else if (hostSetInfoVecs.nodes2Triangles_6[nn]  == nnn){// || hostSetInfoVecs.nodes2Triangles_6[nn]  == nnnn){
				hostSetInfoVecs.nodes2Triangles_6[nn] = -INT_MAX;
			}
			else if (hostSetInfoVecs.nodes2Triangles_7[nn]  == nnn){// || hostSetInfoVecs.nodes2Triangles_7[nn]  == nnnn){
				hostSetInfoVecs.nodes2Triangles_7[nn] = -INT_MAX;
			}
			else if (hostSetInfoVecs.nodes2Triangles_8[nn]  == nnn){// || hostSetInfoVecs.nodes2Triangles_8[nn]  == nnnn){
				hostSetInfoVecs.nodes2Triangles_8[nn] = -INT_MAX;
			}
			else if (hostSetInfoVecs.nodes2Triangles_9[nn]  == nnn){// || hostSetInfoVecs.nodes2Triangles_9[nn]  == nnnn){
				hostSetInfoVecs.nodes2Triangles_9[nn] = -INT_MAX;
			}
		}
		
		for (int j = 0; j < 2; j++){
			int nn, nnn, nnnn;
			if (j == 0){
				nn = n2;
			}
			else if (j == 1){
				nn = n4;
			}
			for (int k = 0; k < 2; k++){
				if (j == 0 && k == 0){
					nnn = coordInfoVecs.num_triangles - 4;
				}
				else if (j == 0 && k == 1){
					nnn = coordInfoVecs.num_triangles - 3;
				}
				else if (j == 1 && k == 0){
					nnn = coordInfoVecs.num_triangles - 2;
				}
				else if (j == 1 && k == 1){
					nnn = coordInfoVecs.num_triangles - 1;
				}
				if (hostSetInfoVecs.nodes2Triangles_1[nn] < 0 ){
					hostSetInfoVecs.nodes2Triangles_1[nn] = nnn;
				}
				else if (hostSetInfoVecs.nodes2Triangles_2[nn] < 0 ){
					hostSetInfoVecs.nodes2Triangles_2[nn] = nnn;
				}
				else if (hostSetInfoVecs.nodes2Triangles_3[nn] < 0 ){
					hostSetInfoVecs.nodes2Triangles_3[nn] = nnn;
				}
				else if (hostSetInfoVecs.nodes2Triangles_4[nn] < 0 ){
					hostSetInfoVecs.nodes2Triangles_4[nn] = nnn;
				}
				else if (hostSetInfoVecs.nodes2Triangles_5[nn] < 0 ){
					hostSetInfoVecs.nodes2Triangles_5[nn] = nnn;
				}
				else if (hostSetInfoVecs.nodes2Triangles_6[nn] < 0 ){
					hostSetInfoVecs.nodes2Triangles_6[nn] = nnn;
				}
				else if (hostSetInfoVecs.nodes2Triangles_7[nn] < 0 ){
					hostSetInfoVecs.nodes2Triangles_7[nn] = nnn;
				}
				else if (hostSetInfoVecs.nodes2Triangles_8[nn] < 0 ){
					hostSetInfoVecs.nodes2Triangles_8[nn] = nnn;
				}
				else if (hostSetInfoVecs.nodes2Triangles_9[nn] < 0 ){
					hostSetInfoVecs.nodes2Triangles_9[nn] = nnn;
				}
			}
		}
		//generalParams.num_of_nodes += 1;

		


		////std::cout<<"GROWTH ERROR 17"<<std::endl;	
		//Now we update the edges2Triangles data structure with new edges.
		//std::cout<<"elem 1 = "<<elem1<<std::endl;
		//std::cout<<"elem 2 = "<<elem2<<std::endl;
		for (int i = 0; i < coordInfoVecs.num_edges; i++){
		//	std::cout<<"edges2triangles"<<" "<< i <<" : "<<hostSetInfoVecs.edges2Triangles_1[i]<<" "<<hostSetInfoVecs.edges2Triangles_2[i]<<std::endl;
		}
		int TRIANGLESIZE = coordInfoVecs.num_triangles;//hostSetInfoVecs.triangles2Nodes_1.size();
		if (hostSetInfoVecs.edges2Triangles_1[t1e1] == elem1){
			hostSetInfoVecs.edges2Triangles_1[t1e1] = TRIANGLESIZE-4;
		}
		else if (hostSetInfoVecs.edges2Triangles_2[t1e1] == elem1){
			hostSetInfoVecs.edges2Triangles_2[t1e1] = TRIANGLESIZE-4;
		}
		else{}
		////std::cout<<"GROWTH ERROR 18"<<std::endl;	
		if (hostSetInfoVecs.edges2Triangles_1[t1e2] == elem1){
			hostSetInfoVecs.edges2Triangles_1[t1e2] = TRIANGLESIZE-3;
		}
		else if (hostSetInfoVecs.edges2Triangles_2[t1e2] == elem1){
			hostSetInfoVecs.edges2Triangles_2[t1e2] = TRIANGLESIZE-3;
		}
		else{}
		////std::cout<<"GROWTH ERROR 19"<<std::endl;	
		if (hostSetInfoVecs.edges2Triangles_1[t2e1] == elem2){
			hostSetInfoVecs.edges2Triangles_1[t2e1] = TRIANGLESIZE-2;
		}
		else if (hostSetInfoVecs.edges2Triangles_2[t2e1] == elem2){
			hostSetInfoVecs.edges2Triangles_2[t2e1] = TRIANGLESIZE-2;
		}
		else{}
		////std::cout<<"GROWTH ERROR 20"<<std::endl;	
		if (hostSetInfoVecs.edges2Triangles_1[t2e2] == elem2){
			hostSetInfoVecs.edges2Triangles_1[t2e2] = TRIANGLESIZE-1;
		}
		else if (hostSetInfoVecs.edges2Triangles_2[t2e2] == elem2){
			hostSetInfoVecs.edges2Triangles_2[t2e2] = TRIANGLESIZE-1;
		}
		else{}
		//std::cout<<"t1e1 "<<t1e1<<std::endl;
		//std::cout<<"t1e2 "<<t1e2<<std::endl;
		//std::cout<<"t1e3 "<<t1e3<<std::endl;
		//std::cout<<"t2e1 "<<t2e1<<std::endl;
		//std::cout<<"t2e2 "<<t2e2<<std::endl;
		//std::cout<<"t2e3 "<<t2e3<<std::endl;

		//for (int i = 0; i < coordInfoVecs.num_edges; i++){
		//	std::cout<<"edges2triangles"<<" "<< i <<" : "<<hostSetInfoVecs.edges2Triangles_1[i]<<" "<<hostSetInfoVecs.edges2Triangles_2[i]<<std::endl;
		//}
		//The above change the existing edges2Triangles data structure to accomodate new triangles added.

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		
		//Now we will take care of the last unedited data structure "triangles2Edges".
		//int aa, bb;
		int EDGESIZE = coordInfoVecs.num_edges;//hostSetInfoVecs.edges2Nodes_1.size();
		for (int j = 0; j < 4; j++){
		//	//std::cout<<"GROWTH ERROR 21"<<std::endl;	
			if (j == 0){
				hostSetInfoVecs.triangles2Edges_1[coordInfoVecs.num_triangles - 4] = (EDGESIZE-4);
				hostSetInfoVecs.triangles2Edges_2[coordInfoVecs.num_triangles - 4] = (t1e1);
				hostSetInfoVecs.triangles2Edges_3[coordInfoVecs.num_triangles - 4] = (EDGESIZE-3);   
			}
			else if (j == 1){
				hostSetInfoVecs.triangles2Edges_1[coordInfoVecs.num_triangles - 3] = (EDGESIZE-3);
				hostSetInfoVecs.triangles2Edges_2[coordInfoVecs.num_triangles - 3] = (t1e2);
				hostSetInfoVecs.triangles2Edges_3[coordInfoVecs.num_triangles - 3] = (EDGESIZE-2);   
			}
			else if (j ==2){
				hostSetInfoVecs.triangles2Edges_1[coordInfoVecs.num_triangles - 2] = (EDGESIZE-2);
				hostSetInfoVecs.triangles2Edges_2[coordInfoVecs.num_triangles - 2] = (t2e1);
				hostSetInfoVecs.triangles2Edges_3[coordInfoVecs.num_triangles - 2] = (EDGESIZE-1);   
			}
			else if (j ==3){
				hostSetInfoVecs.triangles2Edges_1[coordInfoVecs.num_triangles - 1] = (EDGESIZE-1);
				hostSetInfoVecs.triangles2Edges_2[coordInfoVecs.num_triangles - 1] = (t2e2);
				hostSetInfoVecs.triangles2Edges_3[coordInfoVecs.num_triangles - 1] = (EDGESIZE-4);   
			}
			
		}
	
		
		if (hostSetInfoVecs.nodes_in_upperhem[hostSetInfoVecs.edges2Nodes_1[iedge]] == 1 && generalParams.nodes_in_upperhem[hostSetInfoVecs.edges2Nodes_2[iedge]] == 1){
			//generalParams.nodes_in_upperhem.push_back(1);
            hostSetInfoVecs.nodes_in_upperhem[generalParams.maxNodeCount - 1] = 1;

		}
		else if (hostSetInfoVecs.nodes_in_upperhem[hostSetInfoVecs.edges2Nodes_1[iedge]] == 1 || generalParams.nodes_in_upperhem[hostSetInfoVecs.edges2Nodes_2[iedge]] == 1){
			//generalParams.nodes_in_upperhem.push_back(0);
            hostSetInfoVecs.nodes_in_upperhem[generalParams.maxNodeCount - 1] = 0;
		}
		else{
			//generalParams.nodes_in_upperhem.push_back(-1);
            hostSetInfoVecs.nodes_in_upperhem[generalParams.maxNodeCount - 1] = -1;
		}
		//Finally, we will fill the edge data chosen for growth (expansion) with INT_MAX so its data is no longer relevant to the computation
		////std::cout<<"GROWTH ERROR 22"<<std::endl;	
		hostSetInfoVecs.edges2Nodes_1[iedge] = INT_MAX;
		hostSetInfoVecs.edges2Nodes_2[iedge] = INT_MAX;
		for (int i = 0; i < coordInfoVecs.num_triangles; i++){
		//	//std::cout<<"GROWTH ERROR 23"<<std::endl;	
			if (hostSetInfoVecs.triangles2Edges_1[i] == iedge){
				hostSetInfoVecs.triangles2Edges_1[i] = INT_MAX;
			}
			if (hostSetInfoVecs.triangles2Edges_2[i] == iedge){
				hostSetInfoVecs.triangles2Edges_2[i] = INT_MAX;
			}
			if (hostSetInfoVecs.triangles2Edges_3[i] == iedge){
				hostSetInfoVecs.triangles2Edges_3[i] = INT_MAX;
			}
		}
		hostSetInfoVecs.edges2Triangles_1[iedge] = INT_MAX;
		hostSetInfoVecs.edges2Triangles_2[iedge] = INT_MAX;
		
		////std::cout<<"GROWTH ERROR 24"<<std::endl;	
		
			hostSetInfoVecs.triangles2Nodes_1[elem1] = INT_MAX;
			hostSetInfoVecs.triangles2Nodes_2[elem1] = INT_MAX;
			hostSetInfoVecs.triangles2Nodes_3[elem1] = INT_MAX;
			hostSetInfoVecs.triangles2Nodes_1[elem2] = INT_MAX;
			hostSetInfoVecs.triangles2Nodes_2[elem2] = INT_MAX;
			hostSetInfoVecs.triangles2Nodes_3[elem2] = INT_MAX;
			
			//Delete the associated vertices information of the selected triangle.
			//Since we delete the chosen triangles, any triangle indexed lower than the deleted one will have its index reduced (or moved up) by 1.
			//Hence, we need to sweep through all data structures using the triangle index to change the index accordingly.
			for (int i = 0; i < coordInfoVecs.num_edges; i++){
		//		//std::cout<<"GROWTH ERROR 25"<<std::endl;	
				if (hostSetInfoVecs.edges2Triangles_1[i] == elem1 || hostSetInfoVecs.edges2Triangles_1[i] == elem2){
					hostSetInfoVecs.edges2Triangles_1[i] = INT_MAX;
				}
				if (hostSetInfoVecs.edges2Triangles_2[i] == elem1 || hostSetInfoVecs.edges2Triangles_2[i] == elem2){
					hostSetInfoVecs.edges2Triangles_2[i] = INT_MAX;
				}
			if (hostSetInfoVecs.edges2Triangles_1[i] != INT_MAX && hostSetInfoVecs.edges2Triangles_2[i] == INT_MAX){
				std::cout<<"modified edges2Triangles "<<hostSetInfoVecs.edges2Triangles_1[i]<<" "<<hostSetInfoVecs.edges2Triangles_2[i]<<std::endl;
				}
				else if (hostSetInfoVecs.edges2Triangles_1[i] == INT_MAX && hostSetInfoVecs.edges2Triangles_2[i] != INT_MAX){
					std::cout<<"modified edges2Triangles "<<hostSetInfoVecs.edges2Triangles_1[i]<<" "<<hostSetInfoVecs.edges2Triangles_2[i]<<std::endl;
					}
			}
			//This completes the sweep. After this, the indices of triangle used in edges2Triangles data structure should be the correct one.
		//	//std::cout<<"GROWTH ERROR 26"<<std::endl;	
			hostSetInfoVecs.triangles2Edges_1[elem1] = INT_MAX;
			hostSetInfoVecs.triangles2Edges_2[elem1] = INT_MAX;
			hostSetInfoVecs.triangles2Edges_3[elem1] = INT_MAX;
			hostSetInfoVecs.triangles2Edges_1[elem2] = INT_MAX;
			hostSetInfoVecs.triangles2Edges_2[elem2] = INT_MAX;
			hostSetInfoVecs.triangles2Edges_3[elem2] = INT_MAX;
			for (int i = 0; i < coordInfoVecs.num_triangles; i++){
		//		//std::cout<<"GROWTH ERROR 27"<<std::endl;	
				if (hostSetInfoVecs.triangles2Edges_1[i] == iedge){
					hostSetInfoVecs.triangles2Edges_1[i] = INT_MAX;
				}
				if (hostSetInfoVecs.triangles2Edges_2[i] == iedge ){
					hostSetInfoVecs.triangles2Edges_2[i] = INT_MAX;
				}
				if (hostSetInfoVecs.triangles2Edges_3[i] == iedge ){
					hostSetInfoVecs.triangles2Edges_3[i] = INT_MAX;
				}
			}
		
		//Erase the edge infomation related to the deleted triangle. Note the deletion should always start with the largest index.

		//Before we delete the edge, determine whether the newly added node is part of nodes_in_upperhem or not.
		

		
						//Erase the edge infomation related to the deleted triangle.

						//Now we update the nodes_in_upperhem and edges_in_upperhem data structures.
						//This ensures that the newly created edges will have the correct associated spring constant.
//std::cout<<"ERROR HERE?"<<std::endl;
		generalParams.edges_in_upperhem[iedge] = INT_MAX;
		for (int i = 0; i < coordInfoVecs.num_edges; i++){
			if (generalParams.edges_in_upperhem_list[i] == iedge){
				generalParams.edges_in_upperhem_list[i] = INT_MAX;
				//break;
			}
		}

		for (int q = 0; q < 4; q++){
		//	//std::cout<<"GROWTH ERROR 30"<<std::endl;	
			int nodeP = hostSetInfoVecs.triangles2Nodes_1[coordInfoVecs.num_triangles - (4-q)]; 
			int nodeQ = hostSetInfoVecs.triangles2Nodes_2[coordInfoVecs.num_triangles - (4-q)];
			int nodeR = hostSetInfoVecs.triangles2Nodes_3[coordInfoVecs.num_triangles - (4-q)];

			if (hostSetInfoVecs.nodes_in_upperhem[nodeP]==1 && hostSetInfoVecs.nodes_in_upperhem[nodeQ] ==1 && hostSetInfoVecs.nodes_in_upperhem[nodeR] ==1){
				//generalParams.triangles_in_upperhem.push_back(1);
                hostSetInfoVecs.triangles_in_upperhem[coordInfoVecs.num_triangles - (4-q)] = 1;
			}
			else if (hostSetInfoVecs.nodes_in_upperhem[nodeP]==1 && hostSetInfoVecs.nodes_in_upperhem[nodeQ] ==1){
				//generalParams.triangles_in_upperhem.push_back(0);
                hostSetInfoVecs.triangles_in_upperhem[coordInfoVecs.num_triangles - (4-q)] = 0;
			}
			else if (hostSetInfoVecs.nodes_in_upperhem[nodeP]==1 && hostSetInfoVecs.nodes_in_upperhem[nodeR] ==1){
				//generalParams.triangles_in_upperhem.push_back(0);
                hostSetInfoVecs.triangles_in_upperhem[coordInfoVecs.num_triangles - (4-q)] = 0;
			}
			else if (hostSetInfoVecs.nodes_in_upperhem[nodeQ] ==1 && hostSetInfoVecs.nodes_in_upperhem[nodeR] ==1){
				//generalParams.triangles_in_upperhem.push_back(0);
                hostSetInfoVecs.triangles_in_upperhem[coordInfoVecs.num_triangles - (4-q)] = 0;
			}
			else{
				//generalParams.triangles_in_upperhem.push_back(INT_MAX);
                hostSetInfoVecs.triangles_in_upperhem[coordInfoVecs.num_triangles - (4-q)] = INT_MAX;
			}
		}
		//std::cout<<"edges2Triangles size"<<""<<hostSetInfoVecs.edges2Triangles_1.size()<<" "<<hostSetInfoVecs.edges2Triangles_2.size()<<std::endl;
		//std::cout<<"triangles_in_upperhem size "<<generalParams.triangles_in_upperhem.size()<<std::endl;	
		//std::cout<<"GROWTH ERROR 29"<<std::endl;	
		for (int q = 0; q < 4; q++){
			//std::cout<<"GROWTH ERROR 31"<<std::endl;	
			int elem_1 = hostSetInfoVecs.edges2Triangles_1[coordInfoVecs.num_edges-(4 - q)];
			//std::cout<<coordInfoVecs.num_edges-(4 - q)<<std::endl;
			//std::cout<<"elem_1 "<<elem_1<<std::endl;
			//std::cout<<generalParams.nodes_in_upperhem[nodeP]<<std::endl;
			int elem_2 = hostSetInfoVecs.edges2Triangles_2[coordInfoVecs.num_edges-(4 - q)];
			//std::cout<<"elem_2"<<elem_2<<std::endl;
			//std::cout<<generalParams.nodes_in_upperhem[nodeQ]<<std::endl;
			//std::cout<<"GROWTH ERROR 31.5"<<std::endl;
			if (hostSetInfoVecs.triangles_in_upperhem[elem_1] == 1 && hostSetInfoVecs.triangles_in_upperhem[elem_2] == 1){
				
				//generalParams.edges_in_upperhem.push_back(1);
                hostSetInfoVecs.edges_in_upperhem[coordInfoVecs.num_edges - (4-q)] = 1;
				//generalParams.edges_in_upperhem_index.push_back(generalParams.num_of_edges - (4 - q));
				//generalParams.edges_in_upperhem_list.push_back(coordInfoVecs.num_edges - (4 - q));
                hostSetInfoVecs.edges_in_upperhem_list[coordInfoVecs.num_edges - (4-q)] = coordInfoVecs.num_edges - (4-q);
                generalParams.edges_in_upperhem_list_length += 1;
			}
			
			else if (generalParams.triangles_in_upperhem[elem_1] == 1 || generalParams.triangles_in_upperhem[elem_2] == 1){
				
				//generalParams.edges_in_upperhem.push_back(0);
                hostSetInfoVecs.edges_in_upperhem[coordInfoVecs.num_edges - (4-q)] = 0;
                hostSetInfoVecs.edges_in_upperhem_list[coordInfoVecs.num_edges - (4-q)] = -INT_MAX;
                generalParams.edges_in_upperhem_list_length += 1;
				//generalParams.edges_in_upperhem_index.push_back(generalParams.num_of_edges - (4 - q));
			}
			else{
				
				//generalParams.edges_in_upperhem.push_back(-1);
                hostSetInfoVecs.edges_in_upperhem[coordInfoVecs.num_edges - (4-q)] = -1;
                hostSetInfoVecs.edges_in_upperhem_list[coordInfoVecs.num_edges - (4-q)] = -INT_MAX;
                generalParams.edges_in_upperhem_list_length += 1;
			}

			//generalParams.boundaries_in_upperhem.push_back(-1);
            hostSetInfoVecs.boundaries_in_upperhem[coordInfoVecs.num_edges - (4-q)] = -1;		
			
		}
		generalParams.triangles_in_upperhem[elem1] = INT_MAX;
		generalParams.triangles_in_upperhem[elem2] = INT_MAX;
      }
    }  

    if (GROWTH_ACCEPTED == true){
        double midpoint_x = (hostSetInfoVecs.nodeLocX[edge_start] + hostSetInfoVecs.nodeLocX[edge_end])/2.0;
        // std::cout<<"wwwww"<<std::endl;
        double midpoint_y = (hostSetInfoVecs.nodeLocY[edge_start] + hostSetInfoVecs.nodeLocY[edge_end])/2.0;
        // std::cout<<"wwwwww"<<std::endl;
        double midpoint_z = (hostSetInfoVecs.nodeLocZ[edge_start] + hostSetInfoVecs.nodeLocZ[edge_end])/2.0;
        std::cout<<"midpoint of edge chosen for growth : "<<midpoint_x<<", "<<midpoint_y<<", "<<midpoint_z<<std::endl;
        alpha = iedge;//1;
    }
    else{
        // alpha = -1;
        // generalParams.triangle_undergoing_growth.clear();
        // alpha.push_back(-1);
        alpha = -1;
    }
    // return alpha;
    return alpha;
}

int Utilities::nodes2Triangles_host_vecs(
    int inode,
    HostSetInfoVecs& hostSetInfoVecs,
    CoordInfoVecs& coordInfoVecs,
	GeneralParams& generalParams,
    AuxVecs& auxVecs
){

    hostSetInfoVecs.nodes2Triangles_1[inode] = -INT_MAX;
    hostSetInfoVecs.nodes2Triangles_2[inode] = -INT_MAX;
    hostSetInfoVecs.nodes2Triangles_3[inode] = -INT_MAX;
    hostSetInfoVecs.nodes2Triangles_4[inode] = -INT_MAX;
    hostSetInfoVecs.nodes2Triangles_5[inode] = -INT_MAX;
    hostSetInfoVecs.nodes2Triangles_6[inode] = -INT_MAX;
    hostSetInfoVecs.nodes2Triangles_7[inode] = -INT_MAX;
    hostSetInfoVecs.nodes2Triangles_8[inode] = -INT_MAX;
    hostSetInfoVecs.nodes2Triangles_9[inode] = -INT_MAX;
    //int coordInfoVecs.nodes2Triangles_10 = -INT_MAX;
    //int coordInfoVecs.nodes2Triangles_11 = -INT_MAX;
    //int coordInfoVecs.nodes2Triangles_12 = -INT_MAX;

    //for now iterate through all membrane id's
    // int begin = keyBegin[bucketId];
    // int end = keyEnd[bucketId];

    //First, create a new vector detailing all the neighboring nodes of node i//
    
        int COUNT = 0;
        for (int k = 0; k < coordInfoVecs.num_triangles; k++){
            int v1 = hostSetInfoVecs.triangles2Nodes_1[k];
            int v2 = hostSetInfoVecs.triangles2Nodes_2[k];
            int v3 = hostSetInfoVecs.triangles2Nodes_3[k];
            if (v1 < 0 || v2 < 0 || v3 < 0){
                continue;
            }
            else if (v1 == INT_MAX || v2 == INT_MAX || v3 == INT_MAX){
                continue;
            }

            if (v1 == inode || v2 == inode || v3 == inode){
                
                COUNT += 1;
                if (COUNT == 1){
                    hostSetInfoVecs.nodes2Triangles_1[inode] = k;
                }
                else if (COUNT == 2){
                    hostSetInfoVecs.nodes2Triangles_2[inode] = k;
                }
                else if (COUNT == 3){
                    hostSetInfoVecs.nodes2Triangles_3[inode] = k;
                }
                else if (COUNT == 4){
                    hostSetInfoVecs.nodes2Triangles_4[inode] = k;
                }
                else if (COUNT == 5){
                    hostSetInfoVecs.nodes2Triangles_5[inode] = k;
                }
                else if (COUNT == 6){
                    hostSetInfoVecs.nodes2Triangles_6[inode] = k;
                }
                else if (COUNT == 7){
                    hostSetInfoVecs.nodes2Triangles_7[inode] = k;
                }
                else if (COUNT == 8){
                    hostSetInfoVecs.nodes2Triangles_8[inode] = k;
                }
                else if (COUNT == 9){
                    hostSetInfoVecs.nodes2Triangles_9[inode] = k;
                }
                /*else if (COUNT == 10){
                    coordInfoVecs.nodes2Triangles_10 = k;
                }*/
                /* else if (COUNT == 11){
                    coordInfoVecs.nodes2Triangles_11 = k;
                }
                else if (COUNT == 12){
                    coordInfoVecs.nodes2Triangles_12 = k;
                } */

            }
        }

}
void Utilities::triangles2Triangles_host_vecs(
    int elem,
    HostSetInfoVecs& hostSetInfoVecs,
    CoordInfoVecs& coordInfoVecs,
	GeneralParams& generalParams,
    AuxVecs& auxVecs){
        if (hostSetInfoVecs.triangles2Edges_1[elem] >= (INT_MAX-1000) ||
            hostSetInfoVecs.triangles2Edges_1[elem] <= (-INT_MAX+1000) ||
            hostSetInfoVecs.triangles2Edges_1[elem] < 0 ||
            hostSetInfoVecs.triangles2Edges_2[elem] >= (INT_MAX-1000)  ||
            hostSetInfoVecs.triangles2Edges_2[elem] <= (-INT_MAX+1000) ||
            hostSetInfoVecs.triangles2Edges_2[elem] < 0 ||
            hostSetInfoVecs.triangles2Edges_1[elem] >= (INT_MAX-1000) ||
            hostSetInfoVecs.triangles2Edges_1[elem] <= (-INT_MAX+1000)||
            hostSetInfoVecs.triangles2Edges_1[elem] < 0
        ){
            hostSetInfoVecs.triangles2Triangles_1[elem] = -1.0;//-INT_MAX;
            hostSetInfoVecs.triangles2Triangles_2[elem] = -1.0;//-INT_MAX;
            hostSetInfoVecs.triangles2Triangles_3[elem] = -1.0;//-INT_MAX;
        }
        else{

            int edge1_n1 = hostSetInfoVecs.triangles2Nodes_1[elem];
            int edge1_n2 = hostSetInfoVecs.triangles2Nodes_2[elem];
            int edge2_n1 = hostSetInfoVecs.triangles2Nodes_2[elem];
            int edge2_n2 = hostSetInfoVecs.triangles2Nodes_3[elem];
            int edge3_n1 = hostSetInfoVecs.triangles2Nodes_3[elem];
            int edge3_n2 = hostSetInfoVecs.triangles2Nodes_1[elem];
            if (edge1_n1 != edge3_n2){
                std::cout<<"error: triangles2Triangles_host_vecs() => SOMETHING WENT WRONG SETTING UP THE NODE ORDER OF EDGES!"<<std::endl;
            }
            int edge;
            for (int i = 0; i < 3; i++){
                //std::cout<<"i = "<<i<<std::endl;
                if (i == 0){
                    edge = hostSetInfoVecs.triangles2Edges_1[elem];
                }
                else if (i == 1){
                    edge = hostSetInfoVecs.triangles2Edges_2[elem];
                }
                else if (i == 2){
                    edge = hostSetInfoVecs.triangles2Edges_3[elem];
                }
                
                if ((hostSetInfoVecs.edges2Nodes_1[edge] == edge1_n1 && hostSetInfoVecs.edges2Nodes_2[edge] == edge1_n2) ||
                    (hostSetInfoVecs.edges2Nodes_1[edge] == edge1_n2 && hostSetInfoVecs.edges2Nodes_2[edge] == edge1_n1)){
                        if (hostSetInfoVecs.edges2Triangles_1[edge] == elem){
                            hostSetInfoVecs.triangles2Triangles_1[elem] = hostSetInfoVecs.edges2Triangles_2[edge];
                        }
                        else if (hostSetInfoVecs.edges2Triangles_2[edge] == elem){
                            hostSetInfoVecs.triangles2Triangles_1[elem] = hostSetInfoVecs.edges2Triangles_1[edge];
                        }
                        else{
                            std::cout<<"error: triangles2Triangles_host_vecs() => SOMETHING WENT WRONG SETTING CORRESPONDING ELEM TO EACH EDGE!"<<std::endl;
                        }
                    }
                else if ((hostSetInfoVecs.edges2Nodes_1[edge] == edge2_n1 && hostSetInfoVecs.edges2Nodes_2[edge] == edge2_n2) ||
                    (hostSetInfoVecs.edges2Nodes_1[edge] == edge2_n2 && hostSetInfoVecs.edges2Nodes_2[edge] == edge2_n1)){
                        if (hostSetInfoVecs.edges2Triangles_1[edge] == elem){
                            hostSetInfoVecs.triangles2Triangles_2[elem] = hostSetInfoVecs.edges2Triangles_2[edge];
                        }
                        else if (hostSetInfoVecs.edges2Triangles_2[edge] == elem){
                            hostSetInfoVecs.triangles2Triangles_2[elem] = hostSetInfoVecs.edges2Triangles_1[edge];
                        }
                        else{
                            std::cout<<"error: triangles2Triangles_host_vecs() => SOMETHING WENT WRONG SETTING CORRESPONDING ELEM TO EACH EDGE!"<<std::endl;
                        }
                    }
                else if ((hostSetInfoVecs.edges2Nodes_1[edge] == edge3_n1 && hostSetInfoVecs.edges2Nodes_2[edge] == edge3_n2) ||
                    (hostSetInfoVecs.edges2Nodes_1[edge] == edge3_n2 && hostSetInfoVecs.edges2Nodes_2[edge] == edge3_n1)){
                        if (hostSetInfoVecs.edges2Triangles_1[edge] == elem){
                            hostSetInfoVecs.triangles2Triangles_3[elem] = hostSetInfoVecs.edges2Triangles_2[edge];
                        }
                        else if (hostSetInfoVecs.edges2Triangles_2[edge] == elem){
                            hostSetInfoVecs.triangles2Triangles_3[elem] = hostSetInfoVecs.edges2Triangles_1[edge];
                        }
                        else{
                            std::cout<<"error: triangles2Triangles_host_vecs() => SOMETHING WENT WRONG SETTING CORRESPONDING ELEM TO EACH EDGE!"<<std::endl;
                        }
                    }
                else{
                    std::cout<<"error: triangles2Triangles_host_vecs() => edges2Nodes info did not agree with any nodes on the triangles listing"<<std::endl;
                }
            }
            if (hostSetInfoVecs.triangles2Triangles_1[elem] == hostSetInfoVecs.triangles2Triangles_2[elem]){
                std::cout<<"error: triangles2Triangles_host_vecs() => overlapping triangles2Triangles info for a given elem ->"<<elem<<std::endl;
            }
            else if (hostSetInfoVecs.triangles2Triangles_1[elem] == hostSetInfoVecs.triangles2Triangles_3[elem]){
                std::cout<<"error: triangles2Triangles_host_vecs() => overlapping triangles2Triangles info for a given elem ->"<<elem<<std::endl;
            }
            else if (hostSetInfoVecs.triangles2Triangles_2[elem] == hostSetInfoVecs.triangles2Triangles_3[elem]){
                std::cout<<"error: triangles2Triangles_host_vecs() => overlapping triangles2Triangles info for a given elem ->"<<elem<<std::endl;
            }
        }
       // std::cout<<"triangles2Triangles["<<elem<<"] -> "<<hostSetInfoVecs.triangles2Triangles_1[elem]<<" "<<hostSetInfoVecs.triangles2Triangles_2[elem]<<" "<<hostSetInfoVecs.triangles2Triangles_3[elem]<<std::endl;
    }

// void Utilities::triangles2Triangles_host_vecs(
//     int elem,
//     HostSetInfoVecs& hostSetInfoVecs,
//     CoordInfoVecs& coordInfoVecs,
// 	GeneralParams& generalParams,
//     AuxVecs& auxVecs){
//         if (hostSetInfoVecs.triangles2Edges_1[elem] == INT_MAX ||
//             hostSetInfoVecs.triangles2Edges_1[elem] == -INT_MAX ||
//             hostSetInfoVecs.triangles2Edges_1[elem] < 0 ||
//             hostSetInfoVecs.triangles2Edges_2[elem] == INT_MAX ||
//             hostSetInfoVecs.triangles2Edges_2[elem] == -INT_MAX ||
//             hostSetInfoVecs.triangles2Edges_2[elem] < 0 ||
//             hostSetInfoVecs.triangles2Edges_1[elem] == INT_MAX ||
//             hostSetInfoVecs.triangles2Edges_1[elem] == -INT_MAX ||
//             hostSetInfoVecs.triangles2Edges_1[elem] < 0
//         ){
//             hostSetInfoVecs.triangles2Triangles_1[elem] = -1.0;//-INT_MAX;
//             hostSetInfoVecs.triangles2Triangles_2[elem] = -1.0;//-INT_MAX;
//             hostSetInfoVecs.triangles2Triangles_3[elem] = -1.0;//-INT_MAX;
//         }
//         else{

//             int edge1 = hostSetInfoVecs.triangles2Edges_1[elem];
//             int edge2 = hostSetInfoVecs.triangles2Edges_2[elem];
//             int edge3 = hostSetInfoVecs.triangles2Edges_3[elem];
//             int edge;
//             for (int i = 0; i < 3; i++){
//                 //std::cout<<"i = "<<i<<std::endl;
//                 if (i == 0){
//                     edge = edge1;
//                 }
//                 else if (i == 1){
//                     edge = edge2;
//                 }
//                 else if (i == 2){
//                     edge = edge3;
//                 }
//                 int triangle1 = hostSetInfoVecs.edges2Triangles_1[edge];
//                 int triangle2 = hostSetInfoVecs.edges2Triangles_2[edge];
//                 if (i == 0 && triangle1 == elem){
//                     // std::cout<<"elem = "<<elem<<std::endl;
//                     // std::cout<<"triangle1 = "<<triangle1<<std::endl;
//                     hostSetInfoVecs.triangles2Triangles_1[elem] = triangle2;
//                 }
//                 else if (i == 0 && triangle2 == elem){
//                     // std::cout<<"elem = "<<elem<<std::endl;
//                     // std::cout<<"triangle2 = "<<triangle2<<std::endl;
//                     hostSetInfoVecs.triangles2Triangles_1[elem] = triangle1;
//                 }
//                 else if (i == 0 && triangle1 != elem && triangle2 != elem){
//                     // std::cout<<"elem = "<<elem<<std::endl;
//                     // std::cout<<"triangle1 = "<<triangle1<<std::endl;
//                     // std::cout<<"triangle2 = "<<triangle1<<std::endl;
//                     std::cout<<"SOMETHING WENT WRONG CREATING triangles2Triangles DATA STRUCTURE"<<std::endl;
//                 }

//                 if (i == 1 && triangle1 == elem){
//                     // std::cout<<"elem = "<<elem<<std::endl;
//                     // std::cout<<"triangle1 = "<<triangle1<<std::endl;
//                     //std::cout<<"triangle2 = "<<triangle1<<std::endl;
//                     hostSetInfoVecs.triangles2Triangles_2[elem] = triangle2;
//                 }
//                 else if (i == 1 && triangle2 == elem){
//                     // std::cout<<"elem = "<<elem<<std::endl;
//                     // //std::cout<<"triangle1 = "<<triangle1<<std::endl;
//                     // std::cout<<"triangle2 = "<<triangle1<<std::endl;
//                     hostSetInfoVecs.triangles2Triangles_2[elem] = triangle1;
//                 }
//                 else if (i == 1 && triangle1 != elem && triangle2 != elem){
//                     // std::cout<<"elem = "<<elem<<std::endl;
//                     // std::cout<<"triangle1 = "<<triangle1<<std::endl;
//                     // std::cout<<"triangle2 = "<<triangle1<<std::endl;
//                     std::cout<<"SOMETHING WENT WRONG CREATING triangles2Triangles DATA STRUCTURE"<<std::endl;
//                 }

//                 if (i == 2 && triangle1 == elem){
//                     // std::cout<<"elem = "<<elem<<std::endl;
//                     // std::cout<<"triangle1 = "<<triangle1<<std::endl;
//                     //std::cout<<"triangle2 = "<<triangle1<<std::endl;
//                     hostSetInfoVecs.triangles2Triangles_3[elem] = triangle2;
//                 }
//                 else if (i == 2 && triangle2 == elem){
//                     // std::cout<<"elem = "<<elem<<std::endl;
//                     // //std::cout<<"triangle1 = "<<triangle1<<std::endl;
//                     // std::cout<<"triangle2 = "<<triangle1<<std::endl;
//                     hostSetInfoVecs.triangles2Triangles_3[elem] = triangle1;
//                 }
//                 else if (i == 2 && triangle1 != elem && triangle2 != elem){
//                     // std::cout<<"elem = "<<elem<<std::endl;
//                     // std::cout<<"triangle1 = "<<triangle1<<std::endl;
//                     // std::cout<<"triangle2 = "<<triangle1<<std::endl;
//                     std::cout<<"SOMETHING WENT WRONG CREATING triangles2Triangles DATA STRUCTURE"<<std::endl;
//                 }
//             }
//         }
//     }

int Utilities::edge_swap_host_vecs(
    int iedge, 
    GeneralParams& generalParams,
    HostSetInfoVecs& hostSetInfoVecs,
    LinearSpringInfoVecs& linearSpringInfoVecs,
    BendingTriangleInfoVecs& bendingTriangleInfoVecs,
    AreaTriangleInfoVecs& areaTriangleInfoVecs) {

    // double   generalParams.scaling_pow = 4.0;
    int alpha = 0;
        
    int HEAD,TAIL;
    int H0, T0,H1,H2,T1,T2;
    int edge_start, edge_end;
    int a1, b1, c1, a2, b2, c2;
    double temp_bend = 0.0;
    double linear_spring_constant;
    double bend_spring_constant;
    double preferred_angle;
    double vol_0, vol_1;
    double P0x_vol1, P0y_vol1, P0z_vol1, P0x_vol2, P0y_vol2, P0z_vol2;
    double N1x_vol, N1y_vol, N1z_vol, N2x_vol, N2y_vol, N2z_vol;
    // std::cout<<"iedge = "<<iedge<<std::endl;
    //std::cout<<hostSetInfoVecs.edges2Triangles_1[iedge] <<" "<< hostSetInfoVecs.edges2Triangles_2[iedge]<<std::endl;   
    if ( hostSetInfoVecs.edges2Triangles_1[iedge] != hostSetInfoVecs.edges2Triangles_2[iedge]){
        H0 = hostSetInfoVecs.edges2Triangles_1[iedge];//index of the 1st triangle to i-th edge
        //std::cout<<"H0 = "<<H0<<std::endl;
        T0 = hostSetInfoVecs.edges2Triangles_2[iedge];//index of the 2nd triangle to i-th edge
        //std::cout<<"T0 = "<<T0<<std::endl;
        edge_start = hostSetInfoVecs.edges2Nodes_1[iedge];//index of the 1st node of i-th edge
        edge_end = hostSetInfoVecs.edges2Nodes_2[iedge];//index of the 2nd node of i-th edge

        a1 = hostSetInfoVecs.triangles2Edges_1[H0];//index of the 1st node of triangle H0
        //std::cout<<"a1 = "<<a1<<std::endl;
        b1 = hostSetInfoVecs.triangles2Edges_2[H0];//index of the 2nd node of triangle H0
        //std::cout<<"b1 = "<<b1<<std::endl;
        c1 = hostSetInfoVecs.triangles2Edges_3[H0];//index of the 3rd node of triangle H0
        //std::cout<<"c1 = "<<c1<<std::endl;
        
        a2 = hostSetInfoVecs.triangles2Edges_1[T0];//index of the 1st node of triangle T0
        //std::cout<<"a2 = "<<a2<<std::endl;
        b2 = hostSetInfoVecs.triangles2Edges_2[T0];
        //std::cout<<"b2 = "<<b2<<std::endl;
        c2 = hostSetInfoVecs.triangles2Edges_3[T0];
        //std::cout<<"c2 = "<<c2<<std::endl;
        
        //Now we identify the edge indices associated with the small subsystem.
        //This gives us the indices for H1, H2, T1, T2 (see the figure below).
        if (a1 != iedge && hostSetInfoVecs.edges2Nodes_1[a1] == edge_start){H1 = a1;}
        else if (a1 != iedge && hostSetInfoVecs.edges2Nodes_2[a1] == edge_start){H1 = a1;}
        else if (a1 != iedge && hostSetInfoVecs.edges2Nodes_1[a1] == edge_end){H2 = a1;}
        else if (a1 != iedge && hostSetInfoVecs.edges2Nodes_2[a1] == edge_end){H2 = a1;}

        if (b1 != iedge && hostSetInfoVecs.edges2Nodes_1[b1] == edge_start){H1 = b1;}
        else if (b1 != iedge && hostSetInfoVecs.edges2Nodes_2[b1] == edge_start){H1 = b1;}
        else if (b1 != iedge && hostSetInfoVecs.edges2Nodes_1[b1] == edge_end){H2 = b1;}
        else if (b1 != iedge && hostSetInfoVecs.edges2Nodes_2[b1] == edge_end){H2 = b1;}

        if (c1 != iedge && hostSetInfoVecs.edges2Nodes_1[c1] == edge_start){H1 = c1;}
        else if (c1 != iedge && hostSetInfoVecs.edges2Nodes_2[c1] == edge_start){H1 = c1;}
        else if (c1 != iedge && hostSetInfoVecs.edges2Nodes_1[c1] == edge_end){H2 = c1;}
        else if (c1 != iedge && hostSetInfoVecs.edges2Nodes_2[c1] == edge_end){H2 = c1;}
        
        if (a2 != iedge && hostSetInfoVecs.edges2Nodes_1[a2] == edge_start){T1 = a2;}
        else if (a2 != iedge && hostSetInfoVecs.edges2Nodes_2[a2] == edge_start){T1 = a2;}
        else if (a2 != iedge && hostSetInfoVecs.edges2Nodes_1[a2] == edge_end){T2 = a2;}
        else if (a2 != iedge && hostSetInfoVecs.edges2Nodes_2[a2] == edge_end){T2 = a2;}

        if (b2 != iedge && hostSetInfoVecs.edges2Nodes_1[b2] == edge_start){T1 = b2;}
        else if (b2 != iedge && hostSetInfoVecs.edges2Nodes_2[b2] == edge_start){T1 = b2;}
        else if (b2 != iedge && hostSetInfoVecs.edges2Nodes_1[b2] == edge_end){T2 = b2;}
        else if (b2 != iedge && hostSetInfoVecs.edges2Nodes_2[b2] == edge_end){T2 = b2;}

        if (c2 != iedge && hostSetInfoVecs.edges2Nodes_1[c2] == edge_start){T1 = c2;}
        else if (c2 != iedge && hostSetInfoVecs.edges2Nodes_2[c2] == edge_start){T1 = c2;}
        else if (c2 != iedge && hostSetInfoVecs.edges2Nodes_1[c2] == edge_end){T2 = c2;}
        else if (c2 != iedge && hostSetInfoVecs.edges2Nodes_2[c2] == edge_end){T2 = c2;}

        //std::cout<<"H1 = "<<H1<<std::endl;
        //std::cout<<"H2 = "<<H2<<std::endl;
        //std::cout<<"T1 = "<<T1<<std::endl;
        //std::cout<<"T2 = "<<T2<<std::endl;

        //Now search for the associated 

        int CANDIDATE1_1 = hostSetInfoVecs.triangles2Nodes_1[H0];
        int CANDIDATE1_2 = hostSetInfoVecs.triangles2Nodes_2[H0];
        int CANDIDATE1_3 = hostSetInfoVecs.triangles2Nodes_3[H0];
        int CANDIDATE2_1 = hostSetInfoVecs.triangles2Nodes_1[T0];
        int CANDIDATE2_2 = hostSetInfoVecs.triangles2Nodes_2[T0];
        int CANDIDATE2_3 = hostSetInfoVecs.triangles2Nodes_3[T0];
        
        if ((CANDIDATE1_1 != edge_start) 
            && (CANDIDATE1_1 != edge_end)) {
            HEAD = CANDIDATE1_1;
        }
        else if ((CANDIDATE1_2 != edge_start) && (CANDIDATE1_2 != edge_end)){HEAD = CANDIDATE1_2;}
        else if (CANDIDATE1_3 != edge_start && CANDIDATE1_3 != edge_end){HEAD = CANDIDATE1_3;}
        else {std::cout<<"head not set" <<std::endl;}

        if (CANDIDATE2_1 != edge_start && CANDIDATE2_1 != edge_end){TAIL = CANDIDATE2_1;}
        else if (CANDIDATE2_2 != edge_start && CANDIDATE2_2 != edge_end){TAIL = CANDIDATE2_2;}
        else if (CANDIDATE2_3 != edge_start && CANDIDATE2_3 != edge_end){TAIL = CANDIDATE2_3;}
        else {std::cout<<"tail not set" <<std::endl;}

        bool BAD_CHOICE = false;
        for (int q = 0; q < 4; q++){
            double qq;
            if (q == 0){
                qq = edge_start;
            }
            else if (q == 1){
                qq = edge_end;
            }
            else if (q == 2){
                qq = HEAD;
            }
            else if (q == 3){
                qq = TAIL;
            }
            int safe_flip1 = 0;
            if (hostSetInfoVecs.nndata1[qq] >= 0){safe_flip1 += 1;        }
            if (hostSetInfoVecs.nndata2[qq] >= 0){safe_flip1 += 1;        }
            if (hostSetInfoVecs.nndata3[qq] >= 0){safe_flip1 += 1;        }
            if (hostSetInfoVecs.nndata4[qq] >= 0){safe_flip1 += 1;        }
            if (hostSetInfoVecs.nndata5[qq] >= 0){safe_flip1 += 1;        }
            if (hostSetInfoVecs.nndata6[qq] >= 0){safe_flip1 += 1;        }
            if (hostSetInfoVecs.nndata7[qq] >= 0){safe_flip1 += 1;        }
            if (hostSetInfoVecs.nndata8[qq] >= 0){safe_flip1 += 1;        }
            if (hostSetInfoVecs.nndata9[qq] >= 0){safe_flip1 += 1;        }
           // if (hostSetInfoVecs.nndata10[qq] >= 0){safe_flip1 += 1;        }
           // if (hostSetInfoVecs.nndata11[qq] >= 0){safe_flip1 += 1;        }
           // if (hostSetInfoVecs.nndata12[qq] >= 0){safe_flip1 += 1;        }

            if (q == 0 && safe_flip1 == 4){
                BAD_CHOICE = true;
               // std::cout<<"BAD_CHOICE = "<<BAD_CHOICE<<std::endl;
               // std::cout<<"SAFE_FLIP = "<<safe_flip1<<std::endl;
                break;
            }
            else if (q == 1 && safe_flip1 == 4){
                BAD_CHOICE = true;
               // std::cout<<"BAD_CHOICE = "<<BAD_CHOICE<<std::endl;
        //std::cout<<"SAFE_FLIP = "<<safe_flip1<<std::endl;
                break;
            }
            else if (q == 2 && safe_flip1 == generalParams.safeguardthreshold){
                BAD_CHOICE = true;
               // std::cout<<"BAD_CHOICE = "<<BAD_CHOICE<<std::endl;
        //std::cout<<"SAFE_FLIP = "<<safe_flip1<<std::endl;
                break;
            }
            else if (q == 3 && safe_flip1 == generalParams.safeguardthreshold){
                BAD_CHOICE = true;
               // std::cout<<"BAD_CHOICE = "<<BAD_CHOICE<<std::endl;
        //std::cout<<"SAFE_FLIP = "<<safe_flip1<<std::endl;
                break;
            }

            //WE NEED ONE LAST CHECK TO AVOID OVERLAPPING EDGES TO OCCUR, THIS IS EXCLUSIVELY OBSERVED WITH A NARROW BUDDING NECK SHOWS UP IN THE SYSTEM
            if (hostSetInfoVecs.nndata1[TAIL] == HEAD){BAD_CHOICE = true; break; std::cout<<"OVERLAPPING EDGES OCCUR! AVOID!"<<std::endl;}
            else if (hostSetInfoVecs.nndata2[TAIL] == HEAD){BAD_CHOICE = true; break; std::cout<<"OVERLAPPING EDGES OCCUR! AVOID!"<<std::endl;}
            else if (hostSetInfoVecs.nndata3[TAIL] == HEAD){BAD_CHOICE = true; break; std::cout<<"OVERLAPPING EDGES OCCUR! AVOID!"<<std::endl;}
            else if (hostSetInfoVecs.nndata4[TAIL] == HEAD){BAD_CHOICE = true; break; std::cout<<"OVERLAPPING EDGES OCCUR! AVOID!"<<std::endl;}
            else if (hostSetInfoVecs.nndata5[TAIL] == HEAD){BAD_CHOICE = true; break; std::cout<<"OVERLAPPING EDGES OCCUR! AVOID!"<<std::endl;}
            else if (hostSetInfoVecs.nndata6[TAIL] == HEAD){BAD_CHOICE = true; break; std::cout<<"OVERLAPPING EDGES OCCUR! AVOID!"<<std::endl;}
            else if (hostSetInfoVecs.nndata7[TAIL] == HEAD){BAD_CHOICE = true; break; std::cout<<"OVERLAPPING EDGES OCCUR! AVOID!"<<std::endl;}
            else if (hostSetInfoVecs.nndata8[TAIL] == HEAD){BAD_CHOICE = true; break; std::cout<<"OVERLAPPING EDGES OCCUR! AVOID!"<<std::endl;}
            else if (hostSetInfoVecs.nndata9[TAIL] == HEAD){BAD_CHOICE = true; break; std::cout<<"OVERLAPPING EDGES OCCUR! AVOID!"<<std::endl;}
            if (hostSetInfoVecs.nndata1[HEAD] == TAIL){BAD_CHOICE = true; break; std::cout<<"OVERLAPPING EDGES OCCUR! AVOID!"<<std::endl;}
            else if (hostSetInfoVecs.nndata2[HEAD] == TAIL){BAD_CHOICE = true; break; std::cout<<"OVERLAPPING EDGES OCCUR! AVOID!"<<std::endl;}
            else if (hostSetInfoVecs.nndata3[HEAD] == TAIL){BAD_CHOICE = true; break; std::cout<<"OVERLAPPING EDGES OCCUR! AVOID!"<<std::endl;}
            else if (hostSetInfoVecs.nndata4[HEAD] == TAIL){BAD_CHOICE = true; break; std::cout<<"OVERLAPPING EDGES OCCUR! AVOID!"<<std::endl;}
            else if (hostSetInfoVecs.nndata5[HEAD] == TAIL){BAD_CHOICE = true; break; std::cout<<"OVERLAPPING EDGES OCCUR! AVOID!"<<std::endl;}
            else if (hostSetInfoVecs.nndata6[HEAD] == TAIL){BAD_CHOICE = true; break; std::cout<<"OVERLAPPING EDGES OCCUR! AVOID!"<<std::endl;}
            else if (hostSetInfoVecs.nndata7[HEAD] == TAIL){BAD_CHOICE = true; break; std::cout<<"OVERLAPPING EDGES OCCUR! AVOID!"<<std::endl;}
            else if (hostSetInfoVecs.nndata8[HEAD] == TAIL){BAD_CHOICE = true; break; std::cout<<"OVERLAPPING EDGES OCCUR! AVOID!"<<std::endl;}
            else if (hostSetInfoVecs.nndata9[HEAD] == TAIL){BAD_CHOICE = true; break; std::cout<<"OVERLAPPING EDGES OCCUR! AVOID!"<<std::endl;}

        }
        

        if (BAD_CHOICE == false){//(safe_flip1 < generalParams.safeguardthreshold && safe_flip2 < generalParams.safeguardthreshold){



            //int temp_edges2Nodes_2 = HEAD;
            ////std::cout<<"head tail in loop: "<< HEAD << " "<< TAIL <<std::endl;
            //The small subsystem we will be working with is
            //          
            //           edge_start
            //    T10    *   |    *     H10
            //         T1    |     H1
            //        *      |       *
            //    TAIL   T0  |  H0    HEAD
            //        *      |       *
            //         T2    |     H2
            //    T20    *   v    *     H20
            //            edge_end
            //
            //H10 is the triangle sharing the same edge H1 with triangle H0.

            //energy E_0 calculation
            //Since nodes are NOT moved and area is not changed, we only need to calculate 
            //linear spring energy and bending energy.
            //Furthermore, since linear spring energy will only be nontrivial for the edge swapped,
            //we can condense the linear spring energy computation to only one edge.
            //Bending energy is more complicated due to the need of unit normals.
            
            

            std::vector<int> edges_iteration(5);
            edges_iteration[0] = iedge;
            edges_iteration[1] = H1;
            edges_iteration[2] = H2;
            edges_iteration[3] = T1;
            edges_iteration[4] = T2;

            double rep_0;
            double rep_energy_TAIL_HEAD;
            double R_TAIL_HEAD = sqrt((hostSetInfoVecs.nodeLocX[TAIL] - hostSetInfoVecs.nodeLocX[HEAD])*(hostSetInfoVecs.nodeLocX[TAIL] - hostSetInfoVecs.nodeLocX[HEAD]) +
                                    (hostSetInfoVecs.nodeLocY[TAIL] - hostSetInfoVecs.nodeLocY[HEAD])*(hostSetInfoVecs.nodeLocY[TAIL] - hostSetInfoVecs.nodeLocY[HEAD]) +
                                    (hostSetInfoVecs.nodeLocZ[TAIL] - hostSetInfoVecs.nodeLocZ[HEAD])*(hostSetInfoVecs.nodeLocZ[TAIL] - hostSetInfoVecs.nodeLocZ[HEAD]));
            if (R_TAIL_HEAD < generalParams.abs_Rmin){
                //rep_energy_TAIL_HEAD = linearSpringInfoVecs.spring_constant_rep1*(1.0-exp(-linearSpringInfoVecs.spring_constant_rep2*(R_TAIL_HEAD - generalParams.abs_Rmin)))*(1.0-exp(-linearSpringInfoVecs.spring_constant_rep2*(R_TAIL_HEAD - generalParams.abs_Rmin)));
                if (hostSetInfoVecs.nodes_in_upperhem[TAIL] == 1 && hostSetInfoVecs.nodes_in_upperhem[HEAD] == 1){
                    rep_energy_TAIL_HEAD = (linearSpringInfoVecs.spring_constant_weak/2.0)*(R_TAIL_HEAD - generalParams.abs_Rmin)*(R_TAIL_HEAD - generalParams.abs_Rmin);
                }
                else{
                    rep_energy_TAIL_HEAD = (linearSpringInfoVecs.spring_constant/2.0)*(R_TAIL_HEAD - generalParams.abs_Rmin)*(R_TAIL_HEAD - generalParams.abs_Rmin);
                }
                }
            else{
                rep_energy_TAIL_HEAD = 0.0;
                }
            int node_index_H1, node_index_H2, node_index_T1, node_index_T2;
            rep_0 = rep_energy_TAIL_HEAD;
            /*if (generalParams.edges_in_upperhem[edges_iteration[0]] == 1){
                        linear_spring_constant = linearSpringInfoVecs.spring_constant_weak;
                    }
                    else if (generalParams.edges_in_upperhem[edges_iteration[0]] == 0){
                        linear_spring_constant = (linearSpringInfoVecs.spring_constant_weak + linearSpringInfoVecs.spring_constant)/2.0;
                    }
                    else{
                        linear_spring_constant = linearSpringInfoVecs.spring_constant;
            }*/
            if (generalParams.SCALE_TYPE == 0){
                linear_spring_constant = linearSpringInfoVecs.spring_constant*(1.0 - ((1.0/sqrt(2*3.14159*generalParams.gausssigma))*exp(-pow(hostSetInfoVecs.scaling_per_edge[edges_iteration[0]],2.0)/generalParams.gausssigma)));
                if (linear_spring_constant < linearSpringInfoVecs.spring_constant_weak){linear_spring_constant = linearSpringInfoVecs.spring_constant_weak;};
            }
            else if (generalParams.SCALE_TYPE == 1){
                // linear_spring_constant = (linearSpringInfoVecs.spring_constant_weak*2.0)*pow(hostSetInfoVecs.scaling_per_edge[edges_iteration[0]], generalParams.scaling_pow) + 
                //                         linearSpringInfoVecs.spring_constant_weak*(1.0 - pow(hostSetInfoVecs.scaling_per_edge[edges_iteration[0]],generalParams.scaling_pow));
                linear_spring_constant = linearSpringInfoVecs.spring_constant*pow(hostSetInfoVecs.scaling_per_edge[edges_iteration[0]], generalParams.scaling_pow) + 
                                        linearSpringInfoVecs.spring_constant_weak*(1.0 - pow(hostSetInfoVecs.scaling_per_edge[edges_iteration[0]],generalParams.scaling_pow));
            }
            else if (generalParams.SCALE_TYPE == 2){
                linear_spring_constant = linearSpringInfoVecs.spring_constant - 
                                    (linearSpringInfoVecs.spring_constant - linearSpringInfoVecs.spring_constant_weak)*
                                    hostSetInfoVecs.scaling_per_edge[edges_iteration[0]];
                if (linear_spring_constant < linearSpringInfoVecs.spring_constant_weak){linear_spring_constant = linearSpringInfoVecs.spring_constant_weak;}
            }
            else if (generalParams.SCALE_TYPE == 3){
                if (hostSetInfoVecs.edges_in_upperhem[edges_iteration[0]] == 1){// && hostSetInfoVecs.edges_in_tip[edges_iteration[0]] == 1){
                        linear_spring_constant = linearSpringInfoVecs.spring_constant_weak;
                    }
                    else if (hostSetInfoVecs.edges_in_upperhem[edges_iteration[0]] == 0){// && hostSetInfoVecs.edges_in_tip[edges_iteration[0]] == 0){
                        // linear_spring_constant = (linearSpringInfoVecs.spring_constant_weak*2.0);
                        linear_spring_constant = (linearSpringInfoVecs.spring_constant_weak + linearSpringInfoVecs.spring_constant)/2.0;
                    }
                    else{
                        linear_spring_constant = linearSpringInfoVecs.spring_constant;
                    }
            }
            else if (generalParams.SCALE_TYPE == 4){
                //double scaling = 0.0;//linearSpringInfoVecs.spring_constant_weak/linearSpringInfoVecs.spring_constant;
                //linear_spring_constant = linearSpringInfoVecs.spring_constant*((1.0/(1.0+pow(generalParams.hilleqnconst/hostSetInfoVecs.scaling_per_edge[edges_iteration[0]], generalParams.hilleqnpow)))*(1-scaling) + scaling);
			    double spectrum = linearSpringInfoVecs.spring_constant - linearSpringInfoVecs.spring_constant_weak;
                linear_spring_constant = linearSpringInfoVecs.spring_constant_weak + ((1.0/(1.0+pow(generalParams.hilleqnconst/hostSetInfoVecs.scaling_per_edge[edges_iteration[0]], generalParams.hilleqnpow)))*spectrum);
                if (linear_spring_constant < linearSpringInfoVecs.spring_constant_weak){linear_spring_constant = linearSpringInfoVecs.spring_constant_weak;}
		    }
            
                int wrong1, wrong2, wrong3;
                for (int j = 0; j < 5; j++){
                   /* if (generalParams.edges_in_upperhem[edges_iteration[j]] == 1){
                        bend_spring_constant = bendingTriangleInfoVecs.spring_constant_weak;
                    }
                    else if (generalParams.edges_in_upperhem[edges_iteration[j]] == 0){
                        bend_spring_constant = (bendingTriangleInfoVecs.spring_constant_weak + bendingTriangleInfoVecs.spring_constant)/2.0;
                    }
                    else{
                        bend_spring_constant = bendingTriangleInfoVecs.spring_constant;
                    }*/
                    if (generalParams.SCALE_TYPE == 0){
                        bend_spring_constant = bendingTriangleInfoVecs.spring_constant*(1.0 - ((1.0/sqrt(2*3.14159*generalParams.gausssigma))*exp(-pow(hostSetInfoVecs.scaling_per_edge[edges_iteration[j]],2.0)/generalParams.gausssigma)));
                        if (bend_spring_constant < bendingTriangleInfoVecs.spring_constant_weak){bend_spring_constant = bendingTriangleInfoVecs.spring_constant_weak;};
                       /* if (generalParams.edges_in_upperhem[edges_iteration[j]] == 1){
                            bend_spring_constant = bendingTriangleInfoVecs.spring_constant_weak;
                        }
                        else if (generalParams.edges_in_upperhem[edges_iteration[j]] == 0){
                            bend_spring_constant = (bendingTriangleInfoVecs.spring_constant_weak + bendingTriangleInfoVecs.spring_constant)/2.0;
                        }
                        else{
                            bend_spring_constant = bendingTriangleInfoVecs.spring_constant;
                        }*/
                    }
                    else if (generalParams.SCALE_TYPE == 1){
                        bend_spring_constant = (bendingTriangleInfoVecs.spring_constant_weak*4.0)*pow(hostSetInfoVecs.scaling_per_edge[edges_iteration[j]],generalParams.scaling_pow) +
                                            bendingTriangleInfoVecs.spring_constant_weak*(1.0 - pow(hostSetInfoVecs.scaling_per_edge[edges_iteration[j]],generalParams.scaling_pow));
                        // bend_spring_constant = bendingTriangleInfoVecs.spring_constant*pow(hostSetInfoVecs.scaling_per_edge[edges_iteration[j]],generalParams.scaling_pow) +
                        //                     bendingTriangleInfoVecs.spring_constant_weak*(1.0 - pow(hostSetInfoVecs.scaling_per_edge[edges_iteration[j]],generalParams.scaling_pow));
                    }   
                    else if (generalParams.SCALE_TYPE == 2){
                        bend_spring_constant = bendingTriangleInfoVecs.spring_constant_weak - 
                        (bendingTriangleInfoVecs.spring_constant - bendingTriangleInfoVecs.spring_constant_weak)*
                                            hostSetInfoVecs.scaling_per_edge[edges_iteration[j]];
                        if (bend_spring_constant < bendingTriangleInfoVecs.spring_constant_weak){bend_spring_constant = bendingTriangleInfoVecs.spring_constant_weak;}
                    }
                    else if (generalParams.SCALE_TYPE == 3){
                        if (hostSetInfoVecs.edges_in_upperhem[edges_iteration[j]] == 1 ){//&& hostSetInfoVecs.edges_in_tip[edges_iteration[j]] == 1){
                        bend_spring_constant = bendingTriangleInfoVecs.spring_constant_weak;
                        preferred_angle = bendingTriangleInfoVecs.initial_angle_bud;
                        }
                        else if (hostSetInfoVecs.edges_in_upperhem[edges_iteration[j]] == 0){// 1 && hostSetInfoVecs.edges_in_tip[edges_iteration[j]] != 1){
                            //bend_spring_constant = (bendingTriangleInfoVecs.spring_constant_weak*2.0);
                            bend_spring_constant = (bendingTriangleInfoVecs.spring_constant_weak + bendingTriangleInfoVecs.spring_constant)/2.0;
                            preferred_angle = (bendingTriangleInfoVecs.initial_angle_bud + bendingTriangleInfoVecs.initial_angle)/2.0;
                        }
                        else{
                            bend_spring_constant = bendingTriangleInfoVecs.spring_constant;
                            preferred_angle = bendingTriangleInfoVecs.initial_angle;
                        }
                    }
                    else if (generalParams.SCALE_TYPE == 4){
                        //double scaling = 0.0;//bendingTriangleInfoVecs.spring_constant_weak/bendingTriangleInfoVecs.spring_constant;
                        //bend_spring_constant = bendingTriangleInfoVecs.spring_constant*((1.0/(1.0+pow(generalParams.hilleqnconst/hostSetInfoVecs.scaling_per_edge[edges_iteration[j]], generalParams.hilleqnpow)))*(1-scaling) + scaling);
                        double spectrum = bendingTriangleInfoVecs.spring_constant - bendingTriangleInfoVecs.spring_constant_weak;
                        bend_spring_constant = bendingTriangleInfoVecs.spring_constant_weak + ((1.0/(1.0+pow(generalParams.hilleqnconst/hostSetInfoVecs.scaling_per_edge[edges_iteration[j]], generalParams.hilleqnpow)))*spectrum);
                        if (bend_spring_constant < bendingTriangleInfoVecs.spring_constant_weak){bend_spring_constant = bendingTriangleInfoVecs.spring_constant_weak;}
                    }
                    int Tri1 = hostSetInfoVecs.edges2Triangles_1[edges_iteration[j]];//index of the 1st triangle
                    int Tri2 = hostSetInfoVecs.edges2Triangles_2[edges_iteration[j]];
                    //int id_k = hostSetInfoVecs.edges2Nodes_1[edges_iteration[j]];
                    //int id_i = hostSetInfoVecs.edges2Nodes_2[edges_iteration[j]];

                    double N1vec1x, N1vec1y, N1vec1z, N1vec2x, N1vec2y, N1vec2z;
                    double N2vec1x, N2vec1y, N2vec1z, N2vec2x, N2vec2y, N2vec2z;
                    double N1_x, N1_y, N1_z, N2_x, N2_y, N2_z;
                    
                    if (Tri1 != Tri2) {
                        int Tri1_n1 = hostSetInfoVecs.triangles2Nodes_1[Tri1];
                        if (j == 0){
                            P0x_vol1 = hostSetInfoVecs.nodeLocX[Tri1_n1];
                            P0y_vol1 = hostSetInfoVecs.nodeLocY[Tri1_n1];
                            P0z_vol1 = hostSetInfoVecs.nodeLocZ[Tri1_n1];
                        }
                        int Tri1_n2 = hostSetInfoVecs.triangles2Nodes_2[Tri1];
                        int Tri1_n3 = hostSetInfoVecs.triangles2Nodes_3[Tri1];
                        N1vec1x = hostSetInfoVecs.nodeLocX[Tri1_n2] - hostSetInfoVecs.nodeLocX[Tri1_n1];
                        N1vec1y = hostSetInfoVecs.nodeLocY[Tri1_n2] - hostSetInfoVecs.nodeLocY[Tri1_n1];
                        N1vec1z = hostSetInfoVecs.nodeLocZ[Tri1_n2] - hostSetInfoVecs.nodeLocZ[Tri1_n1];
                        N1vec2x = hostSetInfoVecs.nodeLocX[Tri1_n3] - hostSetInfoVecs.nodeLocX[Tri1_n1];
                        N1vec2y = hostSetInfoVecs.nodeLocY[Tri1_n3] - hostSetInfoVecs.nodeLocY[Tri1_n1];
                        N1vec2z = hostSetInfoVecs.nodeLocZ[Tri1_n3] - hostSetInfoVecs.nodeLocZ[Tri1_n1];
                        //std::vector<double> N1(3);
                        N1_x = N1vec1y*N1vec2z - N1vec2y*N1vec1z;
                        N1_y = -(N1vec1x*N1vec2z - N1vec2x*N1vec1z);
                        N1_z = N1vec1x*N1vec2y - N1vec2x*N1vec1y;
                        double nN1 = sqrt(pow(N1_x,2)+pow(N1_y,2)+pow(N1_z,2));
                        ////std::cout<<"nN1 = "<<nN1<<std::endl;

                        int Tri2_n1 = hostSetInfoVecs.triangles2Nodes_1[Tri2];
                        int Tri2_n2 = hostSetInfoVecs.triangles2Nodes_2[Tri2];
                        int Tri2_n3 = hostSetInfoVecs.triangles2Nodes_3[Tri2];
                        if (j == 0){
                            P0x_vol2 = hostSetInfoVecs.nodeLocX[Tri2_n1];
                            P0y_vol2 = hostSetInfoVecs.nodeLocY[Tri2_n1];
                            P0z_vol2 = hostSetInfoVecs.nodeLocZ[Tri2_n1];
                        }
                        N2vec1x = hostSetInfoVecs.nodeLocX[Tri2_n2] - hostSetInfoVecs.nodeLocX[Tri2_n1];
                        N2vec1y = hostSetInfoVecs.nodeLocY[Tri2_n2] - hostSetInfoVecs.nodeLocY[Tri2_n1];
                        N2vec1z = hostSetInfoVecs.nodeLocZ[Tri2_n2] - hostSetInfoVecs.nodeLocZ[Tri2_n1];
                        N2vec2x = hostSetInfoVecs.nodeLocX[Tri2_n3] - hostSetInfoVecs.nodeLocX[Tri2_n1];
                        N2vec2y = hostSetInfoVecs.nodeLocY[Tri2_n3] - hostSetInfoVecs.nodeLocY[Tri2_n1];
                        N2vec2z = hostSetInfoVecs.nodeLocZ[Tri2_n3] - hostSetInfoVecs.nodeLocZ[Tri2_n1];
                        //std::vector<double> N2(3);
                        N2_x = N2vec1y*N2vec2z - N2vec2y*N2vec1z;
                        N2_y = -(N2vec1x*N2vec2z - N2vec2x*N2vec1z);
                        N2_z = N2vec1x*N2vec2y - N2vec2x*N2vec1y; 
                        double nN2 = sqrt(pow(N2_x,2)+pow(N2_y,2)+pow(N2_z,2));;
                        ////std::cout<<"nN2 = "<<nN2<<std::endl;

                        double cosAngle = (N1_x*N2_x + N1_y*N2_y + N1_z*N2_z)/(nN1*nN2);
                        ////std::cout<<"dotproduct = "<<cosAngle<<std::endl;
                    
        
                        if (cosAngle > 1.0) {
                            cosAngle = 1.0;
                        }
                        else if (cosAngle < -1.0){
                            cosAngle = -1.0;
                        }

                        double theta_current = acos( cosAngle );
                        
                        
                        double local_energy = bend_spring_constant * (1 - cos(theta_current - preferred_angle) );
                        
                        //bendingTriangleInfoVecs.spring_constant * (1 - cos(theta_current - bendingTriangleInfoVecs.initial_angle) );
                        temp_bend = temp_bend + local_energy;
                        if (j == 1){
                            wrong1 = HEAD;
                            wrong2 = edge_start;
                            wrong3 = edge_end;
                        }
                        if (j == 2){
                            wrong1 = edge_end;
                            wrong2 = HEAD;
                            wrong3 = edge_start;
                        }
                        if (j == 3){
                            wrong1 = edge_start;
                            wrong2 = TAIL;
                            wrong3 = edge_end;
                        }
                        if (j == 4){
                            wrong1 = TAIL;
                            wrong2 = edge_end;
                            wrong3 = edge_start;
                        }
                        if (Tri1_n1 != wrong1 && Tri1_n1 != wrong2 && Tri1_n1 != wrong3){
                            if (j == 1){
                                node_index_H1 = Tri1_n1;
                            }
                            else if (j == 2){
                                node_index_H2 = Tri1_n1;
                            }
                            else if (j == 3){
                                node_index_T1 = Tri1_n1;
                            }
                            else if (j == 4){
                                node_index_T2 = Tri1_n1;
                            }
                        }
                        else if (Tri1_n2 != wrong1 && Tri1_n2 != wrong2 && Tri1_n2 != wrong3){
                            if (j == 1){
                                node_index_H1 = Tri1_n2;
                            }
                            else if (j == 2){
                                node_index_H2 = Tri1_n2;
                            }
                            else if (j == 3){
                                node_index_T1 = Tri1_n2;
                            }
                            else if (j == 4){
                                node_index_T2 = Tri1_n2;
                            }
                        }
                        else if (Tri1_n3 != wrong1 && Tri1_n3 != wrong2 && Tri1_n3 != wrong3){

                            if (j == 1){
                                node_index_H1 = Tri1_n3;
                            }
                            else if (j == 2){
                                node_index_H2 = Tri1_n3;
                            }
                            else if (j == 3){
                                node_index_T1 = Tri1_n3;
                            }
                            else if (j == 4){
                                node_index_T2 = Tri1_n3;
                            }
                        }
                        else if (Tri2_n1 != wrong1 && Tri2_n1 != wrong2 && Tri2_n1 != wrong3){
                            
                            if (j == 1){
                                node_index_H1 = Tri2_n1;
                            }
                            else if (j == 2){
                                node_index_H2 = Tri2_n1;
                            }
                            else if (j == 3){
                                node_index_T1 = Tri2_n1;
                            }
                            else if (j == 4){
                                node_index_T2 = Tri2_n1;
                            }
                        }
                        else if (Tri2_n2 != wrong1 && Tri2_n2 != wrong2 && Tri2_n2 != wrong3){
                            
                            if (j == 1){
                                node_index_H1 = Tri2_n2;
                            }
                            else if (j == 2){
                                node_index_H2 = Tri2_n2;
                            }
                            else if (j == 3){
                                node_index_T1 = Tri2_n2;
                            }
                            else if (j == 4){
                                node_index_T2 = Tri2_n2;
                            }
                        }
                        else if (Tri2_n3 != wrong1 && Tri2_n3 != wrong2 && Tri2_n3 != wrong3){
                            
                            if (j == 1){
                                node_index_H1 = Tri2_n3;
                            }
                            else if (j == 2){
                                node_index_H2 = Tri2_n3;
                            }
                            else if (j == 3){
                                node_index_T1 = Tri2_n3;
                            }
                            else if (j == 4){
                                node_index_T2 = Tri2_n3;
                            }
                        }

                        
                        /*//std::cout<<"bending energy "<<local_energy<<std::endl;
                        for (int COUNT = 0; COUNT < 3; COUNT++){
                        //std::cout<<"unit normal 1 = "<<N1[COUNT]<<std::endl;
                        //std::cout<<"unit normal 2 = "<<N2[COUNT]<<std::endl;}
                        //std::cout<<"angle "<<theta_current<<std::endl;*/
                        if (j == 0){
                            N1x_vol = N1_x/nN1;
                            N1y_vol = N1_y/nN1;
                            N1z_vol = N1_z/nN1;
                            N2x_vol = N2_x/nN2;
                            N2y_vol = N2_y/nN2;
                            N2z_vol = N2_z/nN2;
                        }
                    }
                }

            int nb1, nb2, nb3, nb4, nb5, nb6, nb7, nb8, nb9;
            int midpt, endpt1, endpt2, endpt3, falsept1, falsept2;
            int true_endpoints;
            for (int t = 0; t < 4; t++){    

                if (t == 0){
                    midpt = edge_start;
                    endpt1 = HEAD;
                    endpt2 = TAIL;
                    endpt3 = edge_end;
                    true_endpoints = 3;
                    falsept1 = node_index_H1;
                    falsept2 = node_index_T1;
                }
                else if (t == 1){
                    midpt = HEAD;
                    endpt1 = edge_start;
                    endpt2 = edge_end;
                    endpt3 = TAIL;
                    true_endpoints = 2;
                    falsept1 = node_index_H1;
                    falsept2 = node_index_H2;
                }
                else if (t == 2){
                    midpt = edge_end;
                    endpt1 = HEAD;
                    endpt2 = edge_start;
                    endpt3 = TAIL;
                    true_endpoints = 3;
                    falsept1 = node_index_H2;
                    falsept2 = node_index_T2;
                }
                else if (t == 3){
                    midpt = TAIL;
                    endpt1 = edge_start;
                    endpt2 = edge_end;
                    endpt3 = HEAD;
                    true_endpoints = 2;
                    falsept1 = node_index_T1;
                    falsept2 = node_index_T2;
                }
                nb1 = hostSetInfoVecs.nndata1[midpt];
                nb2 = hostSetInfoVecs.nndata2[midpt];
                nb3 = hostSetInfoVecs.nndata3[midpt];
                nb4 = hostSetInfoVecs.nndata4[midpt];
                nb5 = hostSetInfoVecs.nndata5[midpt];
                nb6 = hostSetInfoVecs.nndata6[midpt];
                nb7 = hostSetInfoVecs.nndata7[midpt];
                nb8 = hostSetInfoVecs.nndata8[midpt];
                nb9 = hostSetInfoVecs.nndata9[midpt];
                for (int y = 0; y < 9; y++){
                    int startpt;
                    if (y == 0 && nb1>= 0 && nb1 != endpt1 && nb1 != endpt2 && nb1 != endpt3 && nb1 != falsept1 && nb1 != falsept2){startpt = nb1;}
                    else if (y == 1 && nb2>= 0 && nb2 != endpt1 && nb2 != endpt2 && nb2 != endpt3 && nb2 != falsept1 && nb2 != falsept2){startpt = nb2;}
                    else if (y == 2 && nb3>= 0 && nb3 != endpt1 && nb3 != endpt2 && nb3 != endpt3 && nb3 != falsept1 && nb3 != falsept2){startpt = nb3;}
                    else if (y == 3 && nb4>= 0 && nb4 != endpt1 && nb4 != endpt2 && nb4 != endpt3 && nb4 != falsept1 && nb4 != falsept2){startpt = nb4;}
                    else if (y == 4 && nb5>= 0 && nb5 != endpt1 && nb5 != endpt2 && nb5 != endpt3 && nb5 != falsept1 && nb5 != falsept2){startpt = nb5;}
                    else if (y == 5 && nb6>= 0 && nb6 != endpt1 && nb6 != endpt2 && nb6 != endpt3 && nb6 != falsept1 && nb6 != falsept2){startpt = nb6;}
                    else if (y == 6 && nb7>= 0 && nb7 != endpt1 && nb7 != endpt2 && nb7 != endpt3 && nb7 != falsept1 && nb7 != falsept2){startpt = nb7;}
                    else if (y == 7 && nb8>= 0 && nb8 != endpt1 && nb8 != endpt2 && nb8 != endpt3 && nb8 != falsept1 && nb8 != falsept2){startpt = nb8;}
                    else if (y == 8 && nb9>= 0 && nb9 != endpt1 && nb9 != endpt2 && nb9 != endpt3 && nb9 != falsept1 && nb9 != falsept2){startpt = nb9;}
                    else{continue;}

                    if (true_endpoints == 3){
                        for (int h = 0; h < 3; h++){
                            int aa;
                            if (h==0){aa = endpt1;}
                            else if (h==1){aa = endpt2;}
                            else if (h==2){aa = endpt3;}
                            R_TAIL_HEAD = sqrt((hostSetInfoVecs.nodeLocX[startpt] - hostSetInfoVecs.nodeLocX[aa])*(hostSetInfoVecs.nodeLocX[startpt] - hostSetInfoVecs.nodeLocX[aa]) +
                                    (hostSetInfoVecs.nodeLocY[startpt] - hostSetInfoVecs.nodeLocY[aa])*(hostSetInfoVecs.nodeLocY[startpt] - hostSetInfoVecs.nodeLocY[aa]) +
                                    (hostSetInfoVecs.nodeLocZ[startpt] - hostSetInfoVecs.nodeLocZ[aa])*(hostSetInfoVecs.nodeLocZ[startpt] - hostSetInfoVecs.nodeLocZ[aa]));
                            if (R_TAIL_HEAD < generalParams.abs_Rmin){
                                //rep_energy_TAIL_HEAD = linearSpringInfoVecs.spring_constant_rep1*(1.0-exp(-linearSpringInfoVecs.spring_constant_rep2*(R_TAIL_HEAD - generalParams.abs_Rmin)))*
                                //                        (1.0-exp(-linearSpringInfoVecs.spring_constant_rep2*(R_TAIL_HEAD - generalParams.abs_Rmin)));
                                if (hostSetInfoVecs.nodes_in_upperhem[startpt] == 1 && hostSetInfoVecs.nodes_in_upperhem[aa] == 1){
                                    rep_energy_TAIL_HEAD = (linearSpringInfoVecs.spring_constant_weak/2.0)*(R_TAIL_HEAD - generalParams.abs_Rmin)*(R_TAIL_HEAD - generalParams.abs_Rmin);
                                }
                                else{
                                    rep_energy_TAIL_HEAD = (linearSpringInfoVecs.spring_constant/2.0)*(R_TAIL_HEAD - generalParams.abs_Rmin)*(R_TAIL_HEAD - generalParams.abs_Rmin);
                                }
                            }
                            else{rep_energy_TAIL_HEAD = 0.0;}
                            rep_0 += rep_energy_TAIL_HEAD;
                        }
                    }
                    else if (true_endpoints == 2){
                        for (int h = 0; h < 2; h++){
                            int aa;
                            if (h==0){aa = endpt1;}
                            else if (h==1){aa = endpt2;}
                            R_TAIL_HEAD = sqrt((hostSetInfoVecs.nodeLocX[startpt] - hostSetInfoVecs.nodeLocX[aa])*(hostSetInfoVecs.nodeLocX[startpt] - hostSetInfoVecs.nodeLocX[aa]) +
                                    (hostSetInfoVecs.nodeLocY[startpt] - hostSetInfoVecs.nodeLocY[aa])*(hostSetInfoVecs.nodeLocY[startpt] - hostSetInfoVecs.nodeLocY[aa]) +
                                    (hostSetInfoVecs.nodeLocZ[startpt] - hostSetInfoVecs.nodeLocZ[aa])*(hostSetInfoVecs.nodeLocZ[startpt] - hostSetInfoVecs.nodeLocZ[aa]));
                            if (R_TAIL_HEAD < generalParams.abs_Rmin){
                               // rep_energy_TAIL_HEAD = linearSpringInfoVecs.spring_constant_rep1*(1.0-exp(-linearSpringInfoVecs.spring_constant_rep2*(R_TAIL_HEAD - generalParams.abs_Rmin)))*
                               //                         (1.0-exp(-linearSpringInfoVecs.spring_constant_rep2*(R_TAIL_HEAD - generalParams.abs_Rmin)));
                               if (hostSetInfoVecs.nodes_in_upperhem[startpt] == 1 && hostSetInfoVecs.nodes_in_upperhem[aa] == 1){
                                    rep_energy_TAIL_HEAD = (linearSpringInfoVecs.spring_constant_weak/2.0)*(R_TAIL_HEAD - generalParams.abs_Rmin)*(R_TAIL_HEAD - generalParams.abs_Rmin);
                                }
                                else{
                                    rep_energy_TAIL_HEAD = (linearSpringInfoVecs.spring_constant/2.0)*(R_TAIL_HEAD - generalParams.abs_Rmin)*(R_TAIL_HEAD - generalParams.abs_Rmin);
                                }
                            }
                            else{rep_energy_TAIL_HEAD = 0.0;}
                            rep_0 += rep_energy_TAIL_HEAD;
                        }
                    }
                }
                
                //NOW WE UTILIZE THIS LOOP TO CALCULATE THE REPULSION ENERGY ASSOCIATED WITH THE FALSEPTS (I.E. THE ONES ASSOCIATED WITH H10, T10, H20, T20 BUT NOT IN THE SUBSYSTEM DEPECTED ABOVE)
                //int oiu1, oiu2, oiu3;
                // if (t == 0){oiu1 = node_index_H1; oiu2 = edge_end; oiu3 = TAIL;}
                // else if (t == 1){oiu1 = node_index_H2; oiu2 = edge_start; oiu3 = TAIL;}
                // else if (t == 2){oiu1 = node_index_T1; oiu2 = edge_end; oiu3 = HEAD;}
                // else if (t == 3){oiu1 = node_index_T2; oiu2 = edge_start; oiu3 = HEAD;}
                // R_TAIL_HEAD = sqrt((hostSetInfoVecs.nodeLocX[oiu1] - hostSetInfoVecs.nodeLocX[oiu2])*(hostSetInfoVecs.nodeLocX[oiu1] - hostSetInfoVecs.nodeLocX[oiu2]) +
                //         (hostSetInfoVecs.nodeLocY[oiu1] - hostSetInfoVecs.nodeLocY[oiu2])*(hostSetInfoVecs.nodeLocY[oiu1] - hostSetInfoVecs.nodeLocY[oiu2]) +
                //         (hostSetInfoVecs.nodeLocZ[oiu1] - hostSetInfoVecs.nodeLocZ[oiu2])*(hostSetInfoVecs.nodeLocZ[oiu1] - hostSetInfoVecs.nodeLocZ[oiu2]));
                // if (R_TAIL_HEAD < generalParams.abs_Rmin){
                //     rep_energy_TAIL_HEAD = linearSpringInfoVecs.spring_constant_rep1*(1.0-exp(-linearSpringInfoVecs.spring_constant_rep2*(R_TAIL_HEAD - generalParams.abs_Rmin)))*
                //                             (1.0-exp(-linearSpringInfoVecs.spring_constant_rep2*(R_TAIL_HEAD - generalParams.abs_Rmin)));
                // }
                // else{rep_energy_TAIL_HEAD = 0.0;}
                // rep_0 += rep_energy_TAIL_HEAD;

                // R_TAIL_HEAD = sqrt((hostSetInfoVecs.nodeLocX[oiu1] - hostSetInfoVecs.nodeLocX[oiu3])*(hostSetInfoVecs.nodeLocX[oiu1] - hostSetInfoVecs.nodeLocX[oiu3]) +
                //         (hostSetInfoVecs.nodeLocY[oiu1] - hostSetInfoVecs.nodeLocY[oiu3])*(hostSetInfoVecs.nodeLocY[oiu1] - hostSetInfoVecs.nodeLocY[oiu3]) +
                //         (hostSetInfoVecs.nodeLocZ[oiu1] - hostSetInfoVecs.nodeLocZ[oiu3])*(hostSetInfoVecs.nodeLocZ[oiu1] - hostSetInfoVecs.nodeLocZ[oiu3]));
                // if (R_TAIL_HEAD < generalParams.abs_Rmin){
                //     rep_energy_TAIL_HEAD = linearSpringInfoVecs.spring_constant_rep1*(1.0-exp(-linearSpringInfoVecs.spring_constant_rep2*(R_TAIL_HEAD - generalParams.abs_Rmin)))*
                //                             (1.0-exp(-linearSpringInfoVecs.spring_constant_rep2*(R_TAIL_HEAD - generalParams.abs_Rmin)));
                // }
                // else{rep_energy_TAIL_HEAD = 0.0;}
                // rep_0 += rep_energy_TAIL_HEAD;
            }

            

            double bend_0 = temp_bend;
            //
            double linear_0;
            double DISTANCE = sqrt(
                pow(hostSetInfoVecs.nodeLocX[edge_end] - hostSetInfoVecs.nodeLocX[edge_start], 2.0) + 
                pow(hostSetInfoVecs.nodeLocY[edge_end] - hostSetInfoVecs.nodeLocY[edge_start], 2.0) + 
                pow(hostSetInfoVecs.nodeLocZ[edge_end] - hostSetInfoVecs.nodeLocZ[edge_start], 2.0));
            //if (DISTANCE < generalParams.abs_Rmin){
            //    linear_0 = (linearSpringInfoVecs.spring_constant/2)*(DISTANCE - generalParams.Rmin)*
            //	(DISTANCE - generalParams.Rmin);// + 
                //(linearSpringInfoVecs.spring_constant_rep1/2)*(DISTANCE - generalParams.abs_Rmin)*
                //(DISTANCE - generalParams.abs_Rmin);
                //linearSpringInfoVecs.spring_constant_rep1*(1-exp(-linearSpringInfoVecs.spring_constant_rep2*(DISTANCE-generalParams.abs_Rmin)))*
                //(1-exp(-linearSpringInfoVecs.spring_constant_rep2*(DISTANCE-generalParams.abs_Rmin)));
            //}
            //else if (DISTANCE != generalParams.Rmin){
            //    linear_0 = (linearSpringInfoVecs.spring_constant/2)*(DISTANCE - generalParams.Rmin)*
            //	(DISTANCE - generalParams.Rmin);
            //}
            //else{
                // linear_0 = (linear_spring_constant/(2.0*generalParams.Rmin*generalParams.Rmin))*(DISTANCE - generalParams.Rmin)*
                // (DISTANCE - generalParams.Rmin);
                // if (hostSetInfoVecs.edges_in_upperhem[edges_iteration[0]] == 1 && hostSetInfoVecs.boundaries_in_upperhem[edges_iteration[0]] == 0){
                //     linear_0 = (linear_spring_constant/(2.0))*(DISTANCE - generalParams.Rmin_growth)*
                //     (DISTANCE - generalParams.Rmin_growth);
                // }
                // else{
                    linear_0 = (linear_spring_constant/(2.0))*(DISTANCE - generalParams.Rmin)*
                    (DISTANCE - generalParams.Rmin);
                //}
            //}
            
            //else if (DISTANCE < generalParams.Rmin ){
            //    linear_0 = linearSpringInfoVecs.spring_constant_rep1*(1-exp(-linearSpringInfoVecs.spring_constant_rep2*(DISTANCE - generalParams.Rmin)))*(1-exp(-linearSpringInfoVecs.spring_constant_rep2*(DISTANCE - generalParams.Rmin)));
            //}
            
            /*double linear_0 = (linearSpringInfoVecs.spring_constant/2)*(sqrt(
                pow(hostSetInfoVecs.nodeLocX[edge_end] - hostSetInfoVecs.nodeLocX[edge_start], 2.0) + 
                pow(hostSetInfoVecs.nodeLocY[edge_end] - hostSetInfoVecs.nodeLocY[edge_start], 2.0) + 
                pow(hostSetInfoVecs.nodeLocZ[edge_end] - hostSetInfoVecs.nodeLocZ[edge_start], 2.0)) - linearSpringInfoVecs.edge_initial_length[0])*
                (sqrt(
                pow(hostSetInfoVecs.nodeLocX[edge_end] - hostSetInfoVecs.nodeLocX[edge_start], 2.0) + 
                pow(hostSetInfoVecs.nodeLocY[edge_end] - hostSetInfoVecs.nodeLocY[edge_start], 2.0) + 
                pow(hostSetInfoVecs.nodeLocZ[edge_end] - hostSetInfoVecs.nodeLocZ[edge_start], 2.0)) - linearSpringInfoVecs.edge_initial_length[0]);*/
                ////std::cout<<"the energy of this edge is = "<<linear_0<<std::endl;
            
            int H0n1 = edge_start;//hostSetInfoVecs.triangles2Nodes_1[H0];
            int H0n2 = edge_end;//hostSetInfoVecs.triangles2Nodes_2[H0];
            int H0n3 = HEAD;//hostSetInfoVecs.triangles2Nodes_3[H0];
            int T0n1 = edge_start;//hostSetInfoVecs.triangles2Nodes_1[T0];
            int T0n2 = TAIL;//hostSetInfoVecs.triangles2Nodes_2[T0];
            int T0n3 = edge_end;//hostSetInfoVecs.triangles2Nodes_3[T0];
            double a = sqrt(pow((hostSetInfoVecs.nodeLocX[H0n2] - hostSetInfoVecs.nodeLocX[H0n1]),2.0) + 
                        pow((hostSetInfoVecs.nodeLocY[H0n2] - hostSetInfoVecs.nodeLocY[H0n1]),2.0) +
                        pow((hostSetInfoVecs.nodeLocZ[H0n2] - hostSetInfoVecs.nodeLocZ[H0n1]),2.0)
                        );
            double b = sqrt(pow((hostSetInfoVecs.nodeLocX[H0n3] - hostSetInfoVecs.nodeLocX[H0n1]),2.0) + 
                        pow((hostSetInfoVecs.nodeLocY[H0n3] - hostSetInfoVecs.nodeLocY[H0n1]),2.0) +
                        pow((hostSetInfoVecs.nodeLocZ[H0n3] - hostSetInfoVecs.nodeLocZ[H0n1]),2.0)
                        );
            double c = sqrt(pow((hostSetInfoVecs.nodeLocX[H0n3] - hostSetInfoVecs.nodeLocX[H0n2]),2.0) + 
                        pow((hostSetInfoVecs.nodeLocY[H0n3] - hostSetInfoVecs.nodeLocY[H0n2]),2.0) +
                        pow((hostSetInfoVecs.nodeLocZ[H0n3] - hostSetInfoVecs.nodeLocZ[H0n2]),2.0)
                        );
            double mean_abc = (a + b + c)/2;
            double d = sqrt(pow((hostSetInfoVecs.nodeLocX[T0n2] - hostSetInfoVecs.nodeLocX[T0n1]),2.0) + 
                        pow((hostSetInfoVecs.nodeLocY[T0n2] - hostSetInfoVecs.nodeLocY[T0n1]),2.0) +
                        pow((hostSetInfoVecs.nodeLocZ[T0n2] - hostSetInfoVecs.nodeLocZ[T0n1]),2.0)
                        );
            double e = sqrt(pow((hostSetInfoVecs.nodeLocX[T0n3] - hostSetInfoVecs.nodeLocX[T0n1]),2.0) + 
                        pow((hostSetInfoVecs.nodeLocY[T0n3] - hostSetInfoVecs.nodeLocY[T0n1]),2.0) +
                        pow((hostSetInfoVecs.nodeLocZ[T0n3] - hostSetInfoVecs.nodeLocZ[T0n1]),2.0)
                        );
            double f = sqrt(pow((hostSetInfoVecs.nodeLocX[T0n3] - hostSetInfoVecs.nodeLocX[T0n2]),2.0) + 
                        pow((hostSetInfoVecs.nodeLocY[T0n3] - hostSetInfoVecs.nodeLocY[T0n2]),2.0) +
                        pow((hostSetInfoVecs.nodeLocZ[T0n3] - hostSetInfoVecs.nodeLocZ[T0n2]),2.0)
                        );
            double mean_def = (d + e + f)/2.0;
            double area_spring_constant_1, area_spring_constant_2;
            /*if (generalParams.triangles_in_upperhem[H0] == 1){
                area_spring_constant_1 = areaTriangleInfoVecs.spring_constant_weak;
            }
            else if (generalParams.triangles_in_upperhem[H0] == 0){
                area_spring_constant_1 = (areaTriangleInfoVecs.spring_constant_weak + areaTriangleInfoVecs.spring_constant)/2.0;
            }
            else{
                area_spring_constant_1 = areaTriangleInfoVecs.spring_constant;
            }*/
            if (generalParams.SCALE_TYPE == 0){
                area_spring_constant_1 = areaTriangleInfoVecs.spring_constant*((1.0 - ((1.0/sqrt(2*3.14159*generalParams.gausssigma))*exp(-pow(hostSetInfoVecs.scaling_per_edge[iedge],2.0)/generalParams.gausssigma))) +
                                            (1.0 - ((1.0/sqrt(2*3.14159*generalParams.gausssigma))*exp(-pow(hostSetInfoVecs.scaling_per_edge[H1],2.0)/generalParams.gausssigma))) +
                                            (1.0 - ((1.0/sqrt(2*3.14159*generalParams.gausssigma))*exp(-pow(hostSetInfoVecs.scaling_per_edge[H2],2.0)/generalParams.gausssigma))))/3.0;
                                            if (area_spring_constant_1 < areaTriangleInfoVecs.spring_constant_weak){area_spring_constant_1 = areaTriangleInfoVecs.spring_constant_weak;}
            }
            else if(generalParams.SCALE_TYPE == 1){
                area_spring_constant_1 = ((areaTriangleInfoVecs.spring_constant_weak*2.0)*pow(hostSetInfoVecs.scaling_per_edge[iedge],generalParams.scaling_pow) + areaTriangleInfoVecs.spring_constant_weak*(1.0 - pow(hostSetInfoVecs.scaling_per_edge[iedge],generalParams.scaling_pow)) +
                                        (areaTriangleInfoVecs.spring_constant_weak*2.0)*pow(hostSetInfoVecs.scaling_per_edge[H1],generalParams.scaling_pow) + areaTriangleInfoVecs.spring_constant_weak*(1.0 - pow(hostSetInfoVecs.scaling_per_edge[H1],generalParams.scaling_pow)) +
                                        (areaTriangleInfoVecs.spring_constant_weak*2.0)*pow(hostSetInfoVecs.scaling_per_edge[H2],generalParams.scaling_pow) + areaTriangleInfoVecs.spring_constant_weak*(1.0 - pow(hostSetInfoVecs.scaling_per_edge[H2],generalParams.scaling_pow)))/3.0;
            //  area_spring_constant_1 = (areaTriangleInfoVecs.spring_constant*pow(hostSetInfoVecs.scaling_per_edge[iedge],generalParams.scaling_pow) + areaTriangleInfoVecs.spring_constant_weak*(1.0 - pow(hostSetInfoVecs.scaling_per_edge[iedge],generalParams.scaling_pow)) +
            //                             areaTriangleInfoVecs.spring_constant*pow(hostSetInfoVecs.scaling_per_edge[H1],generalParams.scaling_pow) + areaTriangleInfoVecs.spring_constant_weak*(1.0 - pow(hostSetInfoVecs.scaling_per_edge[H1],generalParams.scaling_pow)) +
            //                             areaTriangleInfoVecs.spring_constant*pow(hostSetInfoVecs.scaling_per_edge[H2],generalParams.scaling_pow) + areaTriangleInfoVecs.spring_constant_weak*(1.0 - pow(hostSetInfoVecs.scaling_per_edge[H2],generalParams.scaling_pow)))/3.0;
            }
            else if (generalParams.SCALE_TYPE == 2){
                area_spring_constant_1 = ((areaTriangleInfoVecs.spring_constant - (areaTriangleInfoVecs.spring_constant - areaTriangleInfoVecs.spring_constant_weak)*hostSetInfoVecs.scaling_per_edge[iedge]) +
                                        (areaTriangleInfoVecs.spring_constant - (areaTriangleInfoVecs.spring_constant - areaTriangleInfoVecs.spring_constant_weak)*hostSetInfoVecs.scaling_per_edge[H1]) +
                                        (areaTriangleInfoVecs.spring_constant - (areaTriangleInfoVecs.spring_constant - areaTriangleInfoVecs.spring_constant_weak)*hostSetInfoVecs.scaling_per_edge[H2]))/3.0;
                if (area_spring_constant_1 < areaTriangleInfoVecs.spring_constant_weak){area_spring_constant_1 = areaTriangleInfoVecs.spring_constant_weak;}
            }
            else if (generalParams.SCALE_TYPE == 3){
                if (hostSetInfoVecs.triangles_in_upperhem[H0] == 1){// && hostSetInfoVecs.triangles_in_tip[H0] == 1){
                    area_spring_constant_1 = areaTriangleInfoVecs.spring_constant_weak;
                }
                else if (hostSetInfoVecs.triangles_in_upperhem[H0] == 0){//1 && hostSetInfoVecs.triangles_in_tip[H0] != 1){
                    // area_spring_constant_1 = (areaTriangleInfoVecs.spring_constant_weak*2.0);
                    area_spring_constant_1 = (areaTriangleInfoVecs.spring_constant_weak + areaTriangleInfoVecs.spring_constant)/2.0;
                }
                else{
                    area_spring_constant_1 = areaTriangleInfoVecs.spring_constant;
                }
            }
            else if (generalParams.SCALE_TYPE == 4){
               // double scaling = 0.0;//areaTriangleInfoVecs.spring_constant_weak/areaTriangleInfoVecs.spring_constant;
			    // area_spring_constant_1 = (areaTriangleInfoVecs.spring_constant*((1.0/(1.0+pow(generalParams.hilleqnconst/hostSetInfoVecs.scaling_per_edge[iedge], generalParams.hilleqnpow)))*(1-scaling) + scaling) +
                //                    areaTriangleInfoVecs.spring_constant*((1.0/(1.0+pow(generalParams.hilleqnconst/hostSetInfoVecs.scaling_per_edge[H1], generalParams.hilleqnpow)))*(1-scaling) + scaling) +
                //                    areaTriangleInfoVecs.spring_constant*((1.0/(1.0+pow(generalParams.hilleqnconst/hostSetInfoVecs.scaling_per_edge[H2], generalParams.hilleqnpow)))*(1-scaling) + scaling))/3.0;
			    double spectrum = areaTriangleInfoVecs.spring_constant - areaTriangleInfoVecs.spring_constant_weak;
                area_spring_constant_1 = (areaTriangleInfoVecs.spring_constant_weak + ((1.0/(1.0+pow(generalParams.hilleqnconst/hostSetInfoVecs.scaling_per_edge[iedge], generalParams.hilleqnpow)))*spectrum) +
                   areaTriangleInfoVecs.spring_constant_weak + ((1.0/(1.0+pow(generalParams.hilleqnconst/hostSetInfoVecs.scaling_per_edge[H1], generalParams.hilleqnpow)))*spectrum) +
                   areaTriangleInfoVecs.spring_constant_weak + ((1.0/(1.0+pow(generalParams.hilleqnconst/hostSetInfoVecs.scaling_per_edge[H2], generalParams.hilleqnpow)))*spectrum))/3.0;
                if (area_spring_constant_1 < areaTriangleInfoVecs.spring_constant_weak){area_spring_constant_1 = areaTriangleInfoVecs.spring_constant_weak;}
		       }
            /*if (generalParams.triangles_in_upperhem[T0] == 1){
                area_spring_constant_2 = areaTriangleInfoVecs.spring_constant_weak;
            }
            else if (generalParams.triangles_in_upperhem[T0] == 0){
                area_spring_constant_2 = (areaTriangleInfoVecs.spring_constant_weak + areaTriangleInfoVecs.spring_constant)/2.0;
            }
            else{
                area_spring_constant_2 = areaTriangleInfoVecs.spring_constant;
            }*/

            if (generalParams.SCALE_TYPE == 0){
                area_spring_constant_2 = areaTriangleInfoVecs.spring_constant*((1.0 - ((1.0/sqrt(2*3.14159*generalParams.gausssigma))*exp(-pow(hostSetInfoVecs.scaling_per_edge[iedge],2.0)/generalParams.gausssigma))) +
                                            (1.0 - ((1.0/sqrt(2*3.14159*generalParams.gausssigma))*exp(-pow(hostSetInfoVecs.scaling_per_edge[T1],2.0)/generalParams.gausssigma))) +
                                            (1.0 - ((1.0/sqrt(2*3.14159*generalParams.gausssigma))*exp(-pow(hostSetInfoVecs.scaling_per_edge[T2],2.0)/generalParams.gausssigma))))/3.0;
                                            if (area_spring_constant_2 < areaTriangleInfoVecs.spring_constant_weak){area_spring_constant_2 = areaTriangleInfoVecs.spring_constant_weak;}
            }
            else if (generalParams.SCALE_TYPE == 1){
                area_spring_constant_2 = ((areaTriangleInfoVecs.spring_constant_weak*2.0)*pow(hostSetInfoVecs.scaling_per_edge[iedge],generalParams.scaling_pow) + areaTriangleInfoVecs.spring_constant_weak*(1.0 - pow(hostSetInfoVecs.scaling_per_edge[iedge],generalParams.scaling_pow)) +
                                        (areaTriangleInfoVecs.spring_constant_weak*2.0)*pow(hostSetInfoVecs.scaling_per_edge[T1],generalParams.scaling_pow) + areaTriangleInfoVecs.spring_constant_weak*(1.0 - pow(hostSetInfoVecs.scaling_per_edge[T1],generalParams.scaling_pow)) +
                                        (areaTriangleInfoVecs.spring_constant_weak*2.0)*pow(hostSetInfoVecs.scaling_per_edge[T2],generalParams.scaling_pow) + areaTriangleInfoVecs.spring_constant_weak*(1.0 - pow(hostSetInfoVecs.scaling_per_edge[T2],generalParams.scaling_pow)))/3.0;
                // area_spring_constant_2 = (areaTriangleInfoVecs.spring_constant*pow(hostSetInfoVecs.scaling_per_edge[iedge],generalParams.scaling_pow) + areaTriangleInfoVecs.spring_constant_weak*(1.0 - pow(hostSetInfoVecs.scaling_per_edge[iedge],generalParams.scaling_pow)) +
                //                         areaTriangleInfoVecs.spring_constant*pow(hostSetInfoVecs.scaling_per_edge[T1],generalParams.scaling_pow) + areaTriangleInfoVecs.spring_constant_weak*(1.0 - pow(hostSetInfoVecs.scaling_per_edge[T1],generalParams.scaling_pow)) +
                //                         areaTriangleInfoVecs.spring_constant*pow(hostSetInfoVecs.scaling_per_edge[T2],generalParams.scaling_pow) + areaTriangleInfoVecs.spring_constant_weak*(1.0 - pow(hostSetInfoVecs.scaling_per_edge[T2],generalParams.scaling_pow)))/3.0;
            }
            else if (generalParams.SCALE_TYPE == 2){
                area_spring_constant_2 = ((areaTriangleInfoVecs.spring_constant - (areaTriangleInfoVecs.spring_constant - areaTriangleInfoVecs.spring_constant_weak)*hostSetInfoVecs.scaling_per_edge[iedge]) +
                                        (areaTriangleInfoVecs.spring_constant - (areaTriangleInfoVecs.spring_constant - areaTriangleInfoVecs.spring_constant_weak)*hostSetInfoVecs.scaling_per_edge[T1])+
                                        (areaTriangleInfoVecs.spring_constant - (areaTriangleInfoVecs.spring_constant - areaTriangleInfoVecs.spring_constant_weak)*hostSetInfoVecs.scaling_per_edge[T2]))/3.0;
                if (area_spring_constant_2 < areaTriangleInfoVecs.spring_constant_weak){area_spring_constant_2 = areaTriangleInfoVecs.spring_constant_weak;}
            }
            else if (generalParams.SCALE_TYPE == 3){
                if (hostSetInfoVecs.triangles_in_upperhem[T0] == 1){// && hostSetInfoVecs.triangles_in_tip[T0] == 1){
                area_spring_constant_2 = areaTriangleInfoVecs.spring_constant_weak;
                }
                else if (hostSetInfoVecs.triangles_in_upperhem[T0] == 0){//1 && hostSetInfoVecs.triangles_in_tip[T0] != 1){
                    // area_spring_constant_2 = (areaTriangleInfoVecs.spring_constant_weak*2.0);
                    area_spring_constant_2 = (areaTriangleInfoVecs.spring_constant_weak + areaTriangleInfoVecs.spring_constant)/2.0;
                }
                else{
                    area_spring_constant_2 = areaTriangleInfoVecs.spring_constant;
                }
            }
            else if (generalParams.SCALE_TYPE == 4){
                //double scaling = 0.0;//areaTriangleInfoVecs.spring_constant_weak/areaTriangleInfoVecs.spring_constant;
			    //area_spring_constant_2 = (areaTriangleInfoVecs.spring_constant*((1.0/(1.0+pow(generalParams.hilleqnconst/hostSetInfoVecs.scaling_per_edge[iedge], generalParams.hilleqnpow)))*(1-scaling) + scaling) +
                //                          areaTriangleInfoVecs.spring_constant*((1.0/(1.0+pow(generalParams.hilleqnconst/hostSetInfoVecs.scaling_per_edge[T1], generalParams.hilleqnpow)))*(1-scaling) + scaling) +
                //                          areaTriangleInfoVecs.spring_constant*((1.0/(1.0+pow(generalParams.hilleqnconst/hostSetInfoVecs.scaling_per_edge[T2], generalParams.hilleqnpow)))*(1-scaling) + scaling))/3.0;
			    double spectrum = areaTriangleInfoVecs.spring_constant - areaTriangleInfoVecs.spring_constant_weak;
                area_spring_constant_2 = (areaTriangleInfoVecs.spring_constant_weak + ((1.0/(1.0+pow(generalParams.hilleqnconst/hostSetInfoVecs.scaling_per_edge[iedge], generalParams.hilleqnpow)))*spectrum) +
                                          areaTriangleInfoVecs.spring_constant_weak + ((1.0/(1.0+pow(generalParams.hilleqnconst/hostSetInfoVecs.scaling_per_edge[T1], generalParams.hilleqnpow)))*spectrum) +
                                          areaTriangleInfoVecs.spring_constant_weak + ((1.0/(1.0+pow(generalParams.hilleqnconst/hostSetInfoVecs.scaling_per_edge[T2], generalParams.hilleqnpow)))*spectrum))/3.0;
                if (area_spring_constant_2 < areaTriangleInfoVecs.spring_constant_weak){area_spring_constant_2 = areaTriangleInfoVecs.spring_constant_weak;}
		   }
            double area_H0 = sqrt(mean_abc*(mean_abc - a)*(mean_abc - b)*(mean_abc - c));
            double area_T0 = sqrt(mean_def*(mean_def - d)*(mean_def - e)*(mean_def - f));
            double area_0_energy = area_spring_constant_1*pow((area_H0 - areaTriangleInfoVecs.initial_area),2.0)/(2*areaTriangleInfoVecs.initial_area) +
                                area_spring_constant_2*pow((area_T0 - areaTriangleInfoVecs.initial_area),2.0)/(2*areaTriangleInfoVecs.initial_area);
            
            //double vol_H0 = (1.0/3.0)*(P0x_vol1*N1x_vol + P0y_vol1*N1y_vol + P0z_vol1*N1z_vol)*area_H0;
            //double vol_T0 = (1.0/3.0)*(P0x_vol2*N2x_vol + P0y_vol2*N2y_vol + P0z_vol2*N2z_vol)*area_T0;
            //vol_0 = vol_H0 + vol_T0;
            double E_0 = linear_0 + bend_0 + area_0_energy + rep_0;// + generalParams.volume_energy;
            // std::cout<<"old linear energy: "<<linear_0<<" , old length = "<<DISTANCE<<std::endl;
            // std::cout<<"old bend energy: "<<bend_0<<std::endl;
            // std::cout<<"old area energy: "<<area_0_energy<<std::endl;
            // std::cout<<"old total energy: "<<E_0<<std::endl;

            
            //Flip the edge, build the data structure for the smaller system.
            /*bool BAD_CHOICE = false;
            int temp_edges2Nodes_1 = TAIL;
            int temp_edges2Nodes_2 = HEAD;

            int temp_nndata_HEAD = nndata[HEAD] + 1;
            int temp_nndata_TAIL = nndata[TAIL] + 1;
            int temp_nndata_edge_start = nndata[edge_start] - 1;
            
            int temp_nndata_edge_end = nndata[edge_end] - 1;
            
            if (boundary_node[HEAD] == false && temp_nndata_HEAD < 3){
                BAD_CHOICE = true;
            }
            else if (boundary_node[TAIL] == false && temp_nndata_TAIL < 3){
                BAD_CHOICE = true;
            }
            else if (boundary_node[edge_start] == false && temp_nndata_edge_start < 3){
                BAD_CHOICE = true;
            }
            else if (boundary_node[edge_end] == false && temp_nndata_edge_end < 3){
                BAD_CHOICE = true;
            }
            else if (boundary_node[HEAD] == false && temp_nndata_HEAD > 12){
                BAD_CHOICE = true;
            }
            else if (boundary_node[TAIL] == false && temp_nndata_TAIL > 12){
                BAD_CHOICE = true;
            }
            else if (boundary_node[edge_start] == false && temp_nndata_edge_start > 12){
                BAD_CHOICE = true;
            }
            else if (boundary_node[edge_end] == false && temp_nndata_edge_end > 12){
                BAD_CHOICE = true;
            }
            else {
                BAD_CHOICE = false;
            }*/


            if (BAD_CHOICE == false) {
                /*temp_edges2Nodes_1[iedge] = TAIL;
                temp_edges2Nodes_2[iedge] = HEAD;
                temp_nndata[HEAD] = temp_nndata[HEAD] + 1;
                temp_nndata[TAIL] = temp_nndata[TAIL] + 1;
                temp_nndata[edge_start] = temp_nndata[edge_start] - 1;
                temp_nndata[edge_end] = temp_nndata[edge_end] - 1;*/

                //The labeling of neighboring edge will as follows after swap:
                //          
                //           edge_start
                //           *        *
                //         T1    H0     H1
                //        *              *
                //      TAIL ----------> HEAD
                //        *              *
                //         T2    T0    H2
                //           *        *
                //            edge_end
                //
                //Now we will update the temporary data structure to accomodate the edgeswap
                
                //Update the new triangles2Nodes information
                /*temp_triangles2Nodes_1[H0] = HEAD;
                temp_triangles2Nodes_2[H0] = edge_start;
                temp_triangles2Nodes_3[H0] = TAIL;
                temp_triangles2Nodes_1[T0] = HEAD;
                temp_triangles2Nodes_2[T0] = TAIL;
                temp_triangles2Nodes_3[T0] = edge_end;*/

                
                //Creating vectors to compute the normal vectors under the swapped configuration.
                int H1t1 = hostSetInfoVecs.edges2Triangles_1[H1];
                //std::cout<<"H1t1 = "<<H1t1<<std::endl;
                int H1t2 = hostSetInfoVecs.edges2Triangles_2[H1]; 
                //std::cout<<"H1t2 = "<<H1t2<<std::endl;//These are the associated triangles to edge H1
                //For the following if statement, we identify the triangles that are affected by the edge-swap.
                //Since we do not know the exact index of the affected triangle, we use the if statement to consider possible cases.
                //This gives us the vectors necessary to compute unit normal vectors required for bending energy.
                

        
                double H1t1_vec1x;
                double H1t1_vec1y;
                double H1t1_vec1z;
                double H1t1_vec2x;
                double H1t1_vec2y;
                double H1t1_vec2z;
                double H1t2_vec1x;
                double H1t2_vec1y;
                double H1t2_vec1z;
                double H1t2_vec2x;
                double H1t2_vec2y;
                double H1t2_vec2z;

                double H2t1_vec1x;
                double H2t1_vec1y;
                double H2t1_vec1z;
                double H2t1_vec2x;
                double H2t1_vec2y;
                double H2t1_vec2z;
                double H2t2_vec1x;
                double H2t2_vec1y;
                double H2t2_vec1z;
                double H2t2_vec2x;
                double H2t2_vec2y;
                double H2t2_vec2z;

                double T1t2_vec1x;
                double T1t2_vec1y;
                double T1t2_vec1z;
                double T1t2_vec2x;
                double T1t2_vec2y;
                double T1t2_vec2z;
                double T1t1_vec1x;
                double T1t1_vec1y;
                double T1t1_vec1z;
                double T1t1_vec2x;
                double T1t1_vec2y;
                double T1t1_vec2z;

                double T2t2_vec1x;
                double T2t2_vec1y;
                double T2t2_vec1z;
                double T2t2_vec2x;
                double T2t2_vec2y;
                double T2t2_vec2z;
                double T2t1_vec1x;
                double T2t1_vec1y;
                double T2t1_vec1z;
                double T2t1_vec2x;
                double T2t1_vec2y;
                double T2t1_vec2z;
                
                //           edge_start
                //           *        *
                //         T1    H0     H1
                //        *              *
                //      TAIL ----------> HEAD
                //        *              *
                //         T2    T0    H2
                //           *        *
                //            edge_end
                
                if (H1t1 == H0){H1t1_vec1x = hostSetInfoVecs.nodeLocX[edge_start] - hostSetInfoVecs.nodeLocX[HEAD];
                                H1t1_vec1y = hostSetInfoVecs.nodeLocY[edge_start] - hostSetInfoVecs.nodeLocY[HEAD];
                                H1t1_vec1z = hostSetInfoVecs.nodeLocZ[edge_start] - hostSetInfoVecs.nodeLocZ[HEAD];
                                H1t1_vec2x = hostSetInfoVecs.nodeLocX[TAIL] - hostSetInfoVecs.nodeLocX[HEAD];
                                H1t1_vec2y = hostSetInfoVecs.nodeLocY[TAIL] - hostSetInfoVecs.nodeLocY[HEAD];
                                H1t1_vec2z = hostSetInfoVecs.nodeLocZ[TAIL] - hostSetInfoVecs.nodeLocZ[HEAD];
                                H1t2_vec1x = hostSetInfoVecs.nodeLocX[hostSetInfoVecs.triangles2Nodes_2[H1t2]] - hostSetInfoVecs.nodeLocX[hostSetInfoVecs.triangles2Nodes_1[H1t2]];
                                H1t2_vec1y = hostSetInfoVecs.nodeLocY[hostSetInfoVecs.triangles2Nodes_2[H1t2]] - hostSetInfoVecs.nodeLocY[hostSetInfoVecs.triangles2Nodes_1[H1t2]];
                                H1t2_vec1z = hostSetInfoVecs.nodeLocZ[hostSetInfoVecs.triangles2Nodes_2[H1t2]] - hostSetInfoVecs.nodeLocZ[hostSetInfoVecs.triangles2Nodes_1[H1t2]];
                                H1t2_vec2x = hostSetInfoVecs.nodeLocX[hostSetInfoVecs.triangles2Nodes_3[H1t2]] - hostSetInfoVecs.nodeLocX[hostSetInfoVecs.triangles2Nodes_1[H1t2]];
                                H1t2_vec2y = hostSetInfoVecs.nodeLocY[hostSetInfoVecs.triangles2Nodes_3[H1t2]] - hostSetInfoVecs.nodeLocY[hostSetInfoVecs.triangles2Nodes_1[H1t2]];
                                H1t2_vec2z = hostSetInfoVecs.nodeLocZ[hostSetInfoVecs.triangles2Nodes_3[H1t2]] - hostSetInfoVecs.nodeLocZ[hostSetInfoVecs.triangles2Nodes_1[H1t2]];
                                //std::cout<<"H1t1 = H0"<<std::endl;
                                }
                else if (H1t2 == H0){H1t2_vec1x = hostSetInfoVecs.nodeLocX[edge_start] - hostSetInfoVecs.nodeLocX[HEAD];
                                H1t2_vec1y = hostSetInfoVecs.nodeLocY[edge_start] - hostSetInfoVecs.nodeLocY[HEAD];
                                H1t2_vec1z = hostSetInfoVecs.nodeLocZ[edge_start] - hostSetInfoVecs.nodeLocZ[HEAD];
                                H1t2_vec2x = hostSetInfoVecs.nodeLocX[TAIL] - hostSetInfoVecs.nodeLocX[HEAD];
                                H1t2_vec2y = hostSetInfoVecs.nodeLocY[TAIL] - hostSetInfoVecs.nodeLocY[HEAD];
                                H1t2_vec2z = hostSetInfoVecs.nodeLocZ[TAIL] - hostSetInfoVecs.nodeLocZ[HEAD];
                                H1t1_vec1x = hostSetInfoVecs.nodeLocX[hostSetInfoVecs.triangles2Nodes_2[H1t1]] - hostSetInfoVecs.nodeLocX[hostSetInfoVecs.triangles2Nodes_1[H1t1]];
                                H1t1_vec1y = hostSetInfoVecs.nodeLocY[hostSetInfoVecs.triangles2Nodes_2[H1t1]] - hostSetInfoVecs.nodeLocY[hostSetInfoVecs.triangles2Nodes_1[H1t1]];
                                H1t1_vec1z = hostSetInfoVecs.nodeLocZ[hostSetInfoVecs.triangles2Nodes_2[H1t1]] - hostSetInfoVecs.nodeLocZ[hostSetInfoVecs.triangles2Nodes_1[H1t1]];
                                H1t1_vec2x = hostSetInfoVecs.nodeLocX[hostSetInfoVecs.triangles2Nodes_3[H1t1]] - hostSetInfoVecs.nodeLocX[hostSetInfoVecs.triangles2Nodes_1[H1t1]];
                                H1t1_vec2y = hostSetInfoVecs.nodeLocY[hostSetInfoVecs.triangles2Nodes_3[H1t1]] - hostSetInfoVecs.nodeLocY[hostSetInfoVecs.triangles2Nodes_1[H1t1]];
                                H1t1_vec2z = hostSetInfoVecs.nodeLocZ[hostSetInfoVecs.triangles2Nodes_3[H1t1]] - hostSetInfoVecs.nodeLocZ[hostSetInfoVecs.triangles2Nodes_1[H1t1]];
                                //std::cout<<"H1t2 = H0"<<std::endl;
                                }
                int H2t1 = hostSetInfoVecs.edges2Triangles_1[H2];
                int H2t2 = hostSetInfoVecs.edges2Triangles_2[H2]; //These are the associated triangles to edge H2        
                if (H2t1 == H0){//In this case H2t1 turns into T0.
                                H2t1_vec1x = hostSetInfoVecs.nodeLocX[TAIL] - hostSetInfoVecs.nodeLocX[HEAD];
                                H2t1_vec1y = hostSetInfoVecs.nodeLocY[TAIL] - hostSetInfoVecs.nodeLocY[HEAD];
                                H2t1_vec1z = hostSetInfoVecs.nodeLocZ[TAIL] - hostSetInfoVecs.nodeLocZ[HEAD];
                                H2t1_vec2x = hostSetInfoVecs.nodeLocX[edge_end] - hostSetInfoVecs.nodeLocX[HEAD];
                                H2t1_vec2y = hostSetInfoVecs.nodeLocY[edge_end] - hostSetInfoVecs.nodeLocY[HEAD];
                                H2t1_vec2z = hostSetInfoVecs.nodeLocZ[edge_end] - hostSetInfoVecs.nodeLocZ[HEAD];
                                H2t2_vec1x = hostSetInfoVecs.nodeLocX[hostSetInfoVecs.triangles2Nodes_2[H2t2]] - hostSetInfoVecs.nodeLocX[hostSetInfoVecs.triangles2Nodes_1[H2t2]];
                                H2t2_vec1y = hostSetInfoVecs.nodeLocY[hostSetInfoVecs.triangles2Nodes_2[H2t2]] - hostSetInfoVecs.nodeLocY[hostSetInfoVecs.triangles2Nodes_1[H2t2]];
                                H2t2_vec1z = hostSetInfoVecs.nodeLocZ[hostSetInfoVecs.triangles2Nodes_2[H2t2]] - hostSetInfoVecs.nodeLocZ[hostSetInfoVecs.triangles2Nodes_1[H2t2]];
                                H2t2_vec2x = hostSetInfoVecs.nodeLocX[hostSetInfoVecs.triangles2Nodes_3[H2t2]] - hostSetInfoVecs.nodeLocX[hostSetInfoVecs.triangles2Nodes_1[H2t2]];
                                H2t2_vec2y = hostSetInfoVecs.nodeLocY[hostSetInfoVecs.triangles2Nodes_3[H2t2]] - hostSetInfoVecs.nodeLocY[hostSetInfoVecs.triangles2Nodes_1[H2t2]];
                                H2t2_vec2z = hostSetInfoVecs.nodeLocZ[hostSetInfoVecs.triangles2Nodes_3[H2t2]] - hostSetInfoVecs.nodeLocZ[hostSetInfoVecs.triangles2Nodes_1[H2t2]];
                                //std::cout<<"H2t1 = H0"<<std::endl;
                                }
                else if (H2t2 == H0){//In this case H2t2 tunrs into T0
                                H2t2_vec1x = hostSetInfoVecs.nodeLocX[TAIL] - hostSetInfoVecs.nodeLocX[HEAD];
                                H2t2_vec1y = hostSetInfoVecs.nodeLocY[TAIL] - hostSetInfoVecs.nodeLocY[HEAD];
                                H2t2_vec1z = hostSetInfoVecs.nodeLocZ[TAIL] - hostSetInfoVecs.nodeLocZ[HEAD];
                                H2t2_vec2x = hostSetInfoVecs.nodeLocX[edge_end] - hostSetInfoVecs.nodeLocX[HEAD];
                                H2t2_vec2y = hostSetInfoVecs.nodeLocY[edge_end] - hostSetInfoVecs.nodeLocY[HEAD];
                                H2t2_vec2z = hostSetInfoVecs.nodeLocZ[edge_end] - hostSetInfoVecs.nodeLocZ[HEAD];
                                H2t1_vec1x = hostSetInfoVecs.nodeLocX[hostSetInfoVecs.triangles2Nodes_2[H2t1]] - hostSetInfoVecs.nodeLocX[hostSetInfoVecs.triangles2Nodes_1[H2t1]];
                                H2t1_vec1y = hostSetInfoVecs.nodeLocY[hostSetInfoVecs.triangles2Nodes_2[H2t1]] - hostSetInfoVecs.nodeLocY[hostSetInfoVecs.triangles2Nodes_1[H2t1]];
                                H2t1_vec1z = hostSetInfoVecs.nodeLocZ[hostSetInfoVecs.triangles2Nodes_2[H2t1]] - hostSetInfoVecs.nodeLocZ[hostSetInfoVecs.triangles2Nodes_1[H2t1]];
                                H2t1_vec2x = hostSetInfoVecs.nodeLocX[hostSetInfoVecs.triangles2Nodes_3[H2t1]] - hostSetInfoVecs.nodeLocX[hostSetInfoVecs.triangles2Nodes_1[H2t1]];
                                H2t1_vec2y = hostSetInfoVecs.nodeLocY[hostSetInfoVecs.triangles2Nodes_3[H2t1]] - hostSetInfoVecs.nodeLocY[hostSetInfoVecs.triangles2Nodes_1[H2t1]];
                                H2t1_vec2z = hostSetInfoVecs.nodeLocZ[hostSetInfoVecs.triangles2Nodes_3[H2t1]] - hostSetInfoVecs.nodeLocZ[hostSetInfoVecs.triangles2Nodes_1[H2t1]];
                                //std::cout<<"H2t2 = H0"<<std::endl;
                                }
                int T1t1 = hostSetInfoVecs.edges2Triangles_1[T1];
                //std::cout<<"T1t1 = "<<T1t1<<std::endl;
                int T1t2 = hostSetInfoVecs.edges2Triangles_2[T1];
                //std::cout<<"T1t2 = "<<T1t2<<std::endl;
                if (T1t1 == T0){//In this case T1t1 turns into H0.
                                T1t1_vec1x = hostSetInfoVecs.nodeLocX[HEAD] - hostSetInfoVecs.nodeLocX[TAIL];
                                T1t1_vec1y = hostSetInfoVecs.nodeLocY[HEAD] - hostSetInfoVecs.nodeLocY[TAIL];
                                T1t1_vec1z = hostSetInfoVecs.nodeLocZ[HEAD] - hostSetInfoVecs.nodeLocZ[TAIL];
                                T1t1_vec2x = hostSetInfoVecs.nodeLocX[edge_start] - hostSetInfoVecs.nodeLocX[TAIL];
                                T1t1_vec2y = hostSetInfoVecs.nodeLocY[edge_start] - hostSetInfoVecs.nodeLocY[TAIL];
                                T1t1_vec2z = hostSetInfoVecs.nodeLocZ[edge_start] - hostSetInfoVecs.nodeLocZ[TAIL];
                                T1t2_vec1x = hostSetInfoVecs.nodeLocX[hostSetInfoVecs.triangles2Nodes_2[T1t2]] - hostSetInfoVecs.nodeLocX[hostSetInfoVecs.triangles2Nodes_1[T1t2]];
                                T1t2_vec1y = hostSetInfoVecs.nodeLocY[hostSetInfoVecs.triangles2Nodes_2[T1t2]] - hostSetInfoVecs.nodeLocY[hostSetInfoVecs.triangles2Nodes_1[T1t2]];
                                T1t2_vec1z = hostSetInfoVecs.nodeLocZ[hostSetInfoVecs.triangles2Nodes_2[T1t2]] - hostSetInfoVecs.nodeLocZ[hostSetInfoVecs.triangles2Nodes_1[T1t2]];
                                T1t2_vec2x = hostSetInfoVecs.nodeLocX[hostSetInfoVecs.triangles2Nodes_3[T1t2]] - hostSetInfoVecs.nodeLocX[hostSetInfoVecs.triangles2Nodes_1[T1t2]];
                                T1t2_vec2y = hostSetInfoVecs.nodeLocY[hostSetInfoVecs.triangles2Nodes_3[T1t2]] - hostSetInfoVecs.nodeLocY[hostSetInfoVecs.triangles2Nodes_1[T1t2]];
                                T1t2_vec2z = hostSetInfoVecs.nodeLocZ[hostSetInfoVecs.triangles2Nodes_3[T1t2]] - hostSetInfoVecs.nodeLocZ[hostSetInfoVecs.triangles2Nodes_1[T1t2]];
                                //std::cout<<"T1t1 = T0"<<std::endl;
                                }
                else if (T1t2 == T0){//In this case T1t2 turns into H0.
                                T1t2_vec1x = hostSetInfoVecs.nodeLocX[HEAD] - hostSetInfoVecs.nodeLocX[TAIL];
                                T1t2_vec1y = hostSetInfoVecs.nodeLocY[HEAD] - hostSetInfoVecs.nodeLocY[TAIL];
                                T1t2_vec1z = hostSetInfoVecs.nodeLocZ[HEAD] - hostSetInfoVecs.nodeLocZ[TAIL];
                                T1t2_vec2x = hostSetInfoVecs.nodeLocX[edge_start] - hostSetInfoVecs.nodeLocX[TAIL];
                                T1t2_vec2y = hostSetInfoVecs.nodeLocY[edge_start] - hostSetInfoVecs.nodeLocY[TAIL];
                                T1t2_vec2z = hostSetInfoVecs.nodeLocZ[edge_start] - hostSetInfoVecs.nodeLocZ[TAIL];
                                T1t1_vec1x = hostSetInfoVecs.nodeLocX[hostSetInfoVecs.triangles2Nodes_2[T1t1]] - hostSetInfoVecs.nodeLocX[hostSetInfoVecs.triangles2Nodes_1[T1t1]];
                                T1t1_vec1y = hostSetInfoVecs.nodeLocY[hostSetInfoVecs.triangles2Nodes_2[T1t1]] - hostSetInfoVecs.nodeLocY[hostSetInfoVecs.triangles2Nodes_1[T1t1]];
                                T1t1_vec1z = hostSetInfoVecs.nodeLocZ[hostSetInfoVecs.triangles2Nodes_2[T1t1]] - hostSetInfoVecs.nodeLocZ[hostSetInfoVecs.triangles2Nodes_1[T1t1]];
                                T1t1_vec2x = hostSetInfoVecs.nodeLocX[hostSetInfoVecs.triangles2Nodes_3[T1t1]] - hostSetInfoVecs.nodeLocX[hostSetInfoVecs.triangles2Nodes_1[T1t1]];
                                T1t1_vec2y = hostSetInfoVecs.nodeLocY[hostSetInfoVecs.triangles2Nodes_3[T1t1]] - hostSetInfoVecs.nodeLocY[hostSetInfoVecs.triangles2Nodes_1[T1t1]];
                                T1t1_vec2z = hostSetInfoVecs.nodeLocZ[hostSetInfoVecs.triangles2Nodes_3[T1t1]] - hostSetInfoVecs.nodeLocZ[hostSetInfoVecs.triangles2Nodes_1[T1t1]];
                                //std::cout<<"T1t2 = T0"<<std::endl;
                                }
                int T2t1 = hostSetInfoVecs.edges2Triangles_1[T2];
                int T2t2 = hostSetInfoVecs.edges2Triangles_2[T2];
                if (T2t1 == T0){T2t1_vec1x = hostSetInfoVecs.nodeLocX[edge_end] - hostSetInfoVecs.nodeLocX[TAIL];
                                T2t1_vec1y = hostSetInfoVecs.nodeLocY[edge_end] - hostSetInfoVecs.nodeLocY[TAIL];
                                T2t1_vec1z = hostSetInfoVecs.nodeLocZ[edge_end] - hostSetInfoVecs.nodeLocZ[TAIL];
                                T2t1_vec2x = hostSetInfoVecs.nodeLocX[HEAD] - hostSetInfoVecs.nodeLocX[TAIL];
                                T2t1_vec2y = hostSetInfoVecs.nodeLocY[HEAD] - hostSetInfoVecs.nodeLocY[TAIL];
                                T2t1_vec2z = hostSetInfoVecs.nodeLocZ[HEAD] - hostSetInfoVecs.nodeLocZ[TAIL];
                                T2t2_vec1x = hostSetInfoVecs.nodeLocX[hostSetInfoVecs.triangles2Nodes_2[T2t2]] - hostSetInfoVecs.nodeLocX[hostSetInfoVecs.triangles2Nodes_1[T2t2]];
                                T2t2_vec1y = hostSetInfoVecs.nodeLocY[hostSetInfoVecs.triangles2Nodes_2[T2t2]] - hostSetInfoVecs.nodeLocY[hostSetInfoVecs.triangles2Nodes_1[T2t2]];
                                T2t2_vec1z = hostSetInfoVecs.nodeLocZ[hostSetInfoVecs.triangles2Nodes_2[T2t2]] - hostSetInfoVecs.nodeLocZ[hostSetInfoVecs.triangles2Nodes_1[T2t2]];
                                T2t2_vec2x = hostSetInfoVecs.nodeLocX[hostSetInfoVecs.triangles2Nodes_3[T2t2]] - hostSetInfoVecs.nodeLocX[hostSetInfoVecs.triangles2Nodes_1[T2t2]];
                                T2t2_vec2y = hostSetInfoVecs.nodeLocY[hostSetInfoVecs.triangles2Nodes_3[T2t2]] - hostSetInfoVecs.nodeLocY[hostSetInfoVecs.triangles2Nodes_1[T2t2]];
                                T2t2_vec2z = hostSetInfoVecs.nodeLocZ[hostSetInfoVecs.triangles2Nodes_3[T2t2]] - hostSetInfoVecs.nodeLocZ[hostSetInfoVecs.triangles2Nodes_1[T2t2]];
                                //std::cout<<"T2t1 = T0"<<std::endl;
                                }
                else if (T2t2 == T0){
                                T2t2_vec1x = hostSetInfoVecs.nodeLocX[edge_end] - hostSetInfoVecs.nodeLocX[TAIL];
                                T2t2_vec1y = hostSetInfoVecs.nodeLocY[edge_end] - hostSetInfoVecs.nodeLocY[TAIL];
                                T2t2_vec1z = hostSetInfoVecs.nodeLocZ[edge_end] - hostSetInfoVecs.nodeLocZ[TAIL];
                                T2t2_vec2x = hostSetInfoVecs.nodeLocX[HEAD] - hostSetInfoVecs.nodeLocX[TAIL];
                                T2t2_vec2y = hostSetInfoVecs.nodeLocY[HEAD] - hostSetInfoVecs.nodeLocY[TAIL];
                                T2t2_vec2z = hostSetInfoVecs.nodeLocZ[HEAD] - hostSetInfoVecs.nodeLocZ[TAIL];
                                T2t1_vec1x = hostSetInfoVecs.nodeLocX[hostSetInfoVecs.triangles2Nodes_2[T2t1]] - hostSetInfoVecs.nodeLocX[hostSetInfoVecs.triangles2Nodes_1[T2t1]];
                                T2t1_vec1y = hostSetInfoVecs.nodeLocY[hostSetInfoVecs.triangles2Nodes_2[T2t1]] - hostSetInfoVecs.nodeLocY[hostSetInfoVecs.triangles2Nodes_1[T2t1]];
                                T2t1_vec1z = hostSetInfoVecs.nodeLocZ[hostSetInfoVecs.triangles2Nodes_2[T2t1]] - hostSetInfoVecs.nodeLocZ[hostSetInfoVecs.triangles2Nodes_1[T2t1]];
                                T2t1_vec2x = hostSetInfoVecs.nodeLocX[hostSetInfoVecs.triangles2Nodes_3[T2t1]] - hostSetInfoVecs.nodeLocX[hostSetInfoVecs.triangles2Nodes_1[T2t1]];
                                T2t1_vec2y = hostSetInfoVecs.nodeLocY[hostSetInfoVecs.triangles2Nodes_3[T2t1]] - hostSetInfoVecs.nodeLocY[hostSetInfoVecs.triangles2Nodes_1[T2t1]];
                                T2t1_vec2z = hostSetInfoVecs.nodeLocZ[hostSetInfoVecs.triangles2Nodes_3[T2t1]] - hostSetInfoVecs.nodeLocZ[hostSetInfoVecs.triangles2Nodes_1[T2t1]];
                                //std::cout<<"T2t2 = T0"<<std::endl;
                                }
                
                //First calculate the linear spring energy due to edge-swap.
            double linear_1;
            double DISTANCE = sqrt(
            pow(hostSetInfoVecs.nodeLocX[HEAD] - hostSetInfoVecs.nodeLocX[TAIL], 2.0) + 
            pow(hostSetInfoVecs.nodeLocY[HEAD] - hostSetInfoVecs.nodeLocY[TAIL], 2.0) + 
            pow(hostSetInfoVecs.nodeLocZ[HEAD] - hostSetInfoVecs.nodeLocZ[TAIL], 2.0));
            


            /*if (DISTANCE < generalParams.abs_Rmin){
                linear_1 = (linearSpringInfoVecs.spring_constant/2)*(DISTANCE - generalParams.Rmin)*
                (DISTANCE - generalParams.Rmin) + 
                //(linearSpringInfoVecs.spring_constant_rep1/2)*(DISTANCE - generalParams.abs_Rmin)*
                //(DISTANCE - generalParams.abs_Rmin);
                linearSpringInfoVecs.spring_constant_rep1*(1-exp(-linearSpringInfoVecs.spring_constant_rep2*(DISTANCE-generalParams.abs_Rmin)))*
                (1-exp(-linearSpringInfoVecs.spring_constant_rep2*(DISTANCE-generalParams.abs_Rmin)));
            }
            else if (DISTANCE != generalParams.Rmin){
                linear_1 = (linearSpringInfoVecs.spring_constant/2)*(DISTANCE - generalParams.Rmin)*
                (DISTANCE - generalParams.Rmin);
            }*/
            //else{
                // linear_1 = (linearSpringInfoVecs.spring_constant/(2.0*generalParams.Rmin*generalParams.Rmin))*(DISTANCE - generalParams.Rmin)*
                // (DISTANCE - generalParams.Rmin);
                // if (hostSetInfoVecs.edges_in_upperhem[edges_iteration[0]] == 1 && hostSetInfoVecs.boundaries_in_upperhem[edges_iteration[0]] == 0){
                //     linear_1 = (linearSpringInfoVecs.spring_constant/(2.0))*(DISTANCE - generalParams.Rmin_growth)*
                //     (DISTANCE - generalParams.Rmin_growth);
                // }
                //else{
                    linear_1 = (linearSpringInfoVecs.spring_constant/(2.0))*(DISTANCE - generalParams.Rmin)*
                    (DISTANCE - generalParams.Rmin);
                //}
            //}
            
            //else if (DISTANCE < generalParams.Rmin){
            //   linear_1 =   linearSpringInfoVecs.spring_constant_rep1*(1-exp(-linearSpringInfoVecs.spring_constant_rep2*(DISTANCE - generalParams.Rmin)))*(1-exp(-linearSpringInfoVecs.spring_constant_rep2*(DISTANCE - generalParams.Rmin)));
            //}

            double rep_1;
            rep_energy_TAIL_HEAD = 0.0;
            R_TAIL_HEAD = sqrt((hostSetInfoVecs.nodeLocX[edge_start] - hostSetInfoVecs.nodeLocX[edge_end])*(hostSetInfoVecs.nodeLocX[edge_start] - hostSetInfoVecs.nodeLocX[edge_end]) +
                                    (hostSetInfoVecs.nodeLocY[edge_start] - hostSetInfoVecs.nodeLocY[edge_end])*(hostSetInfoVecs.nodeLocY[edge_start] - hostSetInfoVecs.nodeLocY[edge_end]) +
                                    (hostSetInfoVecs.nodeLocZ[edge_start] - hostSetInfoVecs.nodeLocZ[edge_end])*(hostSetInfoVecs.nodeLocZ[edge_start] - hostSetInfoVecs.nodeLocZ[edge_end]));
            if (R_TAIL_HEAD < generalParams.abs_Rmin){
                //rep_energy_TAIL_HEAD = linearSpringInfoVecs.spring_constant_rep1*(1.0-exp(-linearSpringInfoVecs.spring_constant_rep2*(R_TAIL_HEAD - generalParams.abs_Rmin)))*(1.0-exp(-linearSpringInfoVecs.spring_constant_rep2*(R_TAIL_HEAD - generalParams.abs_Rmin)));
                if (hostSetInfoVecs.nodes_in_upperhem[edge_start] == 1 && hostSetInfoVecs.nodes_in_upperhem[edge_end] == 1){
                    rep_energy_TAIL_HEAD = (linearSpringInfoVecs.spring_constant_weak/2.0)*(R_TAIL_HEAD - generalParams.abs_Rmin)*(R_TAIL_HEAD - generalParams.abs_Rmin);
                }
                else{
                    rep_energy_TAIL_HEAD = (linearSpringInfoVecs.spring_constant/2.0)*(R_TAIL_HEAD - generalParams.abs_Rmin)*(R_TAIL_HEAD - generalParams.abs_Rmin);
                }
                }
            else{
                rep_energy_TAIL_HEAD = 0.0;
                }
            //double rep_energy_T1, rep_energy_T2, rep_energy_H1, rep_energy_H2;
            rep_1 = rep_energy_TAIL_HEAD;
                
            double prob;
            double random_number;
            double Edif;
            if (DISTANCE >= 0.0){
                //WARNING: RESET BENDING COUNTER
                temp_bend = 0.0;


        
                double N1_vec1x, N1_vec1y, N1_vec1z, N1_vec2x, N1_vec2y, N1_vec2z, N2_vec1x, N2_vec1y, N2_vec1z, N2_vec2x, N2_vec2y, N2_vec2z;
                double N1_x, N1_y, N1_z, N2_x, N2_y, N2_z;
                bool THIS_SHOULD_NOT_HAPPEN = false;
                for (int j = 0; j < 5; j++){
                  /*  if (generalParams.edges_in_upperhem[edges_iteration[j]] == 1){
                        bend_spring_constant = bendingTriangleInfoVecs.spring_constant_weak;
                    }
                    else if (generalParams.edges_in_upperhem[edges_iteration[j]] == 0){
                        bend_spring_constant = (bendingTriangleInfoVecs.spring_constant_weak + bendingTriangleInfoVecs.spring_constant)/2.0;
                    }
                    else{
                        bend_spring_constant = bendingTriangleInfoVecs.spring_constant;
                    }*/
                    if (generalParams.SCALE_TYPE == 0){
                        bend_spring_constant = bendingTriangleInfoVecs.spring_constant*(1.0 - ((1.0/sqrt(2*3.14159*generalParams.gausssigma))*exp(-pow(hostSetInfoVecs.scaling_per_edge[edges_iteration[j]],2.0)/generalParams.gausssigma)));
                        if (bend_spring_constant < bendingTriangleInfoVecs.spring_constant_weak){bend_spring_constant = bendingTriangleInfoVecs.spring_constant_weak;};
                       /*if (generalParams.edges_in_upperhem[edges_iteration[j]] == 1){
                            bend_spring_constant = bendingTriangleInfoVecs.spring_constant_weak;
                        }
                        else if (generalParams.edges_in_upperhem[edges_iteration[j]] == 0){
                            bend_spring_constant = (bendingTriangleInfoVecs.spring_constant_weak + bendingTriangleInfoVecs.spring_constant)/2.0;
                        }
                        else{
                            bend_spring_constant = bendingTriangleInfoVecs.spring_constant;
                        }*/
                    }
                    else if (generalParams.SCALE_TYPE == 1){
                        bend_spring_constant = (bendingTriangleInfoVecs.spring_constant_weak*4.0)*pow(hostSetInfoVecs.scaling_per_edge[edges_iteration[j]],generalParams.scaling_pow) +
                                            (bendingTriangleInfoVecs.spring_constant_weak)*(1.0 - pow(hostSetInfoVecs.scaling_per_edge[edges_iteration[j]],generalParams.scaling_pow));
                        // bend_spring_constant = bendingTriangleInfoVecs.spring_constant*pow(hostSetInfoVecs.scaling_per_edge[edges_iteration[j]],generalParams.scaling_pow) +
                        //                     bendingTriangleInfoVecs.spring_constant_weak*(1.0 - pow(hostSetInfoVecs.scaling_per_edge[edges_iteration[j]],generalParams.scaling_pow));
                    }
                    else if (generalParams.SCALE_TYPE == 2){
                        bend_spring_constant = bendingTriangleInfoVecs.spring_constant_weak - 
                                            (bendingTriangleInfoVecs.spring_constant - bendingTriangleInfoVecs.spring_constant_weak)*
                                            hostSetInfoVecs.scaling_per_edge[edges_iteration[j]];
                        if (bend_spring_constant < bendingTriangleInfoVecs.spring_constant_weak){bend_spring_constant = bendingTriangleInfoVecs.spring_constant_weak;}
                    }
                    else if (generalParams.SCALE_TYPE == 3){
                        if (hostSetInfoVecs.edges_in_upperhem[edges_iteration[j]] == 1){// && hostSetInfoVecs.edges_in_tip[edges_iteration[j]] == 1){
                            bend_spring_constant = bendingTriangleInfoVecs.spring_constant_weak;
                            preferred_angle = bendingTriangleInfoVecs.initial_angle_bud;
                        }
                        else if (hostSetInfoVecs.edges_in_upperhem[edges_iteration[j]] == 0){//1 && hostSetInfoVecs.edges_in_tip[edges_iteration[j]] != 1){
                            // bend_spring_constant = (bendingTriangleInfoVecs.spring_constant_weak*2.0);
                            bend_spring_constant = (bendingTriangleInfoVecs.spring_constant_weak + bendingTriangleInfoVecs.spring_constant)/2.0;
                            preferred_angle = (bendingTriangleInfoVecs.initial_angle + bendingTriangleInfoVecs.initial_angle_bud)/2.0;
                        }
                        else{
                            bend_spring_constant = bendingTriangleInfoVecs.spring_constant;
                            preferred_angle = bendingTriangleInfoVecs.initial_angle;
                        }
                    }
                    else if (generalParams.SCALE_TYPE == 4){
                      //double scaling = 0.0;//bendingTriangleInfoVecs.spring_constant_weak/bendingTriangleInfoVecs.spring_constant;
                        //bend_spring_constant = bendingTriangleInfoVecs.spring_constant*((1.0/(1.0+pow(generalParams.hilleqnconst/hostSetInfoVecs.scaling_per_edge[edges_iteration[j]], generalParams.hilleqnpow)))*(1-scaling) + scaling);
                        double spectrum = bendingTriangleInfoVecs.spring_constant - bendingTriangleInfoVecs.spring_constant_weak;
                        bend_spring_constant = bendingTriangleInfoVecs.spring_constant_weak + ((1.0/(1.0+pow(generalParams.hilleqnconst/hostSetInfoVecs.scaling_per_edge[edges_iteration[j]], generalParams.hilleqnpow)))*spectrum);
                        if (bend_spring_constant < bendingTriangleInfoVecs.spring_constant_weak){bend_spring_constant = bendingTriangleInfoVecs.spring_constant_weak;}
                      }

                        if (j == 0){
                            N1_vec1x = hostSetInfoVecs.nodeLocX[HEAD] - hostSetInfoVecs.nodeLocX[TAIL];//x component of the 1st vector to calculate N1
                            N1_vec1y = hostSetInfoVecs.nodeLocY[HEAD] - hostSetInfoVecs.nodeLocY[TAIL];
                            N1_vec1z = hostSetInfoVecs.nodeLocZ[HEAD] - hostSetInfoVecs.nodeLocZ[TAIL];
                            N1_vec2x = hostSetInfoVecs.nodeLocX[edge_start] - hostSetInfoVecs.nodeLocX[TAIL];//x component of the 2nd vector to calculate N1
                            N1_vec2y = hostSetInfoVecs.nodeLocY[edge_start] - hostSetInfoVecs.nodeLocY[TAIL];
                            N1_vec2z = hostSetInfoVecs.nodeLocZ[edge_start] - hostSetInfoVecs.nodeLocZ[TAIL];
                            N2_vec1x = hostSetInfoVecs.nodeLocX[edge_end] - hostSetInfoVecs.nodeLocX[TAIL];
                            N2_vec1y = hostSetInfoVecs.nodeLocY[edge_end] - hostSetInfoVecs.nodeLocY[TAIL];
                            N2_vec1z = hostSetInfoVecs.nodeLocZ[edge_end] - hostSetInfoVecs.nodeLocZ[TAIL];
                            N2_vec2x = hostSetInfoVecs.nodeLocX[HEAD] - hostSetInfoVecs.nodeLocX[TAIL];
                            N2_vec2y = hostSetInfoVecs.nodeLocY[HEAD] - hostSetInfoVecs.nodeLocY[TAIL];
                            N2_vec2z = hostSetInfoVecs.nodeLocZ[HEAD] - hostSetInfoVecs.nodeLocZ[TAIL];

                            P0x_vol1 = hostSetInfoVecs.nodeLocX[HEAD];
                            P0y_vol1 = hostSetInfoVecs.nodeLocY[HEAD];
                            P0z_vol1 = hostSetInfoVecs.nodeLocZ[HEAD];
                            P0x_vol2 = hostSetInfoVecs.nodeLocX[HEAD];
                            P0y_vol2 = hostSetInfoVecs.nodeLocY[HEAD];
                            P0z_vol2 = hostSetInfoVecs.nodeLocZ[HEAD];
                        }
                        else if (j == 1){
                            N1_vec1x = H1t1_vec1x;
                            //std::cout<<"N1_vec1x = "<<N1_vec1x<<std::endl;
                            N1_vec1y = H1t1_vec1y;
                            //std::cout<<"N1_vec1y = "<<N1_vec1y<<std::endl;
                            N1_vec1z = H1t1_vec1z;
                            //std::cout<<"N1_vec1z = "<<N1_vec1z<<std::endl;
                            N1_vec2x = H1t1_vec2x;
                            //std::cout<<"N1_vec2x = "<<N1_vec2x<<std::endl;
                            N1_vec2y = H1t1_vec2y;
                            //std::cout<<"N1_vec2y = "<<N1_vec2y<<std::endl;
                            N1_vec2z = H1t1_vec2z;
                            //std::cout<<"N1_vec2z = "<<N1_vec2z<<std::endl;
                            N2_vec1x = H1t2_vec1x;
                            //std::cout<<"N2_vec1x = "<<N2_vec1x<<std::endl;
                            N2_vec1y = H1t2_vec1y;
                            //std::cout<<"N2_vec1y = "<<N2_vec1y<<std::endl;
                            N2_vec1z = H1t2_vec1z;
                            //std::cout<<"N2_vec1z = "<<N2_vec1z<<std::endl;
                            N2_vec2x = H1t2_vec2x;
                            //std::cout<<"N2_vec2x = "<<N2_vec2x<<std::endl;
                            N2_vec2y = H1t2_vec2y;
                            //std::cout<<"N2_vec2y = "<<N2_vec2y<<std::endl;
                            N2_vec2z = H1t2_vec2z;
                            //std::cout<<"N2_vec2z = "<<N2_vec2z<<std::endl;
                        }
                        else if (j == 2){
                            N1_vec1x = H2t1_vec1x;
                            N1_vec1y = H2t1_vec1y;
                            N1_vec1z = H2t1_vec1z;
                            N1_vec2x = H2t1_vec2x;
                            N1_vec2y = H2t1_vec2y;
                            N1_vec2z = H2t1_vec2z;
                            N2_vec1x = H2t2_vec1x;
                            N2_vec1y = H2t2_vec1y;
                            N2_vec1z = H2t2_vec1z;
                            N2_vec2x = H2t2_vec2x;
                            N2_vec2y = H2t2_vec2y;
                            N2_vec2z = H2t2_vec2z;
                        }
                        else if (j == 3){
                            N1_vec1x = T1t1_vec1x;
                            N1_vec1y = T1t1_vec1y;
                            N1_vec1z = T1t1_vec1z;
                            N1_vec2x = T1t1_vec2x;
                            N1_vec2y = T1t1_vec2y;
                            N1_vec2z = T1t1_vec2z;
                            N2_vec1x = T1t2_vec1x;
                            N2_vec1y = T1t2_vec1y;
                            N2_vec1z = T1t2_vec1z;
                            N2_vec2x = T1t2_vec2x;
                            N2_vec2y = T1t2_vec2y;
                            N2_vec2z = T1t2_vec2z;
                        }
                        else if (j == 4){
                            N1_vec1x = T2t1_vec1x;
                            N1_vec1y = T2t1_vec1y;
                            N1_vec1z = T2t1_vec1z;
                            N1_vec2x = T2t1_vec2x;
                            N1_vec2y = T2t1_vec2y;
                            N1_vec2z = T2t1_vec2z;
                            N2_vec1x = T2t2_vec1x;
                            N2_vec1y = T2t2_vec1y;
                            N2_vec1z = T2t2_vec1z;
                            N2_vec2x = T2t2_vec2x;
                            N2_vec2y = T2t2_vec2y;
                            N2_vec2z = T2t2_vec2z;
                        }                        
                        //std::vector<double> N1(3);
                        N1_x = N1_vec1y*N1_vec2z - N1_vec2y*N1_vec1z;
                        N1_y = -(N1_vec1x*N1_vec2z - N1_vec2x*N1_vec1z);
                        N1_z = N1_vec1x*N1_vec2y - N1_vec2x*N1_vec1y;
                        double nN1 = sqrt(pow(N1_x,2)+pow(N1_y,2)+pow(N1_z,2));
                        ////std::cout<<"newnN1 = "<<nN1<<std::endl;

        
                        //std::vector<double> N2(3);
                        N2_x = N2_vec1y*N2_vec2z - N2_vec2y*N2_vec1z;
                        N2_y = -(N2_vec1x*N2_vec2z - N2_vec2x*N2_vec1z);
                        N2_z = N2_vec1x*N2_vec2y - N2_vec2x*N2_vec1y; 
                        double nN2 = sqrt(pow(N2_x,2)+pow(N2_y,2)+pow(N2_z,2));;
                        ////std::cout<<"newnN2 = "<<nN2<<std::endl;

                        if (j == 0){
                            N1x_vol = N1_x/nN1;
                            N1y_vol = N1_y/nN1;
                            N1z_vol = N1_z/nN1;
                            N2x_vol = N2_x/nN2;
                            N2y_vol = N2_y/nN2;
                            N2z_vol = N2_z/nN2;
                        }

                        double cosAngle = (N1_x*N2_x + N1_y*N2_y + N1_z*N2_z)/(nN1*nN2);
                        ////std::cout<<"cosAngle = "<<cosAngle<<std::endl;
                        
                        if (cosAngle > 1.0) {
                            cosAngle = 1.0;
                        }
                        else if (cosAngle < -1.0){
                            cosAngle = -1.0;
                        }
                        if (cosAngle == -1.0){
                            THIS_SHOULD_NOT_HAPPEN = true;
                        }

                        double theta_current = acos( cosAngle );
                        
                        double local_energy = bend_spring_constant * (1 - cos(theta_current - preferred_angle) );
                        temp_bend = temp_bend + local_energy;
                        //std::cout<<"bending energy "<<local_energy<<std::endl;
                        //for (int COUNT = 0; COUNT < 3; COUNT++){
                        ////std::cout<<"unit normal 1 = "<<N1[COUNT]<<std::endl;
                        ////std::cout<<"unit normal 2 = "<<N2[COUNT]<<std::endl;}
                        //std::cout<<"angle "<<theta_current<<std::endl;
                    }
                double bend_1 = temp_bend;

                for (int t = 0; t < 4; t++){    

                if (t == 0){
                    midpt = edge_start;
                    endpt1 = HEAD;
                    endpt2 = TAIL;
                    endpt3 = edge_end;
                    true_endpoints = 2;
                    falsept1 = node_index_H1;
                    falsept2 = node_index_T1;
                }
                else if (t == 1){
                    midpt = HEAD;
                    endpt1 = edge_start;
                    endpt2 = edge_end;
                    endpt3 = TAIL;
                    true_endpoints = 3;
                    falsept1 = node_index_H1;
                    falsept2 = node_index_H2;
                }
                else if (t == 2){
                    midpt = edge_end;
                    endpt1 = HEAD;
                    endpt2 = edge_start;
                    endpt3 = TAIL;
                    true_endpoints = 2;
                    falsept1 = node_index_H2;
                    falsept2 = node_index_T2;
                }
                else if (t == 3){
                    midpt = TAIL;
                    endpt1 = edge_start;
                    endpt2 = edge_end;
                    endpt3 = HEAD;
                    true_endpoints = 3;
                    falsept1 = node_index_T1;
                    falsept2 = node_index_T2;
                }
                nb1 = hostSetInfoVecs.nndata1[midpt];
                nb2 = hostSetInfoVecs.nndata2[midpt];
                nb3 = hostSetInfoVecs.nndata3[midpt];
                nb4 = hostSetInfoVecs.nndata4[midpt];
                nb5 = hostSetInfoVecs.nndata5[midpt];
                nb6 = hostSetInfoVecs.nndata6[midpt];
                nb7 = hostSetInfoVecs.nndata7[midpt];
                nb8 = hostSetInfoVecs.nndata8[midpt];
                nb9 = hostSetInfoVecs.nndata9[midpt];
                for (int y = 0; y < 9; y++){
                    int startpt;
                    if (y == 0 && nb1>= 0 && nb1 != endpt1 && nb1 != endpt2 && nb1 != endpt3 && nb1 != falsept1 && nb1 != falsept2){startpt = nb1;}
                    else if (y == 1 && nb2>= 0 && nb2 != endpt1 && nb2 != endpt2 && nb2 != endpt3 && nb2 != falsept1 && nb2 != falsept2){startpt = nb2;}
                    else if (y == 2 && nb3>= 0 && nb3 != endpt1 && nb3 != endpt2 && nb3 != endpt3 && nb3 != falsept1 && nb3 != falsept2){startpt = nb3;}
                    else if (y == 3 && nb4>= 0 && nb4 != endpt1 && nb4 != endpt2 && nb4 != endpt3 && nb4 != falsept1 && nb4 != falsept2){startpt = nb4;}
                    else if (y == 4 && nb5>= 0 && nb5 != endpt1 && nb5 != endpt2 && nb5 != endpt3 && nb5 != falsept1 && nb5 != falsept2){startpt = nb5;}
                    else if (y == 5 && nb6>= 0 && nb6 != endpt1 && nb6 != endpt2 && nb6 != endpt3 && nb6 != falsept1 && nb6 != falsept2){startpt = nb6;}
                    else if (y == 6 && nb7>= 0 && nb7 != endpt1 && nb7 != endpt2 && nb7 != endpt3 && nb7 != falsept1 && nb7 != falsept2){startpt = nb7;}
                    else if (y == 7 && nb8>= 0 && nb8 != endpt1 && nb8 != endpt2 && nb8 != endpt3 && nb8 != falsept1 && nb8 != falsept2){startpt = nb8;}
                    else if (y == 8 && nb9>= 0 && nb9 != endpt1 && nb9 != endpt2 && nb9 != endpt3 && nb9 != falsept1 && nb9 != falsept2){startpt = nb9;}
                    else{continue;}

                    if (true_endpoints == 3){
                        for (int h = 0; h < 3; h++){
                            int aa;
                            if (h==0){aa = endpt1;}
                            else if (h==1){aa = endpt2;}
                            else if (h==2){aa = endpt3;}
                            R_TAIL_HEAD = sqrt((hostSetInfoVecs.nodeLocX[startpt] - hostSetInfoVecs.nodeLocX[aa])*(hostSetInfoVecs.nodeLocX[startpt] - hostSetInfoVecs.nodeLocX[aa]) +
                                    (hostSetInfoVecs.nodeLocY[startpt] - hostSetInfoVecs.nodeLocY[aa])*(hostSetInfoVecs.nodeLocY[startpt] - hostSetInfoVecs.nodeLocY[aa]) +
                                    (hostSetInfoVecs.nodeLocZ[startpt] - hostSetInfoVecs.nodeLocZ[aa])*(hostSetInfoVecs.nodeLocZ[startpt] - hostSetInfoVecs.nodeLocZ[aa]));
                            if (R_TAIL_HEAD < generalParams.abs_Rmin){
                                //rep_energy_TAIL_HEAD = linearSpringInfoVecs.spring_constant_rep1*(1.0-exp(-linearSpringInfoVecs.spring_constant_rep2*(R_TAIL_HEAD - generalParams.abs_Rmin)))*
                                //                        (1.0-exp(-linearSpringInfoVecs.spring_constant_rep2*(R_TAIL_HEAD - generalParams.abs_Rmin)));
                                if (hostSetInfoVecs.nodes_in_upperhem[startpt] == 1 && hostSetInfoVecs.nodes_in_upperhem[aa] == 1){
                                    rep_energy_TAIL_HEAD = (linearSpringInfoVecs.spring_constant_weak/2.0)*(R_TAIL_HEAD - generalParams.abs_Rmin)*(R_TAIL_HEAD - generalParams.abs_Rmin);
                                }
                                else{
                                    rep_energy_TAIL_HEAD = (linearSpringInfoVecs.spring_constant/2.0)*(R_TAIL_HEAD - generalParams.abs_Rmin)*(R_TAIL_HEAD - generalParams.abs_Rmin);
                                }
                            }
                            else{rep_energy_TAIL_HEAD = 0.0;}
                            rep_1 += rep_energy_TAIL_HEAD;
                        }
                    }
                    else if (true_endpoints == 2){
                        for (int h = 0; h < 2; h++){
                            int aa;
                            if (h==0){aa = endpt1;}
                            else if (h==1){aa = endpt2;}
                            R_TAIL_HEAD = sqrt((hostSetInfoVecs.nodeLocX[startpt] - hostSetInfoVecs.nodeLocX[aa])*(hostSetInfoVecs.nodeLocX[startpt] - hostSetInfoVecs.nodeLocX[aa]) +
                                    (hostSetInfoVecs.nodeLocY[startpt] - hostSetInfoVecs.nodeLocY[aa])*(hostSetInfoVecs.nodeLocY[startpt] - hostSetInfoVecs.nodeLocY[aa]) +
                                    (hostSetInfoVecs.nodeLocZ[startpt] - hostSetInfoVecs.nodeLocZ[aa])*(hostSetInfoVecs.nodeLocZ[startpt] - hostSetInfoVecs.nodeLocZ[aa]));
                            if (R_TAIL_HEAD < generalParams.abs_Rmin){
                                //rep_energy_TAIL_HEAD = linearSpringInfoVecs.spring_constant_rep1*(1.0-exp(-linearSpringInfoVecs.spring_constant_rep2*(R_TAIL_HEAD - generalParams.abs_Rmin)))*
                                //                        (1.0-exp(-linearSpringInfoVecs.spring_constant_rep2*(R_TAIL_HEAD - generalParams.abs_Rmin)));
                                if (hostSetInfoVecs.nodes_in_upperhem[startpt] == 1 && hostSetInfoVecs.nodes_in_upperhem[aa] == 1){
                                    rep_energy_TAIL_HEAD = (linearSpringInfoVecs.spring_constant_weak/2.0)*(R_TAIL_HEAD - generalParams.abs_Rmin)*(R_TAIL_HEAD - generalParams.abs_Rmin);
                                }
                                else{
                                    rep_energy_TAIL_HEAD = (linearSpringInfoVecs.spring_constant/2.0)*(R_TAIL_HEAD - generalParams.abs_Rmin)*(R_TAIL_HEAD - generalParams.abs_Rmin);
                                }
                            }
                            else{rep_energy_TAIL_HEAD = 0.0;}
                            rep_1 += rep_energy_TAIL_HEAD;
                        }
                    }
                }
                
                //NOW WE UTILIZE THIS LOOP TO CALCULATE THE REPULSION ENERGY ASSOCIATED WITH THE FALSEPTS (I.E. THE ONES ASSOCIATED WITH H10, T10, H20, T20 BUT NOT IN THE SUBSYSTEM DEPECTED ABOVE)
                // int oiu1, oiu2, oiu3;
                // if (t == 0){oiu1 = node_index_H1; oiu2 = edge_end; oiu3 = TAIL;}
                // else if (t == 1){oiu1 = node_index_H2; oiu2 = edge_start; oiu3 = TAIL;}
                // else if (t == 2){oiu1 = node_index_T1; oiu2 = edge_end; oiu3 = HEAD;}
                // else if (t == 3){oiu1 = node_index_T2; oiu2 = edge_start; oiu3 = HEAD;}
                // R_TAIL_HEAD = sqrt((hostSetInfoVecs.nodeLocX[oiu1] - hostSetInfoVecs.nodeLocX[oiu2])*(hostSetInfoVecs.nodeLocX[oiu1] - hostSetInfoVecs.nodeLocX[oiu2]) +
                //         (hostSetInfoVecs.nodeLocY[oiu1] - hostSetInfoVecs.nodeLocY[oiu2])*(hostSetInfoVecs.nodeLocY[oiu1] - hostSetInfoVecs.nodeLocY[oiu2]) +
                //         (hostSetInfoVecs.nodeLocZ[oiu1] - hostSetInfoVecs.nodeLocZ[oiu2])*(hostSetInfoVecs.nodeLocZ[oiu1] - hostSetInfoVecs.nodeLocZ[oiu2]));
                // if (R_TAIL_HEAD < generalParams.abs_Rmin){
                //     rep_energy_TAIL_HEAD = linearSpringInfoVecs.spring_constant_rep1*(1.0-exp(-linearSpringInfoVecs.spring_constant_rep2*(R_TAIL_HEAD - generalParams.abs_Rmin)))*
                //                             (1.0-exp(-linearSpringInfoVecs.spring_constant_rep2*(R_TAIL_HEAD - generalParams.abs_Rmin)));
                // }
                // else{rep_energy_TAIL_HEAD = 0.0;}
                // rep_1 += rep_energy_TAIL_HEAD;
                // R_TAIL_HEAD = sqrt((hostSetInfoVecs.nodeLocX[oiu1] - hostSetInfoVecs.nodeLocX[oiu3])*(hostSetInfoVecs.nodeLocX[oiu1] - hostSetInfoVecs.nodeLocX[oiu3]) +
                //         (hostSetInfoVecs.nodeLocY[oiu1] - hostSetInfoVecs.nodeLocY[oiu3])*(hostSetInfoVecs.nodeLocY[oiu1] - hostSetInfoVecs.nodeLocY[oiu3]) +
                //         (hostSetInfoVecs.nodeLocZ[oiu1] - hostSetInfoVecs.nodeLocZ[oiu3])*(hostSetInfoVecs.nodeLocZ[oiu1] - hostSetInfoVecs.nodeLocZ[oiu3]));
                // if (R_TAIL_HEAD < generalParams.abs_Rmin){
                //     rep_energy_TAIL_HEAD = linearSpringInfoVecs.spring_constant_rep1*(1.0-exp(-linearSpringInfoVecs.spring_constant_rep2*(R_TAIL_HEAD - generalParams.abs_Rmin)))*
                //                             (1.0-exp(-linearSpringInfoVecs.spring_constant_rep2*(R_TAIL_HEAD - generalParams.abs_Rmin)));
                // }
                // else{rep_energy_TAIL_HEAD = 0.0;}
                // rep_1 += rep_energy_TAIL_HEAD;
            }

                int H0n1 = HEAD;
                int H0n2 = edge_start;
                int H0n3 = TAIL;
                int T0n1 = TAIL;
                int T0n2 = edge_end;
                int T0n3 = HEAD;
                double a = sqrt(pow((hostSetInfoVecs.nodeLocX[H0n2] - hostSetInfoVecs.nodeLocX[H0n1]),2.0) + 
                        pow((hostSetInfoVecs.nodeLocY[H0n2] - hostSetInfoVecs.nodeLocY[H0n1]),2.0) +
                        pow((hostSetInfoVecs.nodeLocZ[H0n2] - hostSetInfoVecs.nodeLocZ[H0n1]),2.0)
                        );
                double b = sqrt(pow((hostSetInfoVecs.nodeLocX[H0n3] - hostSetInfoVecs.nodeLocX[H0n1]),2.0) + 
                            pow((hostSetInfoVecs.nodeLocY[H0n3] - hostSetInfoVecs.nodeLocY[H0n1]),2.0) +
                            pow((hostSetInfoVecs.nodeLocZ[H0n3] - hostSetInfoVecs.nodeLocZ[H0n1]),2.0)
                            );
                double c = sqrt(pow((hostSetInfoVecs.nodeLocX[H0n3] - hostSetInfoVecs.nodeLocX[H0n2]),2.0) + 
                            pow((hostSetInfoVecs.nodeLocY[H0n3] - hostSetInfoVecs.nodeLocY[H0n2]),2.0) +
                            pow((hostSetInfoVecs.nodeLocZ[H0n3] - hostSetInfoVecs.nodeLocZ[H0n2]),2.0)
                            );
                double mean_abc = (a + b + c)/2;
                double d = sqrt(pow((hostSetInfoVecs.nodeLocX[T0n2] - hostSetInfoVecs.nodeLocX[T0n1]),2.0) + 
                            pow((hostSetInfoVecs.nodeLocY[T0n2] - hostSetInfoVecs.nodeLocY[T0n1]),2.0) +
                            pow((hostSetInfoVecs.nodeLocZ[T0n2] - hostSetInfoVecs.nodeLocZ[T0n1]),2.0)
                            );
                double e = sqrt(pow((hostSetInfoVecs.nodeLocX[T0n3] - hostSetInfoVecs.nodeLocX[T0n1]),2.0) + 
                            pow((hostSetInfoVecs.nodeLocY[T0n3] - hostSetInfoVecs.nodeLocY[T0n1]),2.0) +
                            pow((hostSetInfoVecs.nodeLocZ[T0n3] - hostSetInfoVecs.nodeLocZ[T0n1]),2.0)
                            );
                double f = sqrt(pow((hostSetInfoVecs.nodeLocX[T0n3] - hostSetInfoVecs.nodeLocX[T0n2]),2.0) + 
                            pow((hostSetInfoVecs.nodeLocY[T0n3] - hostSetInfoVecs.nodeLocY[T0n2]),2.0) +
                            pow((hostSetInfoVecs.nodeLocZ[T0n3] - hostSetInfoVecs.nodeLocZ[T0n2]),2.0)
                            );
                double mean_def = (d + e + f)/2.0;
                /*if (generalParams.triangles_in_upperhem[H0] == 1){
                    area_spring_constant_1 = areaTriangleInfoVecs.spring_constant_weak;
                }
                else if (generalParams.triangles_in_upperhem[H0] == 0){
                    area_spring_constant_1 = (areaTriangleInfoVecs.spring_constant_weak + areaTriangleInfoVecs.spring_constant)/2.0;
                }
                else{
                    area_spring_constant_1 = areaTriangleInfoVecs.spring_constant;
                }*/
                if (generalParams.SCALE_TYPE == 0){
                    area_spring_constant_1 = areaTriangleInfoVecs.spring_constant*((1.0 - ((1.0/sqrt(2*3.14159*generalParams.gausssigma))*exp(-pow(hostSetInfoVecs.scaling_per_edge[iedge],2.0)/generalParams.gausssigma))) +
                                            (1.0 - ((1.0/sqrt(2*3.14159*generalParams.gausssigma))*exp(-pow(hostSetInfoVecs.scaling_per_edge[H1],2.0)/generalParams.gausssigma))) +
                                            (1.0 - ((1.0/sqrt(2*3.14159*generalParams.gausssigma))*exp(-pow(hostSetInfoVecs.scaling_per_edge[T1],2.0)/generalParams.gausssigma))))/3.0;
                                            if (area_spring_constant_1 < areaTriangleInfoVecs.spring_constant_weak){area_spring_constant_1 = areaTriangleInfoVecs.spring_constant_weak;}
                }
                else if (generalParams.SCALE_TYPE == 1 ){
                    area_spring_constant_1 =  ((areaTriangleInfoVecs.spring_constant_weak*2.0)*pow(hostSetInfoVecs.scaling_per_edge[iedge],generalParams.scaling_pow) + areaTriangleInfoVecs.spring_constant_weak*(1.0 - pow(hostSetInfoVecs.scaling_per_edge[iedge],generalParams.scaling_pow)) +
                                       (areaTriangleInfoVecs.spring_constant_weak*2.0)*pow(hostSetInfoVecs.scaling_per_edge[H1],generalParams.scaling_pow) + areaTriangleInfoVecs.spring_constant_weak*(1.0 - pow(hostSetInfoVecs.scaling_per_edge[H1],generalParams.scaling_pow)) +
                                       (areaTriangleInfoVecs.spring_constant_weak*2.0)*pow(hostSetInfoVecs.scaling_per_edge[T1],generalParams.scaling_pow) + areaTriangleInfoVecs.spring_constant_weak*(1.0 - pow(hostSetInfoVecs.scaling_per_edge[T1],generalParams.scaling_pow)))/3.0;
                    // area_spring_constant_1 =  (areaTriangleInfoVecs.spring_constant*pow(hostSetInfoVecs.scaling_per_edge[iedge],generalParams.scaling_pow) + areaTriangleInfoVecs.spring_constant_weak*(1.0 - pow(hostSetInfoVecs.scaling_per_edge[iedge],generalParams.scaling_pow)) +
                    //                     areaTriangleInfoVecs.spring_constant*pow(hostSetInfoVecs.scaling_per_edge[H1],generalParams.scaling_pow) + areaTriangleInfoVecs.spring_constant_weak*(1.0 - pow(hostSetInfoVecs.scaling_per_edge[H1],generalParams.scaling_pow)) +
                    //                     areaTriangleInfoVecs.spring_constant*pow(hostSetInfoVecs.scaling_per_edge[T1],generalParams.scaling_pow) + areaTriangleInfoVecs.spring_constant_weak*(1.0 - pow(hostSetInfoVecs.scaling_per_edge[T1],generalParams.scaling_pow)))/3.0;
                }
                else if (generalParams.SCALE_TYPE == 2){
                    area_spring_constant_1 = ((areaTriangleInfoVecs.spring_constant - (areaTriangleInfoVecs.spring_constant - areaTriangleInfoVecs.spring_constant_weak)*hostSetInfoVecs.scaling_per_edge[iedge]) +
                                        (areaTriangleInfoVecs.spring_constant - (areaTriangleInfoVecs.spring_constant - areaTriangleInfoVecs.spring_constant_weak)*hostSetInfoVecs.scaling_per_edge[H1]) +
                                        (areaTriangleInfoVecs.spring_constant - (areaTriangleInfoVecs.spring_constant - areaTriangleInfoVecs.spring_constant_weak)*hostSetInfoVecs.scaling_per_edge[T1]))/3.0;
                    if (area_spring_constant_1 < areaTriangleInfoVecs.spring_constant_weak){area_spring_constant_1 = areaTriangleInfoVecs.spring_constant_weak;}
                }
                else if (generalParams.SCALE_TYPE == 3){
                    if (hostSetInfoVecs.triangles_in_upperhem[H0] == 1){// && hostSetInfoVecs.triangles_in_tip[H0] == 1){
                    area_spring_constant_1 = areaTriangleInfoVecs.spring_constant_weak;
                    }
                    else if (hostSetInfoVecs.triangles_in_upperhem[H0] == 0){//1 && hostSetInfoVecs.triangles_in_tip[H0] != 1){
                        // area_spring_constant_1 = (areaTriangleInfoVecs.spring_constant_weak*2.0);
                        area_spring_constant_1 = (areaTriangleInfoVecs.spring_constant_weak + areaTriangleInfoVecs.spring_constant)/2.0;
                    }
                    else{
                        area_spring_constant_1 = areaTriangleInfoVecs.spring_constant;
                    }
                }
                else if (generalParams.SCALE_TYPE == 4){
                //double scaling = 0.0;//areaTriangleInfoVecs.spring_constant_weak/areaTriangleInfoVecs.spring_constant;
			    //area_spring_constant_1 = (areaTriangleInfoVecs.spring_constant*((1.0/(1.0+pow(generalParams.hilleqnconst/hostSetInfoVecs.scaling_per_edge[iedge], generalParams.hilleqnpow)))*(1-scaling) + scaling) +
                //                   areaTriangleInfoVecs.spring_constant*((1.0/(1.0+pow(generalParams.hilleqnconst/hostSetInfoVecs.scaling_per_edge[H1], generalParams.hilleqnpow)))*(1-scaling) + scaling) +
                //                   areaTriangleInfoVecs.spring_constant*((1.0/(1.0+pow(generalParams.hilleqnconst/hostSetInfoVecs.scaling_per_edge[T1], generalParams.hilleqnpow)))*(1-scaling) + scaling))/3.0;
			    double spectrum = areaTriangleInfoVecs.spring_constant - areaTriangleInfoVecs.spring_constant_weak;
                area_spring_constant_1 = (areaTriangleInfoVecs.spring_constant_weak + ((1.0/(1.0+pow(generalParams.hilleqnconst/hostSetInfoVecs.scaling_per_edge[iedge], generalParams.hilleqnpow)))*spectrum) +
                                          areaTriangleInfoVecs.spring_constant_weak + ((1.0/(1.0+pow(generalParams.hilleqnconst/hostSetInfoVecs.scaling_per_edge[H1], generalParams.hilleqnpow)))*spectrum) +
                                          areaTriangleInfoVecs.spring_constant_weak + ((1.0/(1.0+pow(generalParams.hilleqnconst/hostSetInfoVecs.scaling_per_edge[T1], generalParams.hilleqnpow)))*spectrum))/3.0;
                if (area_spring_constant_1 < areaTriangleInfoVecs.spring_constant_weak){area_spring_constant_1 = areaTriangleInfoVecs.spring_constant_weak;}
		      }
                /*if (generalParams.triangles_in_upperhem[T0] == 1){
                    area_spring_constant_2 = areaTriangleInfoVecs.spring_constant_weak;
                }
                else if (generalParams.triangles_in_upperhem[T0] == 0){
                    area_spring_constant_2 = (areaTriangleInfoVecs.spring_constant_weak + areaTriangleInfoVecs.spring_constant)/2.0;
                }
                else{
                    area_spring_constant_2 = areaTriangleInfoVecs.spring_constant;
                }*/
                if (generalParams.SCALE_TYPE == 0){
                    area_spring_constant_2 = areaTriangleInfoVecs.spring_constant*((1.0 - ((1.0/sqrt(2*3.14159*generalParams.gausssigma))*exp(-pow(hostSetInfoVecs.scaling_per_edge[iedge],2.0)/generalParams.gausssigma))) +
                                            (1.0 - ((1.0/sqrt(2*3.14159*generalParams.gausssigma))*exp(-pow(hostSetInfoVecs.scaling_per_edge[T2],2.0)/generalParams.gausssigma))) +
                                            (1.0 - ((1.0/sqrt(2*3.14159*generalParams.gausssigma))*exp(-pow(hostSetInfoVecs.scaling_per_edge[H2],2.0)/generalParams.gausssigma))))/3.0;
                                            if (area_spring_constant_2 < areaTriangleInfoVecs.spring_constant_weak){area_spring_constant_2 = areaTriangleInfoVecs.spring_constant_weak;}
                }
                else if (generalParams.SCALE_TYPE == 1){
                    area_spring_constant_2 = ((areaTriangleInfoVecs.spring_constant_weak*2.0)*pow(hostSetInfoVecs.scaling_per_edge[iedge],generalParams.scaling_pow) + areaTriangleInfoVecs.spring_constant_weak*(1.0 - pow(hostSetInfoVecs.scaling_per_edge[iedge],generalParams.scaling_pow)) +
                        (areaTriangleInfoVecs.spring_constant_weak*2.0)*pow(hostSetInfoVecs.scaling_per_edge[T2],generalParams.scaling_pow) + areaTriangleInfoVecs.spring_constant_weak*(1.0 - pow(hostSetInfoVecs.scaling_per_edge[T2],generalParams.scaling_pow)) +
                        (areaTriangleInfoVecs.spring_constant_weak*2.0)*pow(hostSetInfoVecs.scaling_per_edge[H2],generalParams.scaling_pow) + areaTriangleInfoVecs.spring_constant_weak*(1.0 - pow(hostSetInfoVecs.scaling_per_edge[H2],generalParams.scaling_pow)))/3.0;
                    // area_spring_constant_2 = (areaTriangleInfoVecs.spring_constant*pow(hostSetInfoVecs.scaling_per_edge[iedge],generalParams.scaling_pow) + areaTriangleInfoVecs.spring_constant_weak*(1.0 - pow(hostSetInfoVecs.scaling_per_edge[iedge],generalParams.scaling_pow)) +
                    //     areaTriangleInfoVecs.spring_constant*pow(hostSetInfoVecs.scaling_per_edge[T2],generalParams.scaling_pow) + areaTriangleInfoVecs.spring_constant_weak*(1.0 - pow(hostSetInfoVecs.scaling_per_edge[T2],generalParams.scaling_pow)) +
                    //     areaTriangleInfoVecs.spring_constant*pow(hostSetInfoVecs.scaling_per_edge[H2],generalParams.scaling_pow) + areaTriangleInfoVecs.spring_constant_weak*(1.0 - pow(hostSetInfoVecs.scaling_per_edge[H2],generalParams.scaling_pow)))/3.0;
                }
                else if (generalParams.SCALE_TYPE == 2){
                    area_spring_constant_2 = ((areaTriangleInfoVecs.spring_constant - (areaTriangleInfoVecs.spring_constant - areaTriangleInfoVecs.spring_constant_weak)*hostSetInfoVecs.scaling_per_edge[iedge]) +
                                        (areaTriangleInfoVecs.spring_constant - (areaTriangleInfoVecs.spring_constant - areaTriangleInfoVecs.spring_constant_weak)*hostSetInfoVecs.scaling_per_edge[T2] )+
                                        (areaTriangleInfoVecs.spring_constant - (areaTriangleInfoVecs.spring_constant - areaTriangleInfoVecs.spring_constant_weak)*hostSetInfoVecs.scaling_per_edge[H2]))/3.0;
                    if (area_spring_constant_2 < areaTriangleInfoVecs.spring_constant_weak){area_spring_constant_2 = areaTriangleInfoVecs.spring_constant_weak;}
                }
                else if (generalParams.SCALE_TYPE == 3){
                    if (hostSetInfoVecs.triangles_in_upperhem[T0] == 1){// && hostSetInfoVecs.triangles_in_tip[T0] == 1){
                    area_spring_constant_2 = areaTriangleInfoVecs.spring_constant_weak;
                    }
                    else if (hostSetInfoVecs.triangles_in_upperhem[T0] == 0){//1 && hostSetInfoVecs.triangles_in_tip[T0] != 1){
                        // area_spring_constant_2 = (areaTriangleInfoVecs.spring_constant_weak*2.0);
                        area_spring_constant_2 = (areaTriangleInfoVecs.spring_constant_weak + areaTriangleInfoVecs.spring_constant)/2.0;
                    }
                    else{
                        area_spring_constant_2 = areaTriangleInfoVecs.spring_constant;
                    }
                }
                else if (generalParams.SCALE_TYPE == 4){
                  //double scaling = 0.0;// areaTriangleInfoVecs.spring_constant_weak/areaTriangleInfoVecs.spring_constant;
			        //area_spring_constant_2 = (areaTriangleInfoVecs.spring_constant*((1.0/(1.0+pow(generalParams.hilleqnconst/hostSetInfoVecs.scaling_per_edge[iedge], generalParams.hilleqnpow)))*(1-scaling) + scaling) +
                    //               areaTriangleInfoVecs.spring_constant*((1.0/(1.0+pow(generalParams.hilleqnconst/hostSetInfoVecs.scaling_per_edge[T2], generalParams.hilleqnpow)))*(1-scaling) + scaling) +
                    //               areaTriangleInfoVecs.spring_constant*((1.0/(1.0+pow(generalParams.hilleqnconst/hostSetInfoVecs.scaling_per_edge[H2], generalParams.hilleqnpow)))*(1-scaling) + scaling))/3.0;
			        double spectrum = areaTriangleInfoVecs.spring_constant - areaTriangleInfoVecs.spring_constant_weak;
                    area_spring_constant_2 = (areaTriangleInfoVecs.spring_constant_weak + ((1.0/(1.0+pow(generalParams.hilleqnconst/hostSetInfoVecs.scaling_per_edge[iedge], generalParams.hilleqnpow)))*spectrum) +
                                              areaTriangleInfoVecs.spring_constant_weak + ((1.0/(1.0+pow(generalParams.hilleqnconst/hostSetInfoVecs.scaling_per_edge[T2], generalParams.hilleqnpow)))*spectrum) +
                                              areaTriangleInfoVecs.spring_constant_weak + ((1.0/(1.0+pow(generalParams.hilleqnconst/hostSetInfoVecs.scaling_per_edge[H2], generalParams.hilleqnpow)))*spectrum))/3.0;
                    if (area_spring_constant_2 < areaTriangleInfoVecs.spring_constant_weak){area_spring_constant_2 = areaTriangleInfoVecs.spring_constant_weak;}
		   }
                double area_H0 = sqrt(mean_abc*(mean_abc - a)*(mean_abc - b)*(mean_abc - c));
                double area_T0 = sqrt(mean_def*(mean_def - d)*(mean_def - e)*(mean_def - f));
                double area_1_energy = area_spring_constant_1*pow((area_H0 - areaTriangleInfoVecs.initial_area),2.0)/(2*areaTriangleInfoVecs.initial_area) +
                                    area_spring_constant_2*pow((area_T0 - areaTriangleInfoVecs.initial_area),2.0)/(2*areaTriangleInfoVecs.initial_area);
                //double vol_H0 = (1.0/3.0)*(P0x_vol1*N1x_vol + P0y_vol1*N1y_vol + P0z_vol1*N1z_vol)*area_H0;
                //double vol_T0 = (1.0/3.0)*(P0x_vol2*N2x_vol + P0y_vol2*N2y_vol + P0z_vol2*N2z_vol)*area_T0;
                //vol_1 = vol_H0 + vol_T0;
                //double new_vol = generalParams.true_current_total_volume + (vol_1 - vol_0);
                //double new_vol_energy = generalParams.volume_spring_constant*(new_vol - generalParams.eq_total_volume)*(new_vol - generalParams.eq_total_volume)/
                //                        (2.0*generalParams.Rmin*generalParams.Rmin*generalParams.Rmin*generalParams.eq_total_volume);
                double E_1 = linear_1 + bend_1 + area_1_energy + rep_1;// + new_vol_energy;
                // std::cout<<"new linear energy = "<<linear_1<<" , new length = "<<DISTANCE<<std::endl;
                // std::cout<<"new bend energy = "<<bend_1<<std::endl;
                // std::cout<<"new area energy = "<<area_1_energy<<std::endl;
                // std::cout<<"new total energy: "<<E_1<<std::endl;
            
            //Now compute the Boltzmann factor to determine if a swap occurs.
            std::random_device rd;  //Will be used to obtain a seed for the random number engine
            std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
            std::uniform_real_distribution<> dis(0.0, 1.0);
            random_number = dis(gen);
            //double random_number = 0;
            Edif = (E_1 - E_0);
            

            prob = generalParams.tau*exp(-Edif/generalParams.kT);
            
            }
            else{
                prob = -1.0;
            }
            // std::cout<<"P(swap): "<<prob<<std::endl;
            
            bool ACCEPT2;
            if (!isnan(prob)){
                if (prob >= 1){ACCEPT2 = true;}
                else if (prob < 1 && random_number <= prob){ACCEPT2 = true;}
                else if (prob < 1 && random_number > prob){ACCEPT2 = false;}
            }
            else{ACCEPT2 = false;}
            ////std::cout<<"ACCEPT2 = "<<ACCEPT2<<std::endl;
            //Perform real update
            //if (ACCEPT2 == true){
            //if (Edif < 100.0){
            if (ACCEPT2 == true ){//&& THIS_SHOULD_NOT_HAPPEN == false){
                alpha = 1;
                hostSetInfoVecs.triangles2Nodes_1[H0] = HEAD;
                hostSetInfoVecs.triangles2Nodes_2[H0] = edge_start;
                hostSetInfoVecs.triangles2Nodes_3[H0] = TAIL;
                hostSetInfoVecs.triangles2Nodes_1[T0] = HEAD;
                hostSetInfoVecs.triangles2Nodes_2[T0] = TAIL;
                hostSetInfoVecs.triangles2Nodes_3[T0] = edge_end;
                hostSetInfoVecs.triangles2Edges_1[H0] = iedge;
                hostSetInfoVecs.triangles2Edges_2[H0] = H1;
                hostSetInfoVecs.triangles2Edges_3[H0] = T1;
                hostSetInfoVecs.triangles2Edges_1[T0] = iedge;
                hostSetInfoVecs.triangles2Edges_2[T0] = T2;
                hostSetInfoVecs.triangles2Edges_3[T0] = H2;
                hostSetInfoVecs.edges2Nodes_1[iedge] = TAIL;
                hostSetInfoVecs.edges2Nodes_2[iedge] = HEAD;
                H1t1 = hostSetInfoVecs.edges2Triangles_1[H1];
                H1t2 = hostSetInfoVecs.edges2Triangles_2[H1]; //These are the associated triangles to edge H1
                if (H1t1 == H0){hostSetInfoVecs.edges2Triangles_1[H1] = H0;}
                if (H1t2 == H0){hostSetInfoVecs.edges2Triangles_2[H1] = H0;}
                H2t1 = hostSetInfoVecs.edges2Triangles_1[H2];
                H2t2 = hostSetInfoVecs.edges2Triangles_2[H2]; //These are the associated triangles to edge H2        
                if (H2t1 == H0){hostSetInfoVecs.edges2Triangles_1[H2] = T0;}
                if (H2t2 == H0){hostSetInfoVecs.edges2Triangles_2[H2] = T0;}
                T1t1 = hostSetInfoVecs.edges2Triangles_1[T1];
                T1t2 = hostSetInfoVecs.edges2Triangles_2[T1];
                if (T1t1 == T0){hostSetInfoVecs.edges2Triangles_1[T1] = H0;}
                if (T1t2 == T0){hostSetInfoVecs.edges2Triangles_2[T1] = H0;}
                T2t1 = hostSetInfoVecs.edges2Triangles_1[T2];
                T2t2 = hostSetInfoVecs.edges2Triangles_2[T2];
                if (T2t1 == T0){hostSetInfoVecs.edges2Triangles_1[T2] = T0;}
                if (T2t2 == T0){hostSetInfoVecs.edges2Triangles_2[T2] = T0;}
                ////std::cout<<"IS THE ERROR HERE 1?"<<std::endl;
                ////////////////////////////////////////////////////////////////////////////////////
                ///////////// UPDATING NEIGHBORING NODE INFO ///////////////////////////////////////
                ////////////////////////////////////////////////////////////////////////////////////

                ///////// DELETING CONNECTIVITY BETWEEN EDGE_START AND EDGE_END ////////////////////
                //int data_id;
                if (hostSetInfoVecs.nndata1[edge_start] == edge_end){
                    //data_id = 0;
                    hostSetInfoVecs.nndata1[edge_start] = -2;
                }
                else if (hostSetInfoVecs.nndata2[edge_start] == edge_end){
                    //data_id = 1;
                    hostSetInfoVecs.nndata2[edge_start] = -2;
                }
                else if (hostSetInfoVecs.nndata3[edge_start] == edge_end){
                //    data_id = 2;
                    hostSetInfoVecs.nndata3[edge_start] = -2;
                }
                else if (hostSetInfoVecs.nndata4[edge_start] == edge_end){
                //  data_id = 3;
                    hostSetInfoVecs.nndata4[edge_start] = -2;
                }
                else if (hostSetInfoVecs.nndata5[edge_start] == edge_end){
                    //data_id = 4;
                    hostSetInfoVecs.nndata5[edge_start] = -2;
                }
                else if (hostSetInfoVecs.nndata6[edge_start] == edge_end){
                //    data_id = 5;
                    hostSetInfoVecs.nndata6[edge_start] = -2;
                }
                else if (hostSetInfoVecs.nndata7[edge_start] == edge_end){
                //    data_id = 6;
                    hostSetInfoVecs.nndata7[edge_start] = -2;
                }
                else if (hostSetInfoVecs.nndata8[edge_start] == edge_end){
                //   data_id = 7;
                    hostSetInfoVecs.nndata8[edge_start] = -2;
                }
                else if (hostSetInfoVecs.nndata9[edge_start] == edge_end){
                // data_id = 8;
                    hostSetInfoVecs.nndata9[edge_start] = -2;
                }
                /*else if (hostSetInfoVecs.nndata10[edge_start] == edge_end){
                //   data_id = 9;
                    hostSetInfoVecs.nndata10[edge_start] = -2;
                }
                 else if (hostSetInfoVecs.nndata11[edge_start] == edge_end){
                // data_id = 10;
                    hostSetInfoVecs.nndata11[edge_start] = -2;
                }
                else if (hostSetInfoVecs.nndata12[edge_start] == edge_end){
                //   data_id = 11;
                    hostSetInfoVecs.nndata12[edge_start] = -2;
                } */
                else {}

                if (hostSetInfoVecs.nndata1[edge_end] == edge_start){
                // data_id = 0;
                    hostSetInfoVecs.nndata1[edge_end] = -2;
                }
                else if (hostSetInfoVecs.nndata2[edge_end] == edge_start){
                //  data_id = 1;
                    hostSetInfoVecs.nndata2[edge_end] = -2;
                }
                else if (hostSetInfoVecs.nndata3[edge_end] == edge_start){
                //  data_id = 2;
                    hostSetInfoVecs.nndata3[edge_end] = -2;
                }
                else if (hostSetInfoVecs.nndata4[edge_end] == edge_start){
                //  data_id = 3;
                    hostSetInfoVecs.nndata4[edge_end] = -2;
                }
                else if (hostSetInfoVecs.nndata5[edge_end] == edge_start){
                //  data_id = 4;
                    hostSetInfoVecs.nndata5[edge_end] = -2;
                }
                else if (hostSetInfoVecs.nndata6[edge_end] == edge_start){
                //  data_id = 5;
                    hostSetInfoVecs.nndata6[edge_end] = -2;
                }
                else if (hostSetInfoVecs.nndata7[edge_end] == edge_start){
                //  data_id = 6;
                    hostSetInfoVecs.nndata7[edge_end] = -2;
                }
                else if (hostSetInfoVecs.nndata8[edge_end] == edge_start){
                // data_id = 7;
                    hostSetInfoVecs.nndata8[edge_end] = -2;
                }
                else if (hostSetInfoVecs.nndata9[edge_end] == edge_start){
                // data_id = 8;
                    hostSetInfoVecs.nndata9[edge_end] = -2;
                }
                /*else if (hostSetInfoVecs.nndata10[edge_end] == edge_start){
                //  data_id = 9;
                    hostSetInfoVecs.nndata10[edge_end] = -2;
                }
                 else if (hostSetInfoVecs.nndata11[edge_end] == edge_start){
                //   data_id = 10;
                    hostSetInfoVecs.nndata11[edge_end] = -2;
                }
                else if (hostSetInfoVecs.nndata12[edge_end] == edge_start){
                //   data_id = 11;
                    hostSetInfoVecs.nndata12[edge_end] = -2;
                } */
                else {}
////std::cout<<"IS THE ERROR HERE 2?"<<std::endl;

            
                ///////////// ESTABLISHING NEW CONNECTIVITY ////////////////////
            if (hostSetInfoVecs.nndata1[HEAD] < 0){
                //   data_id = 0;
                    hostSetInfoVecs.nndata1[HEAD] = TAIL;
                }
                else if (hostSetInfoVecs.nndata2[HEAD] < 0){
                //   data_id = 1;
                    hostSetInfoVecs.nndata2[HEAD] = TAIL;
                }
                else if (hostSetInfoVecs.nndata3[HEAD] < 0){
                //   data_id = 2;
                    hostSetInfoVecs.nndata3[HEAD] = TAIL;
                }
                else if (hostSetInfoVecs.nndata4[HEAD] < 0){
                //  data_id = 3;
                    hostSetInfoVecs.nndata4[HEAD] = TAIL;
                }
                else if (hostSetInfoVecs.nndata5[HEAD] < 0){
                //   data_id = 4;
                    hostSetInfoVecs.nndata5[HEAD] = TAIL;
                }
                else if (hostSetInfoVecs.nndata6[HEAD] < 0){
                //   data_id = 5;
                    hostSetInfoVecs.nndata6[HEAD] = TAIL;
                }
                else if (hostSetInfoVecs.nndata7[HEAD] < 0){
                //  data_id = 6;
                    hostSetInfoVecs.nndata7[HEAD] = TAIL;
                }
                else if (hostSetInfoVecs.nndata8[HEAD] < 0){
                //  data_id = 7;
                    hostSetInfoVecs.nndata8[HEAD] = TAIL;
                }
                else if (hostSetInfoVecs.nndata9[HEAD] < 0){
                //  data_id = 8;
                    hostSetInfoVecs.nndata9[HEAD] = TAIL;
                }
                /*else if (hostSetInfoVecs.nndata10[HEAD] < 0){
                //  data_id = 9;
                    hostSetInfoVecs.nndata10[HEAD] = TAIL;
                }
                 else if (hostSetInfoVecs.nndata11[HEAD] < 0){
                //  data_id = 10;
                    hostSetInfoVecs.nndata11[HEAD] = TAIL;
                }
                else if (hostSetInfoVecs.nndata12[HEAD] < 0){
                //  data_id = 11;
                    hostSetInfoVecs.nndata12[HEAD] = TAIL;
                } */
                else {}

                ////std::cout<<"IS THE ERROR HERE 3?"<<std::endl;

                if (hostSetInfoVecs.nndata1[TAIL] < 0){
                //  data_id = 0;
                    hostSetInfoVecs.nndata1[TAIL] = HEAD;
                }
                else if (hostSetInfoVecs.nndata2[TAIL] < 0){
                //  data_id = 1;
                    hostSetInfoVecs.nndata2[TAIL] = HEAD;
                }
                else if (hostSetInfoVecs.nndata3[TAIL] < 0){
                //  data_id = 2;
                    hostSetInfoVecs.nndata3[TAIL] = HEAD;
                }
                else if (hostSetInfoVecs.nndata4[TAIL] < 0){
                //  data_id = 3;
                    hostSetInfoVecs.nndata4[TAIL] = HEAD;
                }
                else if (hostSetInfoVecs.nndata5[TAIL] < 0){
                //  data_id = 4;
                    hostSetInfoVecs.nndata5[TAIL] = HEAD;
                }
                else if (hostSetInfoVecs.nndata6[TAIL] < 0){
                //  data_id = 5;
                    hostSetInfoVecs.nndata6[TAIL] = HEAD;
                }
                else if (hostSetInfoVecs.nndata7[TAIL] < 0){
                //  data_id = 6;
                    hostSetInfoVecs.nndata7[TAIL] = HEAD;
                }
                else if (hostSetInfoVecs.nndata8[TAIL] < 0){
                //  data_id = 7;
                    hostSetInfoVecs.nndata8[TAIL] = HEAD;
                }
                else if (hostSetInfoVecs.nndata9[TAIL] < 0){
                //  data_id = 8;
                    hostSetInfoVecs.nndata9[TAIL] = HEAD;
                }
                /*else if (hostSetInfoVecs.nndata10[TAIL] < 0){
                //  data_id = 9;
                    hostSetInfoVecs.nndata10[TAIL] = HEAD;
                }
                 else if (hostSetInfoVecs.nndata11[TAIL] < 0){
                //  data_id = 10;
                    hostSetInfoVecs.nndata11[TAIL] = HEAD;
                }
                else if (hostSetInfoVecs.nndata12[TAIL] < 0){
                //  data_id = 11;
                    hostSetInfoVecs.nndata12[TAIL] = HEAD;
                } */
                else {}
                ////std::cout<<"IS THE ERROR HERE 4?"<<std::endl;

                //nndata[HEAD] += 1;
                //nndata[TAIL] += 1;
                //nndata[edge_start] -= 1;
                //nndata[edge_end] -= 1;
                
                for (int j = 0; j < 2; j++){
                    int nn, nnn, nnnn;
                    if (j == 0){
                        nn = edge_start;
                        nnnn = T0;
                    }
                    else if (j == 1){
                        nn = edge_end;
                        nnnn = H0;
                    }

                    if (hostSetInfoVecs.nodes2Triangles_1[nn] == nnnn){
                        hostSetInfoVecs.nodes2Triangles_1[nn] = -INT_MAX;
                    }
                    else if ( hostSetInfoVecs.nodes2Triangles_2[nn] == nnnn){
                        hostSetInfoVecs.nodes2Triangles_2[nn] = -INT_MAX;
                    }
                    else if (hostSetInfoVecs.nodes2Triangles_3[nn] == nnnn){
                        hostSetInfoVecs.nodes2Triangles_3[nn] = -INT_MAX;
                    }
                    else if (hostSetInfoVecs.nodes2Triangles_4[nn] == nnnn){
                        hostSetInfoVecs.nodes2Triangles_4[nn] = -INT_MAX;
                    }
                    else if (hostSetInfoVecs.nodes2Triangles_5[nn] == nnnn){
                        hostSetInfoVecs.nodes2Triangles_5[nn] = -INT_MAX;
                    }
                    else if ( hostSetInfoVecs.nodes2Triangles_6[nn] == nnnn){
                        hostSetInfoVecs.nodes2Triangles_6[nn] = -INT_MAX;
                    }
                    else if ( hostSetInfoVecs.nodes2Triangles_7[nn] == nnnn){
                        hostSetInfoVecs.nodes2Triangles_7[nn] = -INT_MAX;
                    }
                    else if (hostSetInfoVecs.nodes2Triangles_8[nn] == nnnn){
                        hostSetInfoVecs.nodes2Triangles_8[nn] = -INT_MAX;
                    }
                    else if (hostSetInfoVecs.nodes2Triangles_9[nn] == nnnn){
                        hostSetInfoVecs.nodes2Triangles_9[nn] = -INT_MAX;
                    }
                }

                for (int j = 0; j < 2; j++){
                    int nn, nnn, nnnn;
                    if (j == 0){
                        nn = HEAD;
                        nnnn = T0;
                    }
                    else if (j == 1){
                        nn = TAIL;
                        nnnn = H0;
                    }

                    if (hostSetInfoVecs.nodes2Triangles_1[nn] < 0){
                        hostSetInfoVecs.nodes2Triangles_1[nn] = nnnn;
                    }
                    else if ( hostSetInfoVecs.nodes2Triangles_2[nn] < 0){
                        hostSetInfoVecs.nodes2Triangles_2[nn] = nnnn;
                    }
                    else if (hostSetInfoVecs.nodes2Triangles_3[nn] < 0){
                        hostSetInfoVecs.nodes2Triangles_3[nn] = nnnn;
                    }
                    else if (hostSetInfoVecs.nodes2Triangles_4[nn] < 0){
                        hostSetInfoVecs.nodes2Triangles_4[nn] = nnnn;
                    }
                    else if (hostSetInfoVecs.nodes2Triangles_5[nn] < 0){
                        hostSetInfoVecs.nodes2Triangles_5[nn] = nnnn;
                    }
                    else if ( hostSetInfoVecs.nodes2Triangles_6[nn] < 0){
                        hostSetInfoVecs.nodes2Triangles_6[nn] = nnnn;
                    }
                    else if ( hostSetInfoVecs.nodes2Triangles_7[nn] < 0){
                        hostSetInfoVecs.nodes2Triangles_7[nn] = nnnn;
                    }
                    else if (hostSetInfoVecs.nodes2Triangles_8[nn] < 0){
                        hostSetInfoVecs.nodes2Triangles_8[nn] = nnnn;
                    }
                    else if (hostSetInfoVecs.nodes2Triangles_9[nn] < 0){
                        hostSetInfoVecs.nodes2Triangles_9[nn] = nnnn;
                    }
                }
        
            }
        } 
    }
    };  
    return alpha;
//This completes the update (if necessary) of the following data structures: triangles2Nodes, edges2Nodes, edges2Triangles.
};



//copy configuration from device to host
void Utilities::transferDtoH(GeneralParams& generalParams,
    CoordInfoVecs& coordInfoVecs,
    HostSetInfoVecs& hostSetInfoVecs){
    thrust::copy(coordInfoVecs.nodeLocX.begin(),coordInfoVecs.nodeLocX.end(),hostSetInfoVecs.nodeLocX.begin());
    thrust::copy(coordInfoVecs.nodeLocY.begin(),coordInfoVecs.nodeLocY.end(),hostSetInfoVecs.nodeLocY.begin());
    thrust::copy(coordInfoVecs.nodeLocZ.begin(),coordInfoVecs.nodeLocZ.end(),hostSetInfoVecs.nodeLocZ.begin());

    thrust::copy(generalParams.nodes_in_upperhem.begin(),generalParams.nodes_in_upperhem.end(),hostSetInfoVecs.nodes_in_upperhem.begin());
    thrust::copy(generalParams.triangles_in_upperhem.begin(),generalParams.triangles_in_upperhem.end(),hostSetInfoVecs.triangles_in_upperhem.begin());
    thrust::copy(generalParams.edges_in_upperhem.begin(),generalParams.edges_in_upperhem.end(),hostSetInfoVecs.edges_in_upperhem.begin());
    thrust::copy(generalParams.edges_in_upperhem_list.begin(),generalParams.edges_in_upperhem_list.end(),hostSetInfoVecs.edges_in_upperhem_list.begin());
    thrust::copy(generalParams.boundaries_in_upperhem.begin(),generalParams.boundaries_in_upperhem.end(),hostSetInfoVecs.boundaries_in_upperhem.begin());

    thrust::copy(coordInfoVecs.triangles2Nodes_1.begin(),coordInfoVecs.triangles2Nodes_1.end(),hostSetInfoVecs.triangles2Nodes_1.begin());
    thrust::copy(coordInfoVecs.triangles2Nodes_2.begin(),coordInfoVecs.triangles2Nodes_2.end(),hostSetInfoVecs.triangles2Nodes_2.begin());
    thrust::copy(coordInfoVecs.triangles2Nodes_3.begin(),coordInfoVecs.triangles2Nodes_3.end(),hostSetInfoVecs.triangles2Nodes_3.begin());
    
    thrust::copy(coordInfoVecs.edges2Nodes_1.begin(),coordInfoVecs.edges2Nodes_1.end(),hostSetInfoVecs.edges2Nodes_1.begin());
    thrust::copy(coordInfoVecs.edges2Nodes_2.begin(),coordInfoVecs.edges2Nodes_2.end(),hostSetInfoVecs.edges2Nodes_2.begin());
    
    thrust::copy(coordInfoVecs.edges2Triangles_1.begin(),coordInfoVecs.edges2Triangles_1.end(),hostSetInfoVecs.edges2Triangles_1.begin());
    thrust::copy(coordInfoVecs.edges2Triangles_2.begin(),coordInfoVecs.edges2Triangles_2.end(),hostSetInfoVecs.edges2Triangles_2.begin());
    
    thrust::copy(coordInfoVecs.triangles2Edges_1.begin(),coordInfoVecs.triangles2Edges_1.end(),hostSetInfoVecs.triangles2Edges_1.begin());
    thrust::copy(coordInfoVecs.triangles2Edges_2.begin(),coordInfoVecs.triangles2Edges_2.end(),hostSetInfoVecs.triangles2Edges_2.begin());
    thrust::copy(coordInfoVecs.triangles2Edges_3.begin(),coordInfoVecs.triangles2Edges_3.end(),hostSetInfoVecs.triangles2Edges_3.begin());

    thrust::copy(coordInfoVecs.nndata1.begin(),coordInfoVecs.nndata1.end(),hostSetInfoVecs.nndata1.begin());
    thrust::copy(coordInfoVecs.nndata2.begin(),coordInfoVecs.nndata2.end(),hostSetInfoVecs.nndata2.begin());
    thrust::copy(coordInfoVecs.nndata3.begin(),coordInfoVecs.nndata3.end(),hostSetInfoVecs.nndata3.begin());
    thrust::copy(coordInfoVecs.nndata4.begin(),coordInfoVecs.nndata4.end(),hostSetInfoVecs.nndata4.begin());
    thrust::copy(coordInfoVecs.nndata5.begin(),coordInfoVecs.nndata5.end(),hostSetInfoVecs.nndata5.begin());
    thrust::copy(coordInfoVecs.nndata6.begin(),coordInfoVecs.nndata6.end(),hostSetInfoVecs.nndata6.begin());
    thrust::copy(coordInfoVecs.nndata7.begin(),coordInfoVecs.nndata7.end(),hostSetInfoVecs.nndata7.begin());
    thrust::copy(coordInfoVecs.nndata8.begin(),coordInfoVecs.nndata8.end(),hostSetInfoVecs.nndata8.begin());
    thrust::copy(coordInfoVecs.nndata9.begin(),coordInfoVecs.nndata9.end(),hostSetInfoVecs.nndata9.begin());
    //thrust::copy(coordInfoVecs.nndata10.begin(),coordInfoVecs.nndata10.end(),hostSetInfoVecs.nndata10.begin());
    //thrust::copy(coordInfoVecs.nndata11.begin(),coordInfoVecs.nndata11.end(),hostSetInfoVecs.nndata11.begin());
    //thrust::copy(coordInfoVecs.nndata12.begin(),coordInfoVecs.nndata12.end(),hostSetInfoVecs.nndata12.begin());

    thrust::copy(coordInfoVecs.nodes2Triangles_1.begin(),coordInfoVecs.nodes2Triangles_1.end(),hostSetInfoVecs.nodes2Triangles_1.begin());
    thrust::copy(coordInfoVecs.nodes2Triangles_2.begin(),coordInfoVecs.nodes2Triangles_2.end(),hostSetInfoVecs.nodes2Triangles_2.begin());
    thrust::copy(coordInfoVecs.nodes2Triangles_3.begin(),coordInfoVecs.nodes2Triangles_3.end(),hostSetInfoVecs.nodes2Triangles_3.begin());
    thrust::copy(coordInfoVecs.nodes2Triangles_4.begin(),coordInfoVecs.nodes2Triangles_4.end(),hostSetInfoVecs.nodes2Triangles_4.begin());
    thrust::copy(coordInfoVecs.nodes2Triangles_5.begin(),coordInfoVecs.nodes2Triangles_5.end(),hostSetInfoVecs.nodes2Triangles_5.begin());
    thrust::copy(coordInfoVecs.nodes2Triangles_6.begin(),coordInfoVecs.nodes2Triangles_6.end(),hostSetInfoVecs.nodes2Triangles_6.begin());
    thrust::copy(coordInfoVecs.nodes2Triangles_7.begin(),coordInfoVecs.nodes2Triangles_7.end(),hostSetInfoVecs.nodes2Triangles_7.begin());
    thrust::copy(coordInfoVecs.nodes2Triangles_8.begin(),coordInfoVecs.nodes2Triangles_8.end(),hostSetInfoVecs.nodes2Triangles_8.begin());
    thrust::copy(coordInfoVecs.nodes2Triangles_9.begin(),coordInfoVecs.nodes2Triangles_9.end(),hostSetInfoVecs.nodes2Triangles_9.begin());

    thrust::copy(coordInfoVecs.triangles2Triangles_1.begin(),coordInfoVecs.triangles2Triangles_1.end(),hostSetInfoVecs.triangles2Triangles_1.begin());
    thrust::copy(coordInfoVecs.triangles2Triangles_2.begin(),coordInfoVecs.triangles2Triangles_2.end(),hostSetInfoVecs.triangles2Triangles_2.begin());
    thrust::copy(coordInfoVecs.triangles2Triangles_3.begin(),coordInfoVecs.triangles2Triangles_3.end(),hostSetInfoVecs.triangles2Triangles_3.begin());
};

//copy configuration from host to device
void Utilities::transferHtoD(GeneralParams& generalParams,
    CoordInfoVecs& coordInfoVecs,
    HostSetInfoVecs& hostSetInfoVecs){
    thrust::copy(hostSetInfoVecs.nodeLocX.begin(),hostSetInfoVecs.nodeLocX.end(),coordInfoVecs.nodeLocX.begin());
    thrust::copy(hostSetInfoVecs.nodeLocY.begin(),hostSetInfoVecs.nodeLocY.end(),coordInfoVecs.nodeLocY.begin());
    thrust::copy(hostSetInfoVecs.nodeLocZ.begin(),hostSetInfoVecs.nodeLocZ.end(),coordInfoVecs.nodeLocZ.begin());

    thrust::copy(hostSetInfoVecs.nodes_in_upperhem.begin(),hostSetInfoVecs.nodes_in_upperhem.end(),generalParams.nodes_in_upperhem.begin());
    thrust::copy(hostSetInfoVecs.triangles_in_upperhem.begin(),hostSetInfoVecs.triangles_in_upperhem.end(),generalParams.triangles_in_upperhem.begin());
    thrust::copy(hostSetInfoVecs.edges_in_upperhem.begin(),hostSetInfoVecs.edges_in_upperhem.end(),generalParams.edges_in_upperhem.begin());
    thrust::copy(hostSetInfoVecs.edges_in_upperhem_list.begin(),hostSetInfoVecs.edges_in_upperhem_list.end(),generalParams.edges_in_upperhem_list.begin());
    thrust::copy(hostSetInfoVecs.boundaries_in_upperhem.begin(),hostSetInfoVecs.boundaries_in_upperhem.end(),generalParams.boundaries_in_upperhem.begin());
    
    thrust::copy(hostSetInfoVecs.triangles2Nodes_1.begin(),hostSetInfoVecs.triangles2Nodes_1.end(),coordInfoVecs.triangles2Nodes_1.begin());
    thrust::copy(hostSetInfoVecs.triangles2Nodes_2.begin(),hostSetInfoVecs.triangles2Nodes_2.end(),coordInfoVecs.triangles2Nodes_2.begin());
    thrust::copy(hostSetInfoVecs.triangles2Nodes_3.begin(),hostSetInfoVecs.triangles2Nodes_3.end(),coordInfoVecs.triangles2Nodes_3.begin());
    
    thrust::copy(hostSetInfoVecs.edges2Nodes_1.begin(),hostSetInfoVecs.edges2Nodes_1.end(),coordInfoVecs.edges2Nodes_1.begin());
    thrust::copy(hostSetInfoVecs.edges2Nodes_2.begin(),hostSetInfoVecs.edges2Nodes_2.end(),coordInfoVecs.edges2Nodes_2.begin());
    
    thrust::copy(hostSetInfoVecs.edges2Triangles_1.begin(),hostSetInfoVecs.edges2Triangles_1.end(),coordInfoVecs.edges2Triangles_1.begin());
    thrust::copy(hostSetInfoVecs.edges2Triangles_2.begin(),hostSetInfoVecs.edges2Triangles_2.end(),coordInfoVecs.edges2Triangles_2.begin());
    
    thrust::copy(hostSetInfoVecs.triangles2Edges_1.begin(),hostSetInfoVecs.triangles2Edges_1.end(),coordInfoVecs.triangles2Edges_1.begin());
    thrust::copy(hostSetInfoVecs.triangles2Edges_2.begin(),hostSetInfoVecs.triangles2Edges_2.end(),coordInfoVecs.triangles2Edges_2.begin());
    thrust::copy(hostSetInfoVecs.triangles2Edges_3.begin(),hostSetInfoVecs.triangles2Edges_3.end(),coordInfoVecs.triangles2Edges_3.begin());

    thrust::copy(hostSetInfoVecs.nndata1.begin(),hostSetInfoVecs.nndata1.end(),coordInfoVecs.nndata1.begin());
    thrust::copy(hostSetInfoVecs.nndata2.begin(),hostSetInfoVecs.nndata2.end(),coordInfoVecs.nndata2.begin());
    thrust::copy(hostSetInfoVecs.nndata3.begin(),hostSetInfoVecs.nndata3.end(),coordInfoVecs.nndata3.begin());
    thrust::copy(hostSetInfoVecs.nndata4.begin(),hostSetInfoVecs.nndata4.end(),coordInfoVecs.nndata4.begin());
    thrust::copy(hostSetInfoVecs.nndata5.begin(),hostSetInfoVecs.nndata5.end(),coordInfoVecs.nndata5.begin());
    thrust::copy(hostSetInfoVecs.nndata6.begin(),hostSetInfoVecs.nndata6.end(),coordInfoVecs.nndata6.begin());
    thrust::copy(hostSetInfoVecs.nndata7.begin(),hostSetInfoVecs.nndata7.end(),coordInfoVecs.nndata7.begin());
    thrust::copy(hostSetInfoVecs.nndata8.begin(),hostSetInfoVecs.nndata8.end(),coordInfoVecs.nndata8.begin());
    thrust::copy(hostSetInfoVecs.nndata9.begin(),hostSetInfoVecs.nndata9.end(),coordInfoVecs.nndata9.begin());
    //thrust::copy(hostSetInfoVecs.nndata10.begin(),hostSetInfoVecs.nndata10.end(),coordInfoVecs.nndata10.begin());
    //thrust::copy(hostSetInfoVecs.nndata11.begin(),hostSetInfoVecs.nndata11.end(),coordInfoVecs.nndata11.begin());
    //thrust::copy(hostSetInfoVecs.nndata12.begin(),hostSetInfoVecs.nndata12.end(),coordInfoVecs.nndata12.begin());

    thrust::copy(hostSetInfoVecs.nodes2Triangles_1.begin(),hostSetInfoVecs.nodes2Triangles_1.end(),coordInfoVecs.nodes2Triangles_1.begin());
    thrust::copy(hostSetInfoVecs.nodes2Triangles_2.begin(),hostSetInfoVecs.nodes2Triangles_2.end(),coordInfoVecs.nodes2Triangles_2.begin());
    thrust::copy(hostSetInfoVecs.nodes2Triangles_3.begin(),hostSetInfoVecs.nodes2Triangles_3.end(),coordInfoVecs.nodes2Triangles_3.begin());
    thrust::copy(hostSetInfoVecs.nodes2Triangles_4.begin(),hostSetInfoVecs.nodes2Triangles_4.end(),coordInfoVecs.nodes2Triangles_4.begin());
    thrust::copy(hostSetInfoVecs.nodes2Triangles_5.begin(),hostSetInfoVecs.nodes2Triangles_5.end(),coordInfoVecs.nodes2Triangles_5.begin());
    thrust::copy(hostSetInfoVecs.nodes2Triangles_6.begin(),hostSetInfoVecs.nodes2Triangles_6.end(),coordInfoVecs.nodes2Triangles_6.begin());
    thrust::copy(hostSetInfoVecs.nodes2Triangles_7.begin(),hostSetInfoVecs.nodes2Triangles_7.end(),coordInfoVecs.nodes2Triangles_7.begin());
    thrust::copy(hostSetInfoVecs.nodes2Triangles_8.begin(),hostSetInfoVecs.nodes2Triangles_8.end(),coordInfoVecs.nodes2Triangles_8.begin());
    thrust::copy(hostSetInfoVecs.nodes2Triangles_9.begin(),hostSetInfoVecs.nodes2Triangles_9.end(),coordInfoVecs.nodes2Triangles_9.begin());

    thrust::copy(hostSetInfoVecs.triangles2Triangles_1.begin(),hostSetInfoVecs.triangles2Triangles_1.end(),coordInfoVecs.triangles2Triangles_1.begin());
    thrust::copy(hostSetInfoVecs.triangles2Triangles_2.begin(),hostSetInfoVecs.triangles2Triangles_2.end(),coordInfoVecs.triangles2Triangles_2.begin());
    thrust::copy(hostSetInfoVecs.triangles2Triangles_3.begin(),hostSetInfoVecs.triangles2Triangles_3.end(),coordInfoVecs.triangles2Triangles_3.begin());
};

/* Make_ele_structure_PDE build the data structure necessary for the PDE computation.
    Not all elements of the structure are introduced as some require another function to 
    compute. 
    Contents not updated here: side_vector, side_norm, aff_conorm, aff_adj_ele_conorm,
                                affine_pt, aff_cent */
/*void Utilities::make_ele_structure_PDE(ELEMENT element,
    int id,
    CoordInfoVecs& coordInfoVecs,
    AreaTriangleInfoVecs& areaTriangleInfoVecs){
        int i = id;
        element.Id = id;
        int vt1 = coordInfoVecs.triangles2Nodes_1[i];
        int vt2 = coordInfoVecs.triangles2Nodes_2[i];
        int vt3 = coordInfoVecs.triangles2Nodes_3[i];
        int e1 = coordInfoVecs.triangles2Edges_1[i];
        int e2 = coordInfoVecs.triangles2Edges_2[i];
        int e3 = coordInfoVecs.triangles2Edges_3[i];
        std::vector<int> e1v = {coordInfoVecs.edges2Nodes_1[e1], coordInfoVecs.edges2Nodes_2[e1]};
        std::vector<int> e2v = {coordInfoVecs.edges2Nodes_1[e2], coordInfoVecs.edges2Nodes_2[e2]};
        std::vector<int> e3v = {coordInfoVecs.edges2Nodes_1[e3], coordInfoVecs.edges2Nodes_2[e3]};
        element.vertex[0] = vt1;
        element.vertex[1] = vt2;
        element.vertex[2] = vt3;
        element.adj_ele[0] = coordInfoVecs.triangles2Triangles_1[i];
        element.adj_ele[1] = coordInfoVecs.triangles2Triangles_2[i];
        element.adj_ele[2] = coordInfoVecs.triangles2Triangles_3[i];
        element.area0 = areaTriangleInfoVecs.initial_area;
        element.center[0] = (coordInfoVecs.nodeLocX[vt1] + coordInfoVecs.nodeLocX[vt2] + coordInfoVecs.nodeLocX[vt3])/3.0;
        element.center[1] = (coordInfoVecs.nodeLocY[vt1] + coordInfoVecs.nodeLocY[vt2] + coordInfoVecs.nodeLocY[vt3])/3.0;
        element.center[2] = (coordInfoVecs.nodeLocZ[vt1] + coordInfoVecs.nodeLocZ[vt2] + coordInfoVecs.nodeLocZ[vt3])/3.0;
        element.length_side[0] = sqrt((coordInfoVecs.nodeLocX[e1v[1]] - coordInfoVecs.nodeLocX[e1v[0]])*(coordInfoVecs.nodeLocX[e1v[1]] - coordInfoVecs.nodeLocX[e1v[0]]) +
                                    (coordInfoVecs.nodeLocY[e1v[1]] - coordInfoVecs.nodeLocY[e1v[0]])*(coordInfoVecs.nodeLocY[e1v[1]] - coordInfoVecs.nodeLocY[e1v[0]]) +
                                    (coordInfoVecs.nodeLocZ[e1v[1]] - coordInfoVecs.nodeLocZ[e1v[0]])*(coordInfoVecs.nodeLocZ[e1v[1]] - coordInfoVecs.nodeLocZ[e1v[0]]));
        //std::cout<<"length_side[0] = "<<element.length_side[0]<<std::endl;
        element.length_side[1] = sqrt((coordInfoVecs.nodeLocX[e2v[1]] - coordInfoVecs.nodeLocX[e2v[0]])*(coordInfoVecs.nodeLocX[e2v[1]] - coordInfoVecs.nodeLocX[e2v[0]]) +
                                    (coordInfoVecs.nodeLocY[e2v[1]] - coordInfoVecs.nodeLocY[e2v[0]])*(coordInfoVecs.nodeLocY[e2v[1]] - coordInfoVecs.nodeLocY[e2v[0]]) +
                                    (coordInfoVecs.nodeLocZ[e2v[1]] - coordInfoVecs.nodeLocZ[e2v[0]])*(coordInfoVecs.nodeLocZ[e2v[1]] - coordInfoVecs.nodeLocZ[e2v[0]]));
        //std::cout<<"length_side[1] = "<<element.length_side[1]<<std::endl;
        element.length_side[2] = sqrt((coordInfoVecs.nodeLocX[e3v[1]] - coordInfoVecs.nodeLocX[e3v[0]])*(coordInfoVecs.nodeLocX[e3v[1]] - coordInfoVecs.nodeLocX[e3v[0]]) +
                                    (coordInfoVecs.nodeLocY[e3v[1]] - coordInfoVecs.nodeLocY[e3v[0]])*(coordInfoVecs.nodeLocY[e3v[1]] - coordInfoVecs.nodeLocY[e3v[0]]) +
                                    (coordInfoVecs.nodeLocZ[e3v[1]] - coordInfoVecs.nodeLocZ[e3v[0]])*(coordInfoVecs.nodeLocZ[e3v[1]] - coordInfoVecs.nodeLocZ[e3v[0]]));
        //std::cout<<"length_side[2] = "<<element.length_side[2]<<std::endl;
        double dummy_length = (element.length_side[0] + element.length_side[1] + element.length_side[2])/2.0;
        //std::cout<<"dummy_length = "<<dummy_length<<std::endl;
        
        
        element.area = sqrt(dummy_length*(-element.length_side[0] + dummy_length)*(-element.length_side[1] + dummy_length)*(-element.length_side[2] + dummy_length)) ;
        //std::cout<<"element.area ="<<element.area<<std::endl;
        element.out_norm.x[0] = (coordInfoVecs.nodeLocY[vt2] - coordInfoVecs.nodeLocY[vt1])*
                                    (coordInfoVecs.nodeLocZ[vt3] - coordInfoVecs.nodeLocZ[vt1]) - 
                                    (coordInfoVecs.nodeLocY[vt3] - coordInfoVecs.nodeLocY[vt3])*
                                    (coordInfoVecs.nodeLocZ[vt2] - coordInfoVecs.nodeLocZ[vt1]);
        //std::cout<<"out_norma.x[0] = "<<element.out_norm.x[0]<<std::endl;
        element.out_norm.x[1] = (coordInfoVecs.nodeLocX[vt2] - coordInfoVecs.nodeLocX[vt1])*
                                    (coordInfoVecs.nodeLocZ[vt3] - coordInfoVecs.nodeLocZ[vt1]) - 
                                    (coordInfoVecs.nodeLocX[vt3] - coordInfoVecs.nodeLocX[vt3])*
                                    (coordInfoVecs.nodeLocZ[vt2] - coordInfoVecs.nodeLocZ[vt1]);
        //std::cout<<"out_norma.x[0] = "<<element.out_norm.x[1]<<std::endl;
        element.out_norm.x[2] = (coordInfoVecs.nodeLocX[vt2] - coordInfoVecs.nodeLocX[vt1])*
                                    (coordInfoVecs.nodeLocY[vt3] - coordInfoVecs.nodeLocY[vt1]) - 
                                    (coordInfoVecs.nodeLocX[vt3] - coordInfoVecs.nodeLocX[vt3])*
                                    (coordInfoVecs.nodeLocY[vt2] - coordInfoVecs.nodeLocY[vt1]);
        //std::cout<<"out_norma.x[0] = "<<element.out_norm.x[2]<<std::endl;
        element.vertex[0] = coordInfoVecs.triangles2Nodes_1[i];
        element.vertex[1] = coordInfoVecs.triangles2Nodes_2[i];
        element.vertex[2] = coordInfoVecs.triangles2Nodes_3[i];
        element.affine_trans = NULL;
        element.inv_affine_trans = NULL;
    }*/


void Utilities::make_rbc_structure(RBC *rbc, CoordInfoVecs& coordInfoVecs, GeneralParams& generalParams){
        int      /*i,j,k,*/index,v1,v2,v3,v4,v5, tmp_count;
        int      testforread,restart=0;
        double   dtemp, tmpx, tmpy,tmpz, x0 = 0.0;
        CVector  vt1,vt2,vt[3], tmpcent;
        int      bFind,n,n1,n2,n3,n4,n5,adjtri;
        double   rotA[2][2];
       
        //std::list<int>::const_iterator it;
        std::list<int>  tmp_eles_at_node;
        int        nghbr_ele, nghbr_side;
        double maxx, maxy, maxz, minx, miny, minz;
        double lxcritical, rxcritical; 
        int    Nmaxx, Nminx, Nmaxy,  Nminy, Nmaxz,  Nminz;

		


	//using namespace std;
	// rbc->EleNode=new std::list<int>[rbc->NNum];

        //Loading mesh-------------------------------------
		rbc->avg_radi = 0.0;
        for(int i=0;i<rbc->NNum;i++)
        {
			rbc->node[i].x[0] = coordInfoVecs.nodeLocX[i];
			rbc->node[i].x[1] = coordInfoVecs.nodeLocY[i];
			rbc->node[i].x[2] = coordInfoVecs.nodeLocZ[i];
            tmpcent.x[0] += coordInfoVecs.nodeLocX[i];
            tmpcent.x[1] += coordInfoVecs.nodeLocY[i];
            tmpcent.x[2] += coordInfoVecs.nodeLocZ[i];
            // scan x-, y-, z-coords of nodes
            rbc->nw[i]=0.0;
        }

        tmpcent.x[0] /= rbc->NNum; 
        tmpcent.x[1] /= rbc->NNum; 
        tmpcent.x[2] /= rbc->NNum;
        rbc->Center.x[0] = tmpcent.x[0];
        rbc->Center.x[1] = tmpcent.x[1];
        rbc->Center.x[2] = tmpcent.x[2];         
        for (int i = 0; i < rbc->NNum; i++){
            rbc->avg_radi += sqrt((rbc->node[i].x[0] - rbc->Center.x[0])*(rbc->node[i].x[0] - rbc->Center.x[0]) +
                                    (rbc->node[i].x[1] - rbc->Center.x[1])*(rbc->node[i].x[1] - rbc->Center.x[1]) +
                                    (rbc->node[i].x[2] - rbc->Center.x[2])*(rbc->node[i].x[2] - rbc->Center.x[2]));
        }
        rbc->avg_radi = rbc->avg_radi/rbc->NNum;

        //// TMP
        {
                maxx=minx=rbc->node[0].x[0];
                maxy=miny=rbc->node[0].x[1];
                maxz=minz=rbc->node[0].x[2];
                Nmaxx = Nminx =0;
                Nmaxy = Nminy =0;
                Nmaxz = Nminz =0;

                for(int i=1;i<rbc->NNum;i++)
                {
                    if(rbc->node[i].x[0]>maxx){
                        maxx=rbc->node[i].x[0]; Nmaxx = i;
                    }
                    if(rbc->node[i].x[0]<minx){
                        minx=rbc->node[i].x[0]; Nminx = i;
                    }

                    if(rbc->node[i].x[1]>maxy){
                        maxy=rbc->node[i].x[1]; Nmaxy = i;
                    }
                    if(rbc->node[i].x[1]<miny){
                        miny=rbc->node[i].x[1]; Nminy = i;
                    }

                    if(rbc->node[i].x[2]>maxz){
                        maxz=rbc->node[i].x[2]; Nmaxz = i;
                    }
                    if(rbc->node[i].x[2]<minz){
                        minz=rbc->node[i].x[2]; Nminz = i;
                    }
                }

                // printf("LoadMesh() cell domain x[%g, %g], y[%g, %g], z[%g, %g]\n", 
                //            minx, maxx, miny, maxy,minz, maxz);
        }
        ////END::: TMP
		
        /*The surface mesh "sphere.neu for example" in current neu file has inward normal*/
        /*So exchange orders of node indices when read in*/
        rbc->Curr_it = rbc->Curr_ele->begin();
        for(int i=0; (rbc->Curr_it) != rbc->Curr_ele->end();i++, (rbc->Curr_it)++)
        {
            if (i >= rbc->ENum){
                break;
            }
            rbc->Curr_it->vertex[0] = coordInfoVecs.triangles2Nodes_1[i];
			rbc->Curr_it->vertex[1] = coordInfoVecs.triangles2Nodes_2[i];
			rbc->Curr_it->vertex[2] = coordInfoVecs.triangles2Nodes_3[i];
			// rbc->ele[i].vertex[0] = coordInfoVecs.triangles2Nodes_1[i];
			// rbc->ele[i].vertex[1] = coordInfoVecs.triangles2Nodes_2[i];
			// rbc->ele[i].vertex[2] = coordInfoVecs.triangles2Nodes_3[i];
            
            rbc->Curr_it->adj_ele[0] = coordInfoVecs.triangles2Triangles_1[i];
			rbc->Curr_it->adj_ele[1] = coordInfoVecs.triangles2Triangles_2[i];
			rbc->Curr_it->adj_ele[2] = coordInfoVecs.triangles2Triangles_3[i];

            rbc->Curr_it->u_for_dg[0] = coordInfoVecs.u[i];

		// 	int edge;
		// 	for (int j = 0; j < 3; j++){
		// 		if (j == 0){
		// 			edge = coordInfoVecs.triangles2Edges_1[i];
		// 		}
		// 		else if (j == 1){
		// 			edge = coordInfoVecs.triangles2Edges_2[i];
		// 		}
		// 		else if (j == 2){
		// 			edge = coordInfoVecs.triangles2Edges_3[i];
		// 		}
		// 		int et1 = coordInfoVecs.edges2Triangles_1[edge];
		// 		int et2 = coordInfoVecs.edges2Triangles_2[edge];
		// 		if (coordInfoVecs.edges2Nodes_1[edge] == rbc->ele[i].vertex[0] && coordInfoVecs.edges2Nodes_2[edge] == rbc->ele[i].vertex[1]){
		// 			//note that here adj_ele[0] is the triangles sharing the edge consisted of vertex 0 and vertex 1;
		// 			if (et1 == i){
		// 				rbc->Curr_it->adj_ele[0] = et2;
        //                 // rbc->ele[i].adj_ele[0] = et2;
		// 			}
		// 			else if (et2 == i){
		// 				rbc->Curr_it->adj_ele[0] = et1;
        //                 // rbc->ele[i].adj_ele[0] = et1;
		// 			}
		// 		}
		// 		else if (coordInfoVecs.edges2Nodes_1[edge] == rbc->ele[i].vertex[1] && coordInfoVecs.edges2Nodes_2[edge] == rbc->ele[i].vertex[0]){
		// 			if (et1 == i){
        //                 rbc->Curr_it->adj_ele[0] = et2;
		// 				// rbc->ele[i].adj_ele[0] = et2;
		// 			}
		// 			else if (et2 == i){
        //                 rbc->Curr_it->adj_ele[0] = et1;
		// 				// rbc->ele[i].adj_ele[0] = et1;
		// 			}
		// 		}
		// 		if (coordInfoVecs.edges2Nodes_1[edge] == rbc->ele[i].vertex[1] && coordInfoVecs.edges2Nodes_2[edge] == rbc->ele[i].vertex[2]){
		// 			//note that here adj_ele[0] is the triangles sharing the edge consisted of vertex 0 and vertex 1;
		// 			if (et1 == i){
		// 				rbc->Curr_it->adj_ele[1] = et2;
        //                 // rbc->ele[i].adj_ele[1] = et2;
		// 			}
		// 			else if (et2 == i){
		// 				rbc->Curr_it->adj_ele[1] = et1;
		// 			}
		// 		}
		// 		else if (coordInfoVecs.edges2Nodes_1[edge] == rbc->ele[i].vertex[2] && coordInfoVecs.edges2Nodes_2[edge] == rbc->ele[i].vertex[1]){
		// 			if (et1 == i){
		// 				rbc->Curr_it->adj_ele[1] = et2;
		// 			}
		// 			else if (et2 == i){
		// 				rbc->Curr_it->adj_ele[1] = et1;
		// 			}
		// 		}
		// 		if (coordInfoVecs.edges2Nodes_1[edge] == rbc->ele[i].vertex[2] && coordInfoVecs.edges2Nodes_2[edge] == rbc->ele[i].vertex[0]){
		// 			//note that here adj_ele[0] is the triangles sharing the edge consisted of vertex 0 and vertex 1;
		// 			if (et1 == i){
		// 				rbc->Curr_it->adj_ele[2] = et2;
		// 			}
		// 			else if (et2 == i){
		// 				rbc->Curr_it->adj_ele[2] = et1;
		// 			}
		// 		}
		// 		else if (coordInfoVecs.edges2Nodes_1[edge] == rbc->ele[i].vertex[0] && coordInfoVecs.edges2Nodes_2[edge] == rbc->ele[i].vertex[2]){
		// 			if (et1 == i){
		// 				rbc->Curr_it->adj_ele[2] = et2;
		// 			}
		// 			else if (et2 == i){
		// 				rbc->Curr_it->adj_ele[2] = et1;
		// 			}
		// 		}
		// 	}
        //     // for(int j=0;j<3;j++)
        //     // {
        //     //     (rbc->EleNode[rbc->ele[i].vertex[j]]).push_back(i);  
        //     //           /// EleNode[].push_back --> This node belongs to which element. 
        //     //           /// TO the end, it has elements sharing the same node with node index rbc->ele[i].vertex[j].
        //     // }
        }
		
        
        rbc->Curr_it = rbc->Curr_ele->begin();
        bool deleted_elem;
        for(int i=0;i<rbc->ENum;i++, (rbc->Curr_it)++)
        {
            deleted_elem = false;
            // std::cout<<i<<std::endl;
            for(int j=0;j<3;j++)
            {
                if (rbc->Curr_it->vertex[j] >= (INT_MAX-100)){
                    deleted_elem = true;
                    break;
                }
                else{
                    vt[j]=rbc->node[rbc->Curr_it->vertex[j]];
                }
            }
        
            if (deleted_elem == true){
                continue;
            }
            rbc->Curr_it->out_norm=Cross(vt[1]-vt[0],vt[2]-vt[0]);
            //dtemp=Modul(rbc->ele[i].out_norm);
			dtemp = sqrt((rbc->Curr_it->out_norm.x[0])*(rbc->Curr_it->out_norm.x[0]) + (rbc->Curr_it->out_norm.x[1])*(rbc->Curr_it->out_norm.x[1]) + (rbc->Curr_it->out_norm.x[2])*(rbc->Curr_it->out_norm.x[2]));
			rbc->Curr_it->out_norm=rbc->Curr_it->out_norm/dtemp;
            rbc->Curr_it->area=0.5*dtemp;
            rbc->Curr_it->area0=rbc->Curr_it->area;
            rbc->Curr_it->Id = i; 
            for(int j=0;j<3;j++)
            {
                rbc->nw[rbc->Curr_it->vertex[j]]+=rbc->Curr_it->area/3.0;
			}
			rbc->Curr_it->center[0] = (1.0/3.0)*(vt[0].x[0] + vt[1].x[0] + vt[2].x[0]);
			//std::cout<<rbc->ele[i].center[0]<<std::endl;
			rbc->Curr_it->center[1] = (1.0/3.0)*(vt[0].x[1] + vt[1].x[1] + vt[2].x[1]);
			rbc->Curr_it->center[2] = (1.0/3.0)*(vt[0].x[2] + vt[1].x[2] + vt[2].x[2]);
            for(int j = 0; j < 3; j++){
                rbc->Curr_it->side_vector[j] = vt[(j+1)%3] - vt[j]; 
                rbc->Curr_it->length_side[j] = Modul(rbc->Curr_it->side_vector[j]);
                if (rbc->Curr_it->Id == 0 || rbc->Curr_it->Id == 960){
                    // std::cout<<"Elem["<<rbc->Curr_it->Id<<"], edge["<<j<<"], "<<rbc->Curr_it->length_side[j]<<std::endl;
                } 
            }

            for(int j = 0; j < 3; j++){
                rbc->Curr_it->side_norm[j] = 
                    Cross(rbc->Curr_it->side_vector[j], rbc->Curr_it->out_norm); 
                dtemp = Modul(rbc->Curr_it->side_norm[j]);
                rbc->Curr_it->side_norm[j] = rbc->Curr_it->side_norm[j]/dtemp;
            }
        }
        

        // for(int i=0;i<rbc->NNum;i++)
        //     rbc->EleNode[i].clear();

        rbc->LNum = coordInfoVecs.num_edges;

		

		if (rbc == NULL){
			std::cout<<"rbc is null"<<std::endl;
		}


	
	// for(int tmp_k=0;tmp_k<rbc->ENum;tmp_k++)
	// {
		
	// 	if (rbc->Curr_ele_index == NULL){
	// 		std::cout<<"rbc->Curr_ele_index is NULL"<<std::endl;
	// 	}
	// 	rbc->Curr_ele_index->push_back(tmp_k);
		
	// 	rbc->Curr_ele->push_back(rbc->ele[tmp_k]);
    // }
    
};

void Utilities::Tri_affine_trans(
	double   *pt0,
        double   *pt1,
        double   *pt2,
        double   *norm, 
        double   **affine_trans)
{
        int           i, j, k;
        double        *tri_norm, tmp_vec[4], tang_x[4], tang_y[4], len, tan_tmp[4];
        static double **rot_matr_rows = NULL, **translate;
        double        Eu_crds[3][3] = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
        
    //    std::cout<<"seg fault 1"<<std::endl;
        if(NULL == rot_matr_rows)
        {
            FT_matrix(&rot_matr_rows, 4, 4,sizeof(double));
            FT_matrix(&translate, 4, 4,sizeof(double));
            for(j = 0; j < 4; j++)
                for(k = 0; k < 4; k++)
                    translate[j][k] = 0.0;
        }
     //   std::cout<<"seg fault 2"<<std::endl;
        if(NULL != norm)
            tri_norm = norm;
        else
        {
            printf("ERROR: implement compute tri norm\n");
            //clean_up(ERROR); //20191006
        }
    //    std::cout<<"seg fault 3"<<std::endl;
        for(j=0;j<3;j++)
            tmp_vec[j] = (pt1[j] - pt0[j]);
        len = Mag3d(tmp_vec);
        for(j=0;j<3;j++) tang_x[j] = tmp_vec[j]/len;
        Cross3d(tri_norm, tang_x, tan_tmp);
        len = Mag3d(tan_tmp);
        for(j=0;j<3;j++) tang_y[j] = tan_tmp[j]/len;
  //  std::cout<<"seg fault 4"<<std::endl;
        rot_matr_rows[0][0] = Dot3d(tang_x, Eu_crds[0]);
        rot_matr_rows[0][1] = Dot3d(tang_x, Eu_crds[1]);
        rot_matr_rows[0][2] = Dot3d(tang_x, Eu_crds[2]);
        rot_matr_rows[0][3] = 0.0;
//std::cout<<"seg fault 5"<<std::endl;
        rot_matr_rows[1][0] = Dot3d(tang_y, Eu_crds[0]);
        rot_matr_rows[1][1] = Dot3d(tang_y, Eu_crds[1]);
        rot_matr_rows[1][2] = Dot3d(tang_y, Eu_crds[2]);
        rot_matr_rows[1][3] = 0.0;
//std::cout<<"seg fault 6"<<std::endl;
        rot_matr_rows[2][0] = Dot3d(tri_norm, Eu_crds[0]);
        rot_matr_rows[2][1] = Dot3d(tri_norm, Eu_crds[1]);
        rot_matr_rows[2][2] = Dot3d(tri_norm, Eu_crds[2]);
        rot_matr_rows[2][3] = 0.0;
//std::cout<<"seg fault 7"<<std::endl;
        rot_matr_rows[3][0] = 0.0;
        rot_matr_rows[3][1] = 0.0;
        rot_matr_rows[3][2] = 0.0;
        rot_matr_rows[3][3] = 1.0;
//std::cout<<"seg fault 8"<<std::endl;
        translate[0][0] = 1.0;
        translate[1][1] = 1.0;
        translate[2][2] = 1.0;
        translate[3][3] = 1.0;
        translate[0][3] = -pt0[0];
        translate[1][3] = -pt0[1];
        translate[2][3] = -pt0[2];
//std::cout<<"seg fault 9"<<std::endl;
        matrix_matrix_mult(rot_matr_rows, translate, 4, 4, affine_trans);
  //      std::cout<<"seg fault 10"<<std::endl;
}

void Utilities::grad_vh_loc_basis(
        double  *crds,
        double  *cent,
        double  sqrt_area,
        int     indx,
        double  *val)
{
        val[0] = val[1] = 0.0;
        switch(indx)
        {
        case 0:
        break;
        case 1:
            val[0] = 1.0/sqrt_area;
        break;
        case 2:
            val[1] = 1.0/sqrt_area;
        break;
        case 3:
            val[0] = 2.0*(crds[0]-cent[0])/sqr(sqrt_area);
        break;
        case 4:
            val[0] = (crds[1] - cent[1])/sqr(sqrt_area);
            val[1] = (crds[0] - cent[0])/sqr(sqrt_area);
        break;
        case 5:
            val[1] = 2.0*(crds[1]-cent[1])/sqr(sqrt_area);
        break;
        case 6:
            val[0] = 3.0*sqr(crds[0]-cent[0])/cub(sqrt_area);
        break;
        case 7:
            val[0] = 2.0*(crds[0]-cent[0])*(crds[1]-cent[1])/cub(sqrt_area);
            val[1] = sqr(crds[0]-cent[0])/cub(sqrt_area);
        break;
        case 8:
            val[0] = sqr(crds[1]-cent[1])/cub(sqrt_area);
            val[1] = 2.0*(crds[0]-cent[0])*(crds[1]-cent[1])/cub(sqrt_area);
        break;
        case 9:
            val[1] = 3.0*sqr(crds[1]-cent[1])/cub(sqrt_area);
        break;
        default:
            printf("ERROR grad_vh_loc_basis, implement 2D degree %d\n", indx);
            //clean_up(ERROR); //20191006, not sure how to really use this function so it is commented out for now. hopefully don't really need it.
        }
}

double Utilities::inter_u_integr_7_quad_LDG(
        double    *vt[], // vertices of the tri after affine mapping. 
        double    *cent,  // after affine mapping
        double    area, 
        double    *soln_u, // soln polynomial: u of LDG
        double    crds[][3],
	    double    fluxu[][4])
{
        int       i, j;
        double    *pcrds[3], one_3rd, sqrt_area;
        double    emid0[3], emid1[3], emid2[3];
        // double    *pemid[3];
        double    g_vh[3];
        // double    soln_at_pt[4];

        // one_3rd = 1.0/3.0; 

        // for(j = 0; j < 3; j++)
        //     cent[j] = one_3rd*(vt[0][j] + vt[1][j] + vt[2][j]);

        for(i = 0; i < 3; i++)
        {
            // crds[i][0] = vt[i][0];
            // crds[i][1] = vt[i][1];
            memcpy(crds[i], vt[i], sizeof(double)*3);
        }

        for(i = 0; i < 3; i++)
        {
            emid0[i] = 0.5*(vt[0][i] + vt[1][i]);
            emid1[i] = 0.5*(vt[1][i] + vt[2][i]);
            emid2[i] = 0.5*(vt[2][i] + vt[0][i]);
        }

        // pemid[0] = emid0;
        // pemid[1] = emid1;
        // pemid[2] = emid2;
        memcpy(crds[3], emid0, sizeof(double)*3);
        memcpy(crds[4], emid1, sizeof(double)*3);
        memcpy(crds[5], emid2, sizeof(double)*3);
        memcpy(crds[6], cent, sizeof(double)*3);

        // area = triangle_area_3d(vt[0], vt[1], vt[2]);
        sqrt_area = sqrt(area);

        // Compute u at the point
        for(i = 0; i < 7; i++)
            soln_at_pt(soln_u, crds[i], cent, sqrt_area, fluxu[i]);
}

void Utilities::inter_integr_7_quad_ver2(
	double        *vt[], // vertices of the tri after affine mapping,
    double        *cent,
    double        area,
    int           indx,
	double        crds[][3],
	double        fluxx[][4], //[#quadrature][# eq]
	double        fluxy[][4],
	double        *ans)
{
        int       i, j;
        double    *pcrds[3], one_3rd, sqrt_area;
        double    g_vh[3];
        double    tmpans[8][13];

        for(j = 0; j < N_EQN; j++)
            ans[j] = 0.0;

        if(indx == 0) return;

        // one_3rd = 1.0/3.0;
        // for(j = 0; j < 3; j++)
        //     cent[j] = one_3rd*(vt[0][j] + vt[1][j] + vt[2][j]);

        // area = triangle_area_3d(vt[0], vt[1], vt[2]);
        sqrt_area = sqrt(area);

        for(i = 0; i < 7; i++)
        {
            grad_vh_loc_basis(crds[i], cent, sqrt_area, indx, g_vh);
            for(j = 0; j < N_EQN; j++)
                tmpans[j][i] = (fluxx[i][j]*g_vh[0] + fluxy[i][j]*g_vh[1]);
        }

        for(j = 0; j < N_EQN; j++)
        {
            ans[j] = 0.05*(tmpans[j][0] + tmpans[j][1] + tmpans[j][2]) +
                     2.0/15.0*(tmpans[j][3] + tmpans[j][4] + tmpans[j][5]) +
                     9.0/20.0*tmpans[j][6];
        }

        for(j = 0; j < N_EQN; j++)
            ans[j] *= area;

}

void Utilities::inter_integr_7_quad_ver3(
    double        *prev_soln,
    double        *cent,
    double        area,
    int           indx,
	double        crds[][3],
	double        *ans)
{
        int       i, j;
        double    *pcrds[3], one_3rd, sqrt_area;
        double    g_vh[3];
        double    tmpans[8][13];
        double vh_pt, soln_u[4];

        for(j = 0; j < N_EQN; j++){
            ans[j] = 0.0;
        }

        //if(indx == 0) return;
        sqrt_area = sqrt(area);

        for(i = 0; i < 7; i++)
        {
            // grad_vh_loc_basis(crds[i], cent, sqrt_area, indx, g_vh);
            vh_pt = vh_val_loc_div_free_basis(crds[i], cent, sqrt_area,indx);
            // if (i == 1){
                // std::cout<<"crds["<<i<<"][0] = "<<crds[i][0]<<std::endl;
                // std::cout<<"crds["<<i<<"][1] = "<<crds[i][1]<<std::endl;
                // std::cout<<"crds["<<i<<"][2] = "<<crds[i][2]<<std::endl;
                // std::cout<<"cent = "<<cent[0]<<" "<<cent[1]<<" "<<cent[2]<<std::endl;
                // std::cout<<"sqrt_area = "<<sqrt_area<<std::endl;
            // }
        
            for(j = 0; j < N_EQN; j++){
                // tmpans[j][i] = ((fluxx[i][j]*g_vh[0] + fluxy[i][j]*g_vh[1]));
                soln_at_pt(prev_soln, crds[i], cent, sqrt_area, soln_u);
                tmpans[j][i] = soln_u[0]*vh_pt;
                // std::cout<<"indx = "<<indx<<std::endl;
                // std::cout<<"soln_u[0] = "<<soln_u[0]<<std::endl;
                // std::cout<<"vh_pt = "<<vh_pt<<std::endl;
                // std::cout<<"tmpans["<<j<<"]["<<i<<"]"<<tmpans[j][i]<<std::endl;
            }
        }

        for(j = 0; j < N_EQN; j++)
        {
            ans[j] = 0.05*(tmpans[j][0] + tmpans[j][1] + tmpans[j][2]) +
                     2.0/15.0*(tmpans[j][3] + tmpans[j][4] + tmpans[j][5]) +
                     9.0/20.0*tmpans[j][6];
        }

        for(j = 0; j < N_EQN; j++){
            ans[j] *= area;
        }

}

void Utilities::inter_integr_7_quad_ver4(
    double        source_coef,
    double        *cent,
    double        area,
    int           indx,
	double        crds[][3],
	double        *ans)
{
        int       i, j;
        double    *pcrds[3], one_3rd, sqrt_area;
        double    g_vh[3];
        double    tmpans[8][13];
        double vh_pt, soln_u[4];

        for(j = 0; j < N_EQN; j++){
            ans[j] = 0.0;
        }

        //if(indx == 0) return;
        sqrt_area = sqrt(area);

        for(i = 0; i < 7; i++)
        {
            // grad_vh_loc_basis(crds[i], cent, sqrt_area, indx, g_vh);
            vh_pt = vh_val_loc_div_free_basis(crds[i], cent, sqrt_area,indx);
            // if (i == 1){
                // std::cout<<"crds["<<i<<"][0] = "<<crds[i][0]<<std::endl;
                // std::cout<<"crds["<<i<<"][1] = "<<crds[i][1]<<std::endl;
                // std::cout<<"crds["<<i<<"][2] = "<<crds[i][2]<<std::endl;
                // std::cout<<"cent = "<<cent[0]<<" "<<cent[1]<<" "<<cent[2]<<std::endl;
                // std::cout<<"sqrt_area = "<<sqrt_area<<std::endl;
            // }
        
            for(j = 0; j < N_EQN; j++){
                // tmpans[j][i] = ((fluxx[i][j]*g_vh[0] + fluxy[i][j]*g_vh[1]));
                //soln_at_pt(prev_soln, crds[i], cent, sqrt_area, soln_u);
                tmpans[j][i] = source_coef*vh_pt;
                // std::cout<<"indx = "<<indx<<std::endl;
                // std::cout<<"soln_u[0] = "<<soln_u[0]<<std::endl;
                // std::cout<<"vh_pt = "<<vh_pt<<std::endl;
                // std::cout<<"tmpans["<<j<<"]["<<i<<"]"<<tmpans[j][i]<<std::endl;
            }
        }

        for(j = 0; j < N_EQN; j++)
        {
            ans[j] = 0.05*(tmpans[j][0] + tmpans[j][1] + tmpans[j][2]) +
                     2.0/15.0*(tmpans[j][3] + tmpans[j][4] + tmpans[j][5]) +
                     9.0/20.0*tmpans[j][6];
        }

        for(j = 0; j < N_EQN; j++){
            ans[j] *= area;
        }

}

//This function calculates the integral of k0/(1+(beta*u)^-q)
void Utilities::inter_integr_7_quad_ver5(
    double        k_0,
    double        beta,
    double        q1,
    double        u,
    double        *cent,
    double        area,
    int           indx,
	double        crds[][3],
	double        *ans)
{
        int       i, j;
        double    *pcrds[3], one_3rd, sqrt_area;
        double    g_vh[3];
        double    tmpans[8][13];
        double    vh_pt, soln_u[4];
        double    coef;
        double    q = -q1;

        for(j = 0; j < N_EQN; j++){
            ans[j] = 0.0;
        }

        //if(indx == 0) return;
        sqrt_area = sqrt(area);

        for(i = 0; i < 7; i++)
        {
            vh_pt = vh_val_loc_div_free_basis(crds[i], cent, sqrt_area,indx);
            coef = k_0/(1.0+pow(beta*u, q));
            for(j = 0; j < N_EQN; j++){
                tmpans[j][i] = coef*vh_pt;
            }
        }

        for(j = 0; j < N_EQN; j++)
        {
            ans[j] = 0.05*(tmpans[j][0] + tmpans[j][1] + tmpans[j][2]) +
                     2.0/15.0*(tmpans[j][3] + tmpans[j][4] + tmpans[j][5]) +
                     9.0/20.0*tmpans[j][6];
        }

        for(j = 0; j < N_EQN; j++){
            ans[j] *= area;
        }

}

//This function calculates the integral of k1/(1+(\gamma*p*a)^-h)
void Utilities::inter_integr_7_quad_ver6(
    double        k_1,
    double        beta,
    double        gamma,
    double        q1,
    double        h,
    double        u,
    double        *prev_soln,
    double        *cent,
    double        area,
    int           indx,
	double        crds[][3],
	double        *ans)
{
        int       i, j;
        double    *pcrds[3], one_3rd, sqrt_area;
        double    g_vh[3];
        double    tmpans[8][13];
        double vh_pt, soln_u[4], coef,p;
        double    q = -q1;
        double    h1 = -h;

        for(j = 0; j < N_EQN; j++){
            ans[j] = 0.0;
        }

        //if(indx == 0) return;
        sqrt_area = sqrt(area);

        for(i = 0; i < 7; i++)
        {
            // grad_vh_loc_basis(crds[i], cent, sqrt_area, indx, g_vh);
            vh_pt = vh_val_loc_div_free_basis(crds[i], cent, sqrt_area,indx);
          
        
            for(j = 0; j < N_EQN; j++){
                soln_at_pt(prev_soln, crds[i], cent, sqrt_area, soln_u);
                p = 1.0/(1.0+pow(beta*u,q));
                coef = k_1/(1.0+pow(gamma*p*soln_u[0], h1));
                tmpans[j][i] = coef*vh_pt;
                
            }
        }

        for(j = 0; j < N_EQN; j++)
        {
            ans[j] = 0.05*(tmpans[j][0] + tmpans[j][1] + tmpans[j][2]) +
                     2.0/15.0*(tmpans[j][3] + tmpans[j][4] + tmpans[j][5]) +
                     9.0/20.0*tmpans[j][6];
        }

        for(j = 0; j < N_EQN; j++){
            ans[j] *= area;
        }

}



/*For computing \int_K u_h*div(vec{r}) dx */
/*  inter_u_integr_7_quad_ver2() should be used after
 *  inter_u_integr_7_quad_LDG(). 
 *  crds[][3] constains crds of quadrature pts;
 *  fluxu[][4] has u_h value at quadrature pts.
 */
void Utilities::inter_u_integr_7_quad_ver2(
	double        *vt[], // vertices of the tri after affine mapping,
        double        *cent,
        double        area, 
        int           indx,
	double        crds[][3],
	double        u_flux[][4],
	double        *ans)
{
        int       i, j;
        double    *pcrds[3], one_3rd, sqrt_area;
        double    d_vh[3];
        double    tmpans[4][13];

        for(j = 0; j < N_EQN; j++) ans[j] = 0.0;

        if(indx == 0) return;

        // one_3rd = 1.0/3.0;
        // for(j = 0; j < 3; j++)
        //     cent[j] = one_3rd*(vt[0][j] + vt[1][j] + vt[2][j]);

        // area = triangle_area_3d(vt[0], vt[1], vt[2]);
        sqrt_area = sqrt(area);

        for(i = 0; i < 7; i++)
        {
            // grad_vh_loc_basis(crds[i], cent, sqrt_area, indx, g_vh);
            div_vec_vh_loc_basis(crds[i], cent, sqrt_area, indx, d_vh);
            for(j = 0; j < N_EQN; j++)
                tmpans[j][i] = (u_flux[i][j]*d_vh[0]);
        }

        for(j = 0; j < N_EQN; j++)
        {
            ans[j] = 0.05*(tmpans[j][0] + tmpans[j][1] + tmpans[j][2]) +
                     2.0/15.0*(tmpans[j][3] + tmpans[j][4] + tmpans[j][5]) +
                     9.0/20.0*tmpans[j][6];
        }

        for(j = 0; j < N_EQN; j++)
            ans[j] *= area;
}

void Utilities::div_vec_vh_loc_basis(
        double  *crds,
        double  *cent,
        double  sqrt_area,
        int     indx,
        double  *val)
{
        *val = 0.0;
        switch(indx)
        {
        case 0:
        break;
        case 1:
        break;
        case 2:
            *val = 1.0/sqrt_area;
        break;
        case 3:
        break;
        case 4:
        break;
        case 5:
            *val = 1.0/sqrt_area;
        break;
        case 6:
            *val = 2.0*(crds[0] - cent[0])/sqr(sqrt_area);
        break;
        case 7:
        break;
        case 8:
            *val = (crds[1]-cent[1])/sqr(sqrt_area);
        break;
        case 9:
            *val = (crds[0]-cent[0])/sqr(sqrt_area);
        break;
        case 10:
        break;
        case 11:
            *val = 2.0*(crds[1] - cent[1])/sqr(sqrt_area);
        break;
        default:
            printf("ERROR div_vec_vh_loc_basis(), implement 2D degree %d\n", indx);
            //clean_up(ERROR); //20191006
        }
}

 void Utilities::adv_fw_LDG_update_u(
	std::list<ELEMENT>::iterator &it_ele,
	std::list<ELEMENT>::iterator &n_it_ele,
        std::list<ELEMENT>::iterator *it_adj_eles, 
        std::list<ELEMENT>::iterator *n_it_adj_eles, 
        RBC                     *rbc, 
	double                  dt,
        double                  time, 
        double                  sqrt_diff_coef,
        double                  k_0,
        double                  k_1,
        double                  k_2, //sink_coef1,
        double                  k_3, //sink_coef2,
        double                  k_4,
        double                  k_ss,
        double                  beta,
        double                  gamma,
        double                  q1,
        double                  h,
        // double                  source_coef,
        int                     rk_iter,
        int                     INCLUDE_ADV) /** NEED TO SET INCLUDE_ADV FLAG for different problems **/
{
        double        prev_soln[MAX_N_COEF], prev_nb_solns[3][MAX_N_COEF]; 
                      // from previous rk stepping
        CVector       *vt[3];
        int           i, j, k, indx, side, positive_id;
        double        B[4], affine_pt[3][4], *affine_pt_p[3], aff_cent[4], adj_aff_cent[3][4];
        double        *old_flux_q[2], crds[16][3], 
                      fluxx[16][4], fluxy[16][4], fluxz[16][4];  
        double        *old_nb_flux_q[3][2]; 
        double        q_fluxx[16][4], q_fluxy[16][4], q_fluxz[16][4];  // [][#eq]
        double        i_integr[4], rhs[4][MAX_N_COEF], mulrhs[4][MAX_N_COEF], q_i_integr[4]; 
        double        src_val[16][4], s_integr[4]; 
        double        peint[3][4], q_peint[3][4], u_peint[3][4], area, debug_q_flux_sum[4][MAX_N_COEF];
        // ELEMENT       tmp_ele;
        static double **inv_adj_affine_tran[3] = {NULL,NULL,NULL}, 
                      **aff_inv_adj_aff[3], **inv_affine_tran = NULL, **adj_aff_inv_aff[3];
        // static double **jacobi_x_over_x_hat[3]; 
        double        **adj_ele_affine_trans[3]; 
        int           nbside, side_adj_eles[3];
        double        **Bmass_matrix, **Bmass_inv;
        int           save_at_rk; 
        double        *soln_dg_rho, *old_dg_rho, *prev_rk_dg_rho, *soln_b, *old_b;
        double        u_hat[3][5][10]; //[side][# eq.][# qdr. pt]
        int           debug = NO; 

        bool          COMPUT_aff_cnt_vtx_ON_FLY = false; 
        double        sink_integr[4], source_integr[4], avg_soln[4], k_0_integr[4], k_1_integr[4], unit_soln[4];
        double        u_for_dg, sink_coef1, sink_coef2;

        //DEBUG_ENTER(adv_fw_LDG_update_u)

        /**
        // if(it_ele->Id == 0 || it_ele->Id == 2)
        if(it_ele->Id == 765)
        {
            printf("\n\n--------------------------------------------\n"); 
            printf("\n----------ele[%d] entered adv_fw_LDG_update_u() rk_iter = %d\n", 
                       it_ele->Id, rk_iter);
            debug = YES; 
        }
        **/

       q_i_integr[0] = q_i_integr[1] = q_i_integr[2] = q_i_integr[3] = 0.0;
       q_peint[0][0] = q_peint[0][1] = q_peint[0][2] = q_peint[0][3] = 0.0;
       q_peint[1][0] = q_peint[1][1] = q_peint[1][2] = q_peint[1][3] = 0.0;
       q_peint[2][0] = q_peint[2][1] = q_peint[2][2] = q_peint[2][3] = 0.0;
       sink_integr[0] = sink_integr[1] = sink_integr[2] = sink_integr[3] = 0.0;
       avg_soln[0] = avg_soln[1] = avg_soln[2] = avg_soln[3] = 0.0;
       unit_soln[0] = unit_soln[1] = unit_soln[2] = unit_soln[3] = 0.0;
       k_0_integr[0] = k_0_integr[1] = k_0_integr[2] = k_0_integr[3] = 0.0;
       k_1_integr[0] = k_1_integr[1] = k_1_integr[2] = k_1_integr[3] = 0.0;
        if(NULL == inv_affine_tran)
        {
            FT_matrix(&inv_affine_tran, 4, 4,sizeof(double));
            for(i = 0; i < 3; i++)
                FT_matrix(&inv_adj_affine_tran[i], 4, 4,sizeof(double));
            for(i = 0; i < 3; i++)
                FT_matrix(&aff_inv_adj_aff[i], 4, 4,sizeof(double));
            for(i = 0; i < 3; i++)
                FT_matrix(&adj_aff_inv_aff[i], 4, 4,sizeof(double));
            // FT_matrix(&jacobi_x_over_x_hat, 3, 3,sizeof(double));
        }

        Bmass_matrix = it_ele->mass_matrix;
        Bmass_inv = it_ele->mass_inv;
        area = it_ele->area0;

        save_at_rk = (rk_iter+1)%RK_STEP;
        soln_dg_rho = n_it_ele->dg_rho[save_at_rk]; 
        old_dg_rho = it_ele->dg_rho[0]; 
        soln_b = n_it_ele->b[save_at_rk]; 
        old_b = it_ele->b[0];
        u_for_dg = it_ele->u_for_dg[0];

        if(0 == rk_iter)
        {
            old_flux_q[0] = it_ele->dg_q[0][0]; 
            old_flux_q[1] = it_ele->dg_q[0][1]; 
            // std::cout<<"null_aff_trans_old_flux_q_it_ele = "<<it_ele->dg_q[0][0][0]<<" "<<it_ele->dg_q[0][0][1]<<" "<<it_ele->dg_q[0][0][2]<<std::endl;
            // std::cout<<"null_aff_trans_old_flux_q_it_ele = "<<it_ele->dg_q[0][1][0]<<" "<<it_ele->dg_q[0][1][1]<<" "<<it_ele->dg_q[0][1][2]<<std::endl;
        }
        else
        {
            old_flux_q[0] = n_it_ele->dg_q[rk_iter][0]; 
            old_flux_q[1] = n_it_ele->dg_q[rk_iter][1]; 
            // std::cout<<"null_aff_trans_old_flux_q_n_it_ele = "<<n_it_ele->dg_q[rk_iter][0][0]<<" "<<n_it_ele->dg_q[rk_iter][0][1]<<" "<<n_it_ele->dg_q[rk_iter][0][2]<<std::endl;
            // std::cout<<"null_aff_trans_old_flux_q_n_it_ele = "<<n_it_ele->dg_q[rk_iter][1][0]<<" "<<n_it_ele->dg_q[rk_iter][1][1]<<" "<<n_it_ele->dg_q[rk_iter][1][2]<<std::endl;
        }

        // inverse_matrix(it_ele->affine_trans, 4, inv_affine_tran);
        inv_affine_tran = it_ele->inv_affine_trans;

        // find nb-side numbers
        for(j=0;j<3;j++)
        {
            adj_ele_affine_trans[j] = it_adj_eles[j]->affine_trans; 
            // inverse_matrix(it_adj_eles[j]->affine_trans, 4, inv_adj_affine_tran[j]);
            inv_adj_affine_tran[j] = it_adj_eles[j]->inv_affine_trans;

            // printf("after inverse it_adj_eles[%d] affine_trans = %p\n", 
            //           j, it_adj_eles[j]->affine_trans); fflush(stdout); 
            matrix_matrix_mult(it_adj_eles[j]->affine_trans, inv_affine_tran, 
                               4, 4, adj_aff_inv_aff[j]);
            matrix_matrix_mult(it_ele->affine_trans, inv_adj_affine_tran[j], 
                               4, 4, aff_inv_adj_aff[j]);
            // printf("after it_adj_eles[%d] affine_trans = %p\n", 
            //      j, it_adj_eles[j]->affine_trans);

            for(nbside = 0; nbside < 3; nbside++)
            {
                if(it_adj_eles[j]->adj_ele[nbside] == it_ele->Id)
                {
                    side_adj_eles[j] = nbside;
                    break;
                }
            }    
            if(nbside >= 3)
            {
                printf("ERROR: adv_fw_LDG_update_u(), side[%d]."
                 " adj_ele, ele do not share common edge\n", j);
               // clean_up(ERROR); //20191006
            }
        }

        /*Test code*/
        /*
        if(it_ele->Id == 0)
        {
            double coef[4] = {0.0, 0.0, 0.0, 0.0}; 
            double B[4] = {0.0, 0.0, 1.0, 1.0}, Bprim[4], ans_tmp; 
            coef[0] = 1.0 - (inv_adj_affine_tran[0][0][3] + 
                             inv_adj_affine_tran[0][1][3] + inv_adj_affine_tran[0][2][3]);
            for(j = 0; j < 3; j++)
                coef[1] += inv_adj_affine_tran[0][j][0]; 
            coef[1] *= -1.0;
            for(j = 0; j < 3; j++)
                coef[2] += inv_adj_affine_tran[0][j][1]; 
            coef[2] *= -1.0;
            for(j = 0; j < 3; j++)
                coef[3] += inv_adj_affine_tran[0][j][2]; 
            coef[3] *= -1.0;
            matrix_vec_mult(adj_ele_affine_trans[0], B, 4, 4, Bprim);
            ans_tmp = coef[0] + coef[1]*Bprim[0] + coef[2]*Bprim[1] + coef[3]*Bprim[2];
            printf("----- Result = %g\n", ans_tmp); 
            matrix_vec_mult(inv_affine_tran[0], Bprim, 4, 4, B);
            print_vector("B after tran back", 4, B, "%g "); 
        }
        printf("WARNING adv_fw. stop after test code\n");
        // clean_up(0); 
        **/
        /*END: Test code*/

        if(0 == rk_iter)
        {
            memcpy(prev_soln, it_ele->dg_rho[rk_iter], sizeof(double)*MAX_N_COEF); 
            for(j=0;j<3;j++)
            {
                memcpy(prev_nb_solns[j], it_adj_eles[j]->dg_rho[rk_iter], 
                       sizeof(double)*MAX_N_COEF); 

            }

            for(j=0;j<3;j++)
            {
                old_nb_flux_q[j][0] = it_adj_eles[j]->dg_q[rk_iter][0];  
                old_nb_flux_q[j][1] = it_adj_eles[j]->dg_q[rk_iter][1];
                // std::cout<<"old_nb_flux_q_rk_iter0 = "<<it_adj_eles[j]->dg_q[rk_iter][0][0]<<" "<<it_adj_eles[j]->dg_q[rk_iter][0][1]<<" "<<it_adj_eles[j]->dg_q[rk_iter][0][2]<<std::endl;
                // std::cout<<"old_nb_flux_q_rk_iter0 = "<<it_adj_eles[j]->dg_q[rk_iter][1][0]<<" "<<it_adj_eles[j]->dg_q[rk_iter][1][1]<<" "<<it_adj_eles[j]->dg_q[rk_iter][1][2]<<std::endl;
                
            }
        }
        else
        {
            memcpy(prev_soln, n_it_ele->dg_rho[rk_iter], sizeof(double)*MAX_N_COEF); 
            for(j=0;j<3;j++)
            {
                memcpy(prev_nb_solns[j], n_it_adj_eles[j]->dg_rho[rk_iter], 
                       sizeof(double)*MAX_N_COEF); 
            }

            for(j=0;j<3;j++)
            {
                old_nb_flux_q[j][0] = n_it_adj_eles[j]->dg_q[rk_iter][0];  
                old_nb_flux_q[j][1] = n_it_adj_eles[j]->dg_q[rk_iter][1];
                // std::cout<<"old_nb_flux_q_rk_iter1 = "<<n_it_adj_eles[j]->dg_q[rk_iter][0][0]<<" "<<n_it_adj_eles[j]->dg_q[rk_iter][0][1]<<" "<<n_it_adj_eles[j]->dg_q[rk_iter][0][2]<<std::endl;
                // std::cout<<"old_nb_flux_q_rk_iter1 = "<<n_it_adj_eles[j]->dg_q[rk_iter][1][0]<<" "<<n_it_adj_eles[j]->dg_q[rk_iter][1][1]<<" "<<n_it_adj_eles[j]->dg_q[rk_iter][1][2]<<std::endl;  
            }
        }

        

        /*Transform it_ele vertex coords from global coords to local coords.*/
        if(true == COMPUT_aff_cnt_vtx_ON_FLY){
          
            for(j=0;j<3;j++)
            {
                memcpy(B, vt[j]->x, sizeof(double)*3);
                B[3] = 1.0;
                matrix_vec_mult(it_ele->affine_trans, B, 4, 4, affine_pt[j]);
                affine_pt_p[j] = affine_pt[j]; 
            }

            memcpy(B, it_ele->center, sizeof(double)*3);
            B[3] = 1.0;
            matrix_vec_mult(it_ele->affine_trans, B, 4, 4, aff_cent);
        }
        else{
            memcpy(aff_cent, it_ele->aff_cent, sizeof(double)*4);
            for(j=0;j<3;j++)
                memcpy(affine_pt[j], it_ele->affine_pt[j],sizeof(double)*4);
            for(j=0;j<3;j++) affine_pt_p[j] = affine_pt[j];
        }

        /*Transform adjacent ele center to its own local coords. */
        // B[3] = 1.0;
        if(true == COMPUT_aff_cnt_vtx_ON_FLY){
            for(side = 0; side < 3; side++) {
                memcpy(B, it_adj_eles[side]->center, sizeof(double)*3);  
                B[3] = 1.0;
                matrix_vec_mult(it_adj_eles[side]->affine_trans, B, 4, 4, adj_aff_cent[side]);
            }
        }
        else{
            for(side = 0; side < 3; side++)
                memcpy(adj_aff_cent[side], it_adj_eles[side]->aff_cent, sizeof(double)*4);
        }

        /*Compute \int_K vec{q} \cdot grad(v_h) dx */
        inter_integr_7_quad_flux_LDG(affine_pt_p, aff_cent, area, 
                 it_ele->affine_trans, inv_affine_tran, 
                 old_flux_q, crds, q_fluxx, q_fluxy, q_fluxz); 
        // std::cout<<"crds[1][0] = "<<crds[1][0]<<std::endl;
        // std::cout<<"crds[1][1] = "<<crds[1][1]<<std::endl;
        // std::cout<<"crds[1][2] = "<<crds[1][2]<<std::endl;
        // std::cout<<"affine_pt_p[0][0] = "<<affine_pt_p[0][0]<<std::endl;
        // std::cout<<"affine_pt_p[0][1] = "<<affine_pt_p[0][1]<<std::endl;
        // std::cout<<"affine_pt_p[0][2] = "<<affine_pt_p[0][2]<<std::endl;

        // if(YES == INCLUDE_ADV)
        // {
        //     /*inter_integr_7_quad_flux(affine_pt_p, aff_cent, area, it_ele->affine_trans, 
        //       inv_affine_tran, prev_soln, crds, fluxx, fluxy, fluxz, debug); 

        //     inter_integr_7_quad_source(affine_pt_p, aff_cent, area, rbc->Center.x, time, 
        //       it_ele->affine_trans, inv_affine_tran, prev_soln, crds, src_val, debug); 
        //     */
        // }

        for(indx = 0; indx < MAX_N_COEF; indx++)
        {
            for(i = 0; i < N_EQN; i++) i_integr[i] = 0.0;
            for(i = 0; i < N_EQN; i++) s_integr[i] = 0.0;
            for(i = 0; i < N_EQN; i++) debug_q_flux_sum[i][indx] = 0.0; 

            // if(YES == INCLUDE_ADV)
            // {
            //     inter_integr_7_quad_ver2(affine_pt_p,aff_cent, area,indx,
            //                              crds,fluxx,fluxy,i_integr);
            //     inter_integr_7_quad__source_ver2(affine_pt_p,aff_cent, area,indx,
            //                              crds, src_val,s_integr);
            // }

            inter_integr_7_quad_ver2(affine_pt_p, aff_cent, area, indx,
                                     crds,q_fluxx,q_fluxy,q_i_integr);

            if ((rk_iter)==0){
                inter_integr_7_quad_ver3(old_dg_rho, aff_cent, 
                                    area, indx, crds,sink_integr);
            }
            else{
                inter_integr_7_quad_ver3(//soln_dg_rho, 
                                    n_it_ele->dg_rho[rk_iter], aff_cent, 
                                    area, indx, crds,sink_integr);
            }

            //inter_integr_7_quad_ver4(source_coef, aff_cent, area, indx, crds, source_integr);

            //// Introduce the k0/(1+(beta*u)^-q) term
            inter_integr_7_quad_ver5(k_0, beta, q1, u_for_dg, aff_cent,area,indx,crds,k_0_integr);
            //// Introduce the k1/(1+(gamma*p*a)^-q) term
            if ((rk_iter)==0){
                inter_integr_7_quad_ver6(k_1,beta,gamma,q1,h,u_for_dg,old_dg_rho,aff_cent,area,indx,crds,k_1_integr);
            }
            else{
                inter_integr_7_quad_ver6(k_1,beta,gamma,q1,h,u_for_dg,n_it_ele->dg_rho[rk_iter],aff_cent,area,indx,crds,k_1_integr);
            }
            //// Introduce k2*a term and k3*b*a term as a single sink term
            //sink_coef = k_2 - k_3*b    
            sink_coef1 = k_2;
            sink_coef2 = k_3*old_b[0];
            if (it_ele->Id == 0 || it_ele->Id == 30){
                if (indx==0){
                    // std::cout<<"old_dg_rho of triangle["<<it_ele->Id<<"] = "<<old_dg_rho[0]<<std::endl;
                    // std::cout<<"soln_dg_rho of triangle["<<it_ele->Id<<"] = "<<soln_dg_rho[0]<<std::endl;
                    // std::cout<<"k_0_integr of triangle["<<it_ele->Id<<"] = "<<k_0_integr[0]<<std::endl;
                    // std::cout<<"k_1_integr of triangle["<<it_ele->Id<<"] = "<<k_1_integr[0]<<std::endl;
                }
            }

            for(i = 0; i < N_EQN; i++) {rhs[i][indx] = i_integr[i] + s_integr[i];}

            for(i = 0; i < N_EQN; i++){
                //  std::cout<<"rhs[i][indx] before adding q_i_integr = "<<rhs[i][indx]<<std::endl;
                //  std::cout<<"q_i_integr[0] = "<<q_i_integr[i]<<std::endl;
                rhs[i][indx] -= q_i_integr[i]*sqrt_diff_coef;
                rhs[i][indx] -= sink_integr[i]*sink_coef1;
               rhs[i][indx] -= sink_integr[i]*sink_coef2;
                // rhs[i][indx] += source_integr[i];
            //    rhs[i][indx] += (k_0/(1.0+pow(beta*u_for_dg,-q1)))*it_ele->area;
                  rhs[i][indx] += k_0_integr[i];
               rhs[i][indx] += k_1_integr[i];
            }

            /* Compute diffusion flux on edge */
            for(side = 0; side < 3; side++)
            {
                edge_integr_q_hat_flux(it_ele, it_adj_eles[side], prev_soln, 
                        prev_nb_solns[side], old_flux_q, old_nb_flux_q[side], 
                        side, indx, 
                        q_peint[side], dt, rk_iter, affine_pt_p, aff_cent, 
                        adj_aff_inv_aff[side], aff_inv_adj_aff[side], 
                        adj_aff_cent[side], inv_affine_tran, inv_adj_affine_tran[side], 
                        rbc); 
            }//END::: for(side = 0; side < 3; side++)

            for(i = 0; i < N_EQN; i++)
            {
                for(side = 0; side < 3; side++){
                    rhs[i][indx] += sqrt_diff_coef*q_peint[side][i];
                }

                for(side = 0; side < 3; side++)
                    debug_q_flux_sum[i][indx] += q_peint[side][i];
            }

            // if(YES == debug && indx == 0)
            

        }//END::: for(indx = 0; indx < MAX_N_COEF; indx++)

        /* RK TIME STEPPING */
        if(1 == RK_STEP || 2 == RK_STEP || 3 == RK_STEP)
        {
            for(i = 0; i < N_EQN; i++){
                matrix_vec_mult(Bmass_inv, rhs[i], MAX_N_COEF, MAX_N_COEF, mulrhs[i]);
            }
        }

        /*Here is where computed solution of u is updated into the dg_rho[1][indx] and dg_rho[0][indx] based on the current
        rk_iter*/
        if(1 == RK_STEP)
        {
            for(indx = 0; indx < MAX_N_COEF; indx++)
                soln_dg_rho[indx] = old_dg_rho[indx] + dt*mulrhs[0][indx];

        }
        else if(2 == RK_STEP)
        {   // old_dg_rho = it_ele->dg_rho[0];
            if(0 == rk_iter)
            {
                for(indx = 0; indx < MAX_N_COEF; indx++){
                    if (it_ele->Id == 5){
                        // std::cout<<"dt = "<<dt<<", mulrhs = "<<mulrhs[0][indx]<<std::endl;
                    }
                    soln_dg_rho[indx] = old_dg_rho[indx] + dt*mulrhs[0][indx];
                    // std::cout<<"rk0_old dg rho "<<old_dg_rho[0]<<std::endl;
                    // std::cout<<"rk0_dt "<<dt<<std::endl;
                    // std::cout<<"rk0_mulrhs "<<mulrhs[0][0]<<std::endl;
                    n_it_ele->dg_rho[(rk_iter+1)% RK_STEP][indx] = soln_dg_rho[indx];//added by KT, 06/28/2020
                }
            }
            else if(1 == rk_iter)
            {
                prev_rk_dg_rho = n_it_ele->dg_rho[rk_iter]; // from prev. rk stepping.
                for(indx = 0; indx < MAX_N_COEF; indx++)
                {   
                    soln_dg_rho[indx] = 0.5*old_dg_rho[indx] + 0.5*prev_rk_dg_rho[indx] + 
                                        0.5*dt*mulrhs[0][indx];
                    // std::cout<<"rk1_old_dg_rho : "<<old_dg_rho[0]<<std::endl;
                    // std::cout<<"rk1_prev_rk_dg_rho : "<<prev_rk_dg_rho[0]<<std::endl;
                    // std::cout<<"rk1_dt*mulrhs : "<<dt*mulrhs[0][0]<<std::endl;
                    it_ele->dg_rho[(rk_iter+1)% RK_STEP][indx] = soln_dg_rho[indx];//added by KT, 06/28/2020
                    n_it_ele->dg_rho[(rk_iter+1)% RK_STEP][indx] = soln_dg_rho[indx];
                    //After soln_dg_rho is update succesfully, we update b based on the updated value of dg_rho
                    // inter_integr_7_quad_ver3(soln_dg_rho, aff_cent, area, indx, crds, avg_soln);
                    // inter_integr_7_quad_ver4(1.0, aff_cent, area, indx, crds, unit_soln);
                    // soln_b[0] = old_b[0] + dt*(k_4*(avg_soln[0]/unit_soln[0] - k_ss)*old_b[0]);
                    // // if (it_ele->Id == 0 && indx == 0){std::cout<<"avg_soln = "<<avg_soln[0]<<std::endl;}
                    // // if (it_ele->Id == 0){std::cout<<"unit_soln = "<<unit_soln[0]<<std::endl;}
                    // // if (it_ele->Id == 0 && indx == 0){std::cout<<"avg_soln/unit_soln = "<<avg_soln[0]/it_ele->area<<std::endl;}
                    // it_ele->b[(rk_iter+1)% RK_STEP][indx] = soln_b[indx];
                    // n_it_ele->b[(rk_iter+1)% RK_STEP][indx] = soln_b[indx];
                }
            }
        }//END::: if(RK_STEP == 2)
        else 
        {
            printf("ERROR: adv_fw_LDG_update_u(), implement RK_STEP = %d\n", RK_STEP);
            //clean_up(ERROR); //20191006
        }

        // printf("WARNING adv_fw_LDG_update_u. stop after test ele %d\n", it_ele->Id);
        // clean_up(0); 

        // DEBUG_LEAVE(adv_fw_LDG_update_u)
}

void Utilities::adv_fw_LDG_update_q(
	std::list<ELEMENT>::iterator &it_ele,
	std::list<ELEMENT>::iterator &n_it_ele, // new element
        std::list<ELEMENT>::iterator *it_adj_eles, 
        std::list<ELEMENT>::iterator *n_it_adj_eles, 
        RBC                     *rbc, 
	double                  dt,
        double                  sqrt_diff_coef,
        int                     rk_iter)
{
        double        cur_soln[MAX_N_COEF], cur_nb_solns[3][MAX_N_COEF]; 
        CVector       *vt[3];
        int           i, j, k, indx, side;
        double        B[4], affine_pt[3][4], *affine_pt_p[3], aff_cent[4], 
                      adj_aff_cent[3][4], area;
        double        crds[16][3], fluxx[16][4], fluxy[16][4], fluxz[16][4];  
        // double        *old_nb_flux_q[3][2], *old_flux_q[2]; 
        // double        q_fluxx[16][4], q_fluxy[16][4], q_fluxz[16][4];  
        double        u_flux[16][4];
        double        i_integr[4], rhs[4][64], mulrhs[4][64], u_i_integr[4]; 
        double        peint[3][4], u_peint[3][4];
        static double **inv_adj_affine_tran[3] = {NULL,NULL,NULL}, 
                      **aff_inv_adj_aff[3], **inv_affine_tran = NULL, **adj_aff_inv_aff[3];
        // static double **jacobi_x_over_x_hat[3];
        double        **adj_ele_affine_trans[3]; 
        int           nbside, side_adj_eles[3];
        double        **vec_vh_mass_matrix, **vec_vh_mass_inv;
        int           save_at_rk; 
        double        *soln_dg_rho, *old_dg_rho, *prev_rk_dg_rho;
        double        u_hat[3][5][10]; //[side][# eq.][# qdr. pt]
        int           debug = NO; 
        double        dg_q_soln[3][3]; 
        double        debug_u_flux_sum[4][64]; 

        bool          COMPUT_aff_cnt_vtx_ON_FLY = false; 

        //DEBUG_ENTER(adv_fw_LDG_update_q)
        // printf("\nele %d newID %d entered adv_fw_LDG_update_q()\n", 
        //           it_ele->Id, n_it_ele->Id); 
        // fflush(stdout);

        /**
        if(it_ele->Id == 0 || it_ele->Id == 2)
        // if(it_ele->Id == 0)
        {
            printf("\n\n&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n");
            printf("\n\n&&&&&&&&  ele[%d] entered adv_fw_LDG_update_q()"
              " rk_iter = %d\n\n", it_ele->Id, rk_iter);
            debug = YES; 
        }
        **/

        if(NULL == inv_affine_tran)
        {
            FT_matrix(&inv_affine_tran, 4, 4,sizeof(double));
            for(i = 0; i < 3; i++)
                FT_matrix(&inv_adj_affine_tran[i], 4, 4,sizeof(double));
            for(i = 0; i < 3; i++)
                FT_matrix(&aff_inv_adj_aff[i], 4, 4,sizeof(double));
            for(i = 0; i < 3; i++)
                FT_matrix(&adj_aff_inv_aff[i], 4, 4,sizeof(double));
            // FT_matrix(&jacobi_x_over_x_hat, 3, 3,sizeof(double));
        }

        vec_vh_mass_matrix = it_ele->vec_vh_mass_matrix;
        vec_vh_mass_inv = it_ele->vec_vh_mass_inv;
        area = it_ele->area0;

        /*
        if(YES == debug)
        {
            double ** tmpI; 
            FT_matrix(&tmpI, MAX_N_COEF*2, MAX_N_COEF*2, sizeof(double)); 
            printf("\n------print vec_vh_mass_matrix && its inv----area = %g\n", area); 
            print_ldb_matrix("vec_vh_mass_matrix", 
                MAX_N_COEF*2, MAX_N_COEF*2,vec_vh_mass_matrix,"%9.8g "); 
            print_ldb_matrix("vec_vh_mass_inv", 
                MAX_N_COEF*2, MAX_N_COEF*2,vec_vh_mass_inv,"%9.8g "); 
            matrix_matrix_mult(vec_vh_mass_matrix, vec_vh_mass_inv, 
                               MAX_N_COEF*2, MAX_N_COEF*2, tmpI); 
            print_ldb_matrix("tmpI", MAX_N_COEF*2, MAX_N_COEF*2, tmpI, "%9.8g ");
            free(tmpI); 
        }
        */

        save_at_rk = (rk_iter+1)%RK_STEP;
        soln_dg_rho = n_it_ele->dg_rho[save_at_rk]; 
        // old_dg_rho = it_ele->dg_rho[rk_iter]; 
        // old_flux_q[0] = it_ele->dg_q[rk_iter][0]; 
        // old_flux_q[1] = it_ele->dg_q[rk_iter][1]; 

        // inverse_matrix(it_ele->affine_trans, 4, inv_affine_tran);
        inv_affine_tran = it_ele->inv_affine_trans;

        // find nb-side numbers
        for(j=0;j<3;j++)
        {
            adj_ele_affine_trans[j] = it_adj_eles[j]->affine_trans; 
            // inverse_matrix(it_adj_eles[j]->affine_trans, 4, inv_adj_affine_tran[j]);
            inv_adj_affine_tran[j] = it_adj_eles[j]->inv_affine_trans;

            // printf("after inverse it_adj_eles[%d] affine_trans = %p\n", 
            //           j, it_adj_eles[j]->affine_trans); fflush(stdout); 
            matrix_matrix_mult(it_adj_eles[j]->affine_trans, inv_affine_tran, 4, 4, 
             adj_aff_inv_aff[j]);
            matrix_matrix_mult(it_ele->affine_trans, inv_adj_affine_tran[j], 4, 4, 
             aff_inv_adj_aff[j]);
            // printf("after it_adj_eles[%d] affine_trans = %p\n", 
            //         j, it_adj_eles[j]->affine_trans);

            for(nbside = 0; nbside < 3; nbside++)
            {
                if(it_adj_eles[j]->adj_ele[nbside] == it_ele->Id)
                {
                    side_adj_eles[j] = nbside;
                    break;
                }
            }    
            if(nbside >= 3)
            {
                printf("ERROR: adv_fw_LDG(), side[%d]."
                       " adj_ele, ele do not share common edge\n", j);
                //clean_up(ERROR); //20191006
            }
        }

        memcpy(cur_soln, n_it_ele->dg_rho[(rk_iter+1)%RK_STEP], 
               sizeof(double)*MAX_N_COEF); 
        for(j=0;j<3;j++)
        {
            memcpy(cur_nb_solns[j], n_it_adj_eles[j]->dg_rho[(rk_iter+1)%RK_STEP], 
                       sizeof(double)*MAX_N_COEF); 
        }

        // for(j=0;j<3;j++) {
            // old_nb_flux_q[j][0] = it_adj_eles[j]->dg_q[rk_iter][0];  
            // old_nb_flux_q[j][1] = it_adj_eles[j]->dg_q[rk_iter][1];  
        // }

        /*Transform it_ele vertex coords from global coords to local coords.*/
        if(true == COMPUT_aff_cnt_vtx_ON_FLY){
            /*vt[0]->x[0] = coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_1[i]];
            vt[0]->x[1] = coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_1[i]];
            vt[0]->x[2] = coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_1[i]];
            vt[1]->x[0] = coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_2[i]];
            vt[1]->x[1] = coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_2[i]];
            vt[1]->x[2] = coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_2[i]];
            vt[2]->x[0] = coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_3[i]];
            vt[2]->x[1] = coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_3[i]];
            vt[2]->x[2] = coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_3[i]];*/
            for(j=0;j<3;j++)
                vt[j] = &(rbc->node[it_ele->vertex[j]]);
            for(j=0;j<3;j++)
            {
                memcpy(B, vt[j]->x, sizeof(double)*3);
                B[3] = 1.0;
                matrix_vec_mult(it_ele->affine_trans, B, 4, 4, affine_pt[j]);
                affine_pt_p[j] = affine_pt[j]; 
            }

            // for(j=0;j<3;j++) B[j] = it_ele->center[j];
            memcpy(B, it_ele->center, sizeof(double)*3);
            B[3] = 1.0;
            matrix_vec_mult(it_ele->affine_trans, B, 4, 4, aff_cent);
        }
        else{
            memcpy(aff_cent, it_ele->aff_cent, sizeof(double)*4); 
            for(j=0;j<3;j++)
                memcpy(affine_pt[j], it_ele->affine_pt[j],sizeof(double)*4);
            for(j=0;j<3;j++) affine_pt_p[j] = affine_pt[j]; 
        }

        /*Transform adjacent ele center to its own local coords. */
        if(true == COMPUT_aff_cnt_vtx_ON_FLY){
            // B[3] = 1.0;
            for(side = 0; side < 3; side++) {
                memcpy(B, it_adj_eles[side]->center, sizeof(double)*3);  
                B[3] = 1.0;
                matrix_vec_mult(it_adj_eles[side]->affine_trans, B, 4, 4, adj_aff_cent[side]);
            }
        }
        else{
            for(side = 0; side < 3; side++)
                memcpy(adj_aff_cent[side], it_adj_eles[side]->aff_cent, sizeof(double)*4);
        }

        /* LDG */
        inter_u_integr_7_quad_LDG(affine_pt_p, aff_cent, area, cur_soln, crds, u_flux); 

        for(indx = 0; indx < 2*MAX_N_COEF; indx++)
        {
            for(i = 0; i < N_EQN; i++){ u_i_integr[i] = 0.0;}
            for(i = 0; i < N_EQN; i++){debug_u_flux_sum[i][indx] = 0.0;}         

            inter_u_integr_7_quad_ver2(affine_pt_p, aff_cent, area, 
                         indx, crds, u_flux, u_i_integr);

            if(YES == debug && indx == 0)
            {
                printf("flux integration on volume = %g\n", u_i_integr[0]); 
            }

            for(i = 0; i < N_EQN; i++)
                rhs[i][indx] = -sqrt_diff_coef*u_i_integr[i];

            /* Compute diffusion flux */
            for(side = 0; side < 3; side++)
            {
                edge_integr_new_u_hat_flux(it_ele, it_adj_eles[side], 
                        cur_soln, cur_nb_solns[side], side, indx, 
                        u_peint[side], dt, rk_iter, affine_pt_p, aff_cent, 
                        adj_aff_inv_aff[side], aff_inv_adj_aff[side], 
                        adj_aff_cent[side], inv_affine_tran, inv_adj_affine_tran[side], rbc,
                        it_ele->uhat_f_flag, it_adj_eles[side]->uhat_f_flag,
                        it_ele->u_hat_flux_store, it_adj_eles[side]->u_hat_flux_store); 
            }//END::: for(side = 0; side < 3; side++)

            for(i = 0; i < N_EQN; i++)
            {
                for(side = 0; side < 3; side++)
                    rhs[i][indx] += sqrt_diff_coef*u_peint[side][i];
                for(side = 0; side < 3; side++)
                    debug_u_flux_sum[i][indx] += u_peint[side][i];
            }

            // if(YES == debug && indx == 0)
            if(YES == debug)
            {
                printf("*** rhs of DG eq -q for indx [%d] = %14.12g u_flux_sum = %14.12g\n\n", 
                     indx, rhs[0][indx], debug_u_flux_sum[0][indx]); 
            }

        }//END::: for(indx = 0; indx < MAX_N_COEF; indx++)

        for(i = 0; i < N_EQN; i++)
            matrix_vec_mult(vec_vh_mass_inv, rhs[i], MAX_N_COEF*2, MAX_N_COEF*2, mulrhs[i]);

        for(i = 0; i < N_EQN; i++)
        {
            for(indx = 0; indx < MAX_N_COEF; indx++)
            {
            // n_it_ele->dg_rho[save_at_rk];
                if ((rk_iter+1)%RK_STEP == 1){
                    n_it_ele->dg_q[save_at_rk][0][indx] = mulrhs[i][indx*2]; 
                    n_it_ele->dg_q[save_at_rk][1][indx] = mulrhs[i][indx*2+1]; 
                }
                if ((rk_iter+1)%RK_STEP == 0){
                    it_ele->dg_q[save_at_rk][0][indx] = mulrhs[i][indx*2]; 
                    it_ele->dg_q[save_at_rk][1][indx] = mulrhs[i][indx*2+1]; 
                }
            }
        }

        // printf("WARNING adv_fw_LDG_update_q. stop after test ele %d\n", it_ele->Id);
        // clean_up(0); 

        DEBUG_LEAVE(adv_fw_LDG_update_q)
}

void Utilities::Comp_ele_affine_trans(RBC* rbc, CoordInfoVecs& coordInfoVecs)
{
    //std::cout<<"ERROR 1"<<std::endl;
	int           i,j,k,m;
	CVector       *vt[3];
	double        dtemp, old_center[3], dtemp_max, old_center_max[3];
        int           dir_max, side_max, k_max; 
        double        *tri_norm, tmp_vec[4], tang_x[4], tang_y[4], len, tan_tmp[4];
        static double **rot_matr_rows = NULL, **translate;
        double        Eu_crds[3][3] = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
        double        B[4], soln[3][4], **inv_affine_tran, **jacobi_x_over_x_hat; 
        double        det2, theta; 

        // inverse_matrix(it_ele->affine_trans, 4, inv_affine_tran);
        // one_3rd = 1.0/3.0;
    //std::cout<<"ERROR 2"<<std::endl;
        FT_matrix(&inv_affine_tran, 4, 4,sizeof(double));
        // FT_matrix(&jacobi_x_over_x_hat, 3, 3,sizeof(double)); /* Probably don't need this for the current work.*/

        if(NULL == rot_matr_rows)
        {
            FT_matrix(&rot_matr_rows, 4, 4,sizeof(double));
            FT_matrix(&translate, 4, 4,sizeof(double));
            for(j = 0; j < 4; j++)
                for(k = 0; k < 4; k++)
                    translate[j][k] = 0.0;
        }
    //std::cout<<"ERROR 3"<<std::endl;
	rbc->Curr_it=rbc->Curr_ele->begin();
	rbc->Curr_ind_it=rbc->Curr_ele_index->begin();
    //std::cout<<"ERROR 4"<<std::endl;
	for(i=0;(rbc->Curr_it!=rbc->Curr_ele->end());i++)
	{
        if (i >= rbc->ENum){
            break;
        }
        if (rbc->Curr_it->vertex[0] >= (INT_MAX-100) || rbc->Curr_it->vertex[1] >= (INT_MAX-100) || rbc->Curr_it->vertex[2] >= (INT_MAX-100)){
            (rbc->Curr_it)++;
            continue;
        }
        if (rbc->Curr_it->vertex[0] <= (-INT_MAX+100) || rbc->Curr_it->vertex[1] <= (-INT_MAX+100) || rbc->Curr_it->vertex[2] <= (-INT_MAX+100)){
            (rbc->Curr_it)++;
            continue;
        }
        //std::cout<<"i = "<<i<<std::endl;     
        
	    for(j=0;j<3;j++){
	        vt[j] = &(rbc->node[rbc->Curr_it->vertex[j]]);
            
        }
	        tri_norm = (rbc->Curr_it->out_norm.x);
         // std::cout<<"tri_norm = "<<rbc->Curr_it->out_norm.x[0]<<std::endl;
         //  std::cout<<"tri_norm = "<<rbc->Curr_it->out_norm.x[1]<<std::endl;
         //  std::cout<<"tri_norm = "<<rbc->Curr_it->out_norm.x[2]<<std::endl;
           
        

            /*
            for(j=0;j<3;j++)
                tmp_vec[j] = (vt[1]->x[j] - vt[0]->x[j]);
            len = Mag3d(tmp_vec);
            for(j=0;j<3;j++) tang_x[j] = tmp_vec[j]/len;
            Cross3d(tri_norm, tang_x, tan_tmp);
            len = Mag3d(tan_tmp);
            for(j=0;j<3;j++) tang_y[j] = tan_tmp[j]/len;

            rot_matr_rows[0][0] = Dot3d(tang_x, Eu_crds[0]);
            rot_matr_rows[0][1] = Dot3d(tang_x, Eu_crds[1]);
            rot_matr_rows[0][2] = Dot3d(tang_x, Eu_crds[2]);
            rot_matr_rows[0][3] = 0.0;

            rot_matr_rows[1][0] = Dot3d(tang_y, Eu_crds[0]);
            rot_matr_rows[1][1] = Dot3d(tang_y, Eu_crds[1]);
            rot_matr_rows[1][2] = Dot3d(tang_y, Eu_crds[2]);
            rot_matr_rows[1][3] = 0.0;

            rot_matr_rows[2][0] = Dot3d(tri_norm, Eu_crds[0]);
            rot_matr_rows[2][1] = Dot3d(tri_norm, Eu_crds[1]);
            rot_matr_rows[2][2] = Dot3d(tri_norm, Eu_crds[2]);
            rot_matr_rows[2][3] = 0.0;

            rot_matr_rows[3][0] = 0.0;
            rot_matr_rows[3][1] = 0.0;
            rot_matr_rows[3][2] = 0.0;
            rot_matr_rows[3][3] = 1.0;

            translate[0][0] = 1.0;
            translate[1][1] = 1.0;
            translate[2][2] = 1.0;
            translate[3][3] = 1.0;
            translate[0][3] = -vt[0]->x[0];
            translate[1][3] = -vt[0]->x[1];
            translate[2][3] = -vt[0]->x[2];
            */

            /*
            if(0 == rbc->Curr_it->Id)
            {
                printf("Comp_ele_affine_trans(). elem[%d] affine_trans = %p\n",
                       rbc->Curr_it->Id, rbc->Curr_it->affine_trans); fflush(stdout); 
            }
            */
     //      std::cout<<"IS ERROR HERE 2?"<<std::endl;
            if(NULL == rbc->Curr_it->affine_trans) {
                FT_matrix(&(rbc->Curr_it->affine_trans), 4, 4,sizeof(double));
                FT_matrix(&(rbc->Curr_it->inv_affine_trans), 4, 4,sizeof(double));
            }
        //    std::cout<<"IS ERROR HERE 3?"<<std::endl;
            // matrix_matrix_mult(rot_matr_rows, translate, 4, 4, rbc->Curr_it->affine_trans);
            Tri_affine_trans(vt[0]->x, vt[1]->x, vt[2]->x, tri_norm, rbc->Curr_it->affine_trans);
       //     std::cout<<"IS ERROR HERE 4?"<<std::endl;
            inverse_matrix(rbc->Curr_it->affine_trans, 4, rbc->Curr_it->inv_affine_trans);
        //   std::cout<<"IS ERROR HERE 5?"<<std::endl;
            /**
            // Test code
            if(i < 1)
            {
                inverse_matrix(rbc->Curr_it->affine_trans, 4, inv_affine_tran);
                for(j=0;j<3;j++)
                {
                    for(k =0; k <3; k++)
                        jacobi_x_over_x_hat[j][k] = inv_affine_tran[j][k];
                }
                det2 = Det3d(jacobi_x_over_x_hat[0],jacobi_x_over_x_hat[1],jacobi_x_over_x_hat[2]);

                for(j=0;j<3;j++) B[j] = vt[0]->x[j];
                B[3] = 1.0;
                matrix_vec_mult(rbc->Curr_it->affine_trans, B, 4, 4, soln[0]);
                for(j=0;j<3;j++) B[j] = vt[1]->x[j];
                matrix_vec_mult(rbc->Curr_it->affine_trans, B, 4, 4, soln[1]);
                for(j=0;j<3;j++) B[j] = vt[2]->x[j];
                matrix_vec_mult(rbc->Curr_it->affine_trans, B, 4, 4, soln[2]);
                len = triangle_area_3d(soln[0], soln[1], soln[2]);

                theta = triangle_area_3d(vt[0]->x, vt[1]->x, vt[2]->x);

                printf("\n---3rd-****----Ele[%d], det_of_Jacobi %12.10g, tri-area %12.10g, %12.10g\n",
                          i, det2, len, theta);
                for(j=0;j<3;j++) B[j] = vt[0]->x[j];
                B[3] = 1.0;
                matrix_vec_mult(rbc->Curr_it->affine_trans, B, 4, 4, soln[0]);
                print_vector("vt0", 4, soln[0], "%g ");
                matrix_vec_mult(inv_affine_tran, soln[0], 4, 4, B);
                vt[0]->Printf("vt0-orig");
                print_vector("vt0---", 4, B, "%13.11g ");

                for(j=0;j<3;j++) B[j] = vt[1]->x[j];
                B[3] = 1.0;
                matrix_vec_mult(rbc->Curr_it->affine_trans, B, 4, 4, soln[1]);
                print_vector("vt1", 4, soln[1], "%g ");
                matrix_vec_mult(inv_affine_tran, soln[1], 4, 4, B);
                vt[1]->Printf("vt1-orig");
                print_vector("vt1---", 4, B, "%13.11g ");

                for(j=0;j<3;j++) B[j] = vt[2]->x[j];
                B[3] = 1.0;
                matrix_vec_mult(rbc->Curr_it->affine_trans, B, 4, 4, soln[2]);
                print_vector("vt2", 4, soln[2], "%g ");
                matrix_vec_mult(inv_affine_tran, soln[2], 4, 4, B);
                vt[2]->Printf("vt2-orig");
                print_vector("vt2---", 4, B, "%13.11g ");
            }
            // END::: Test code
            **/

	    (rbc->Curr_it)++;
     //   std::cout<<"IS ERROR HERE 6?"<<std::endl;
	}
    //std::cout<<"ERROR 5"<<std::endl;

	/*for(i=0;i<MAXD;i++)
	    for(j=0;j<2;j++)
		for(k=0;k<rbc->Buf_ele_num[i][j];k++)
		{
		    for(m=0;m<3;m++)
			vt[m] = &(rbc->node[rbc->Buf_ele[i][j][k].vertex[m]]);
	            tri_norm = (rbc->Buf_ele[i][j][k].out_norm.x);

		    if(NULL == rbc->Buf_ele[i][j][k].affine_trans) {
                        FT_matrix(&(rbc->Buf_ele[i][j][k].affine_trans), 4, 4,sizeof(double));
                        FT_matrix(&(rbc->Buf_ele[i][j][k].inv_affine_trans), 4, 4,sizeof(double));
                    }
                    Tri_affine_trans(vt[0]->x, vt[1]->x, vt[2]->x, tri_norm, rbc->Buf_ele[i][j][k].affine_trans);
                    inverse_matrix(rbc->Buf_ele[i][j][k].affine_trans, 4, rbc->Buf_ele[i][j][k].inv_affine_trans);
		}*/
        /* For parallelization purpose? */


        // for(i=0;i<rbc->ENum;i++)
        // {
        //     for(j=0;j<3;j++)
        //         vt[j]= &(rbc->node[rbc->ele[i].vertex[j]]); 
        //     tri_norm = rbc->ele[i].out_norm.x; 

        //     if(NULL == rbc->ele[i].affine_trans) {
        //         FT_matrix(&(rbc->ele[i].affine_trans), 4, 4,sizeof(double));
        //         FT_matrix(&(rbc->ele[i].inv_affine_trans), 4, 4,sizeof(double));
        //     }

        //     Tri_affine_trans(vt[0]->x, vt[1]->x, vt[2]->x, tri_norm, rbc->ele[i].affine_trans); 
        //     inverse_matrix(rbc->ele[i].affine_trans, 4, rbc->ele[i].inv_affine_trans);
        // }
        //std::cout<<"FREE() PROBLEM?"<<std::endl;
        //free(inv_affine_tran); free(jacobi_x_over_x_hat); 
        //std::cout<<"IT IS FREE() PROBLEM"<<std::endl;
}

void Utilities::Comp_ele_mass_matrix(RBC* rbc, CoordInfoVecs& coordInfoVecs)
{
        int           i,j,k,m;
        CVector       *vt[3];
        double        dtemp, one_3rd, affine_center[3], dtemp_max, old_center_max[3];
        int           dir_max, side_max, k_max;
        double        *tri_norm, len, area, sqrt_area;
        static double **rot_matr_rows = NULL, **translate, **inv_affine_tran = NULL, **jacobi_x_over_x_hat;
        double        B[4], affine_pt[3][4];
        double        det2, theta;

        one_3rd = 1.0/3.0;

        if(NULL == inv_affine_tran)
        {
            // FT_matrix(&rot_matr_rows, 4, 4,sizeof(double));
            // FT_matrix(&translate, 4, 4,sizeof(double));
            FT_matrix(&inv_affine_tran, 4, 4,sizeof(double));
            // FT_matrix(&jacobi_x_over_x_hat, 3, 3,sizeof(double));
        }

        rbc->Curr_it=rbc->Curr_ele->begin();
        rbc->Curr_ind_it=rbc->Curr_ele_index->begin();
        for(i=0;(rbc->Curr_it!=rbc->Curr_ele->end());i++)
        {
            if (i >= rbc->ENum){
                break;
            }
            if (rbc->Curr_it->vertex[0] >= (INT_MAX-100) || rbc->Curr_it->vertex[1] >= (INT_MAX-100) || rbc->Curr_it->vertex[2] >= (INT_MAX-100)){
                (rbc->Curr_it)++;
                continue;
            }
            if (rbc->Curr_it->vertex[0] <= (-INT_MAX+100) || rbc->Curr_it->vertex[1] <= (-INT_MAX+100) || rbc->Curr_it->vertex[2] <= (-INT_MAX+100)){
                (rbc->Curr_it)++;
                continue;
            }
            if(NULL == rbc->Curr_it->mass_matrix)
            {
                FT_matrix(&(rbc->Curr_it->mass_matrix), MAX_N_COEF, MAX_N_COEF,sizeof(double));
                FT_matrix(&(rbc->Curr_it->mass_inv), MAX_N_COEF, MAX_N_COEF,sizeof(double));
            }

            if(NULL == rbc->Curr_it->vec_vh_mass_matrix)
            {
                FT_matrix(&(rbc->Curr_it->vec_vh_mass_matrix), MAX_N_COEF*2, MAX_N_COEF*2,sizeof(double));
                FT_matrix(&(rbc->Curr_it->vec_vh_mass_inv), MAX_N_COEF*2, MAX_N_COEF*2,sizeof(double));
            }
        
            inverse_matrix(rbc->Curr_it->affine_trans, 4, inv_affine_tran);

            
            for(j=0;j<3;j++)
                vt[j] = &(rbc->node[rbc->Curr_it->vertex[j]]);

            for(j=0;j<3;j++) B[j] = vt[0]->x[j];
            B[3] = 1.0;
            matrix_vec_mult(rbc->Curr_it->affine_trans, B, 4, 4, affine_pt[0]);
            for(j=0;j<3;j++) B[j] = vt[1]->x[j];
            matrix_vec_mult(rbc->Curr_it->affine_trans, B, 4, 4, affine_pt[1]);
            for(j=0;j<3;j++) B[j] = vt[2]->x[j];
            matrix_vec_mult(rbc->Curr_it->affine_trans, B, 4, 4, affine_pt[2]);

            for(j = 0; j < 3; j++) 
                affine_center[j] = one_3rd*(affine_pt[0][j] + affine_pt[1][j] + affine_pt[2][j]); 
            area = triangle_area_3d(affine_pt[0], affine_pt[1], affine_pt[2]);
            sqrt_area = sqrt(area); 
            comp_mass_matrix(MAX_N_COEF, affine_pt[0], affine_pt[1], affine_pt[2], affine_center,
                      sqrt_area, rbc->Curr_it->mass_matrix);
            inverse_matrix(rbc->Curr_it->mass_matrix,MAX_N_COEF,rbc->Curr_it->mass_inv);

            if(NULL != rbc->Curr_it->vec_vh_mass_matrix)
            {
                comp_mass_matrix_of_vec_vh(MAX_N_COEF*2, affine_pt[0], affine_pt[1], affine_pt[2], 
                      affine_center, sqrt_area, rbc->Curr_it->vec_vh_mass_matrix);
                inverse_matrix(rbc->Curr_it->vec_vh_mass_matrix,MAX_N_COEF*2,
                               rbc->Curr_it->vec_vh_mass_inv);
            }

            (rbc->Curr_it)++;
        }

    /*    for(i=0;i<MAXD;i++)
            for(j=0;j<2;j++)
                for(k=0;k<rbc->Buf_ele_num[i][j];k++)
                {
                    if(NULL == rbc->Buf_ele[i][j][k].mass_matrix)
                    {
                        FT_matrix(&(rbc->Buf_ele[i][j][k].mass_matrix), MAX_N_COEF, MAX_N_COEF,sizeof(double));
                        FT_matrix(&(rbc->Buf_ele[i][j][k].mass_inv), MAX_N_COEF, MAX_N_COEF,sizeof(double));
                    }
                    if(NULL == rbc->Buf_ele[i][j][k].vec_vh_mass_matrix)
                    {
                        FT_matrix(&(rbc->Buf_ele[i][j][k].vec_vh_mass_matrix), MAX_N_COEF*2, MAX_N_COEF*2,sizeof(double));
                        FT_matrix(&(rbc->Buf_ele[i][j][k].vec_vh_mass_inv), MAX_N_COEF*2, MAX_N_COEF*2,sizeof(double));
                    }

                    inverse_matrix(rbc->Buf_ele[i][j][k].affine_trans, 4, inv_affine_tran);

                    for(m=0;m<3;m++)
                        vt[m] = &(rbc->node[rbc->Buf_ele[i][j][k].vertex[m]]);

                    for(m=0;m<3;m++) B[m] = vt[0]->x[m];
                    B[3] = 1.0;
                    matrix_vec_mult(rbc->Buf_ele[i][j][k].affine_trans, B, 4, 4, affine_pt[0]);
                    for(m=0;m<3;m++) B[m] = vt[1]->x[m];
                    matrix_vec_mult(rbc->Buf_ele[i][j][k].affine_trans, B, 4, 4, affine_pt[1]);
                    for(m=0;m<3;m++) B[m] = vt[2]->x[m];
                    matrix_vec_mult(rbc->Buf_ele[i][j][k].affine_trans, B, 4, 4, affine_pt[2]);

                    for(m=0;m<3;m++) 
                        affine_center[m] = one_3rd*(affine_pt[0][m] + affine_pt[1][m] + affine_pt[2][m]);
                    area = triangle_area_3d(affine_pt[0], affine_pt[1], affine_pt[2]);
                    sqrt_area = sqrt(area);
                    comp_mass_matrix(MAX_N_COEF, affine_pt[0], affine_pt[1], affine_pt[2], affine_center,
                      sqrt_area, rbc->Buf_ele[i][j][k].mass_matrix);
                    inverse_matrix(rbc->Buf_ele[i][j][k].mass_matrix, MAX_N_COEF,rbc->Buf_ele[i][j][k].mass_inv);

                    if(NULL != rbc->Buf_ele[i][j][k].vec_vh_mass_matrix)
                    {
                        comp_mass_matrix_of_vec_vh(MAX_N_COEF*2, affine_pt[0], affine_pt[1], affine_pt[2], 
                          affine_center, sqrt_area, rbc->Buf_ele[i][j][k].vec_vh_mass_matrix);
                        inverse_matrix(rbc->Buf_ele[i][j][k].vec_vh_mass_matrix, MAX_N_COEF*2,
                                       rbc->Buf_ele[i][j][k].vec_vh_mass_inv);
                    }
                }
        */
        // for(i=0;i<rbc->ENum;i++)
        // {
            
        //     for(j=0;j<3;j++)
        //         vt[j]= &(rbc->node[rbc->ele[i].vertex[j]]);

        //     if(NULL == rbc->ele[i].mass_matrix)
        //     {
        //         FT_matrix(&(rbc->ele[i].mass_matrix), MAX_N_COEF, MAX_N_COEF,sizeof(double));
        //         FT_matrix(&(rbc->ele[i].mass_inv), MAX_N_COEF, MAX_N_COEF,sizeof(double));
        //     }
        //     if(NULL == rbc->ele[i].vec_vh_mass_matrix)
        //     {
        //         FT_matrix(&(rbc->ele[i].vec_vh_mass_matrix), MAX_N_COEF*2, MAX_N_COEF*2,sizeof(double));
        //         FT_matrix(&(rbc->ele[i].vec_vh_mass_inv), MAX_N_COEF*2, MAX_N_COEF*2,sizeof(double));
        //     }

        //     inverse_matrix(rbc->ele[i].affine_trans, 4, inv_affine_tran);

        //     for(m=0;m<3;m++) B[m] = vt[0]->x[m];
        //     B[3] = 1.0;
        //     matrix_vec_mult(rbc->ele[i].affine_trans, B, 4, 4, affine_pt[0]);
        //     for(m=0;m<3;m++) B[m] = vt[1]->x[m];
        //     matrix_vec_mult(rbc->ele[i].affine_trans, B, 4, 4, affine_pt[1]);
        //     for(m=0;m<3;m++) B[m] = vt[2]->x[m];
        //     matrix_vec_mult(rbc->ele[i].affine_trans, B, 4, 4, affine_pt[2]);

        //     for(m=0;m<3;m++)
        //         affine_center[m] = one_3rd*(affine_pt[0][m] + affine_pt[1][m] + affine_pt[2][m]);
        //     area = triangle_area_3d(affine_pt[0], affine_pt[1], affine_pt[2]);
        //     sqrt_area = sqrt(area);
        //     comp_mass_matrix(MAX_N_COEF, affine_pt[0], affine_pt[1], affine_pt[2], affine_center,
        //         sqrt_area, rbc->ele[i].mass_matrix);
        //     inverse_matrix(rbc->ele[i].mass_matrix, MAX_N_COEF,rbc->ele[i].mass_inv);
        //     if(NULL != rbc->ele[i].vec_vh_mass_matrix)
        //     {
        //         comp_mass_matrix_of_vec_vh(MAX_N_COEF*2, affine_pt[0], affine_pt[1], affine_pt[2], 
        //             affine_center, sqrt_area, rbc->ele[i].vec_vh_mass_matrix);
        //         inverse_matrix(rbc->ele[i].vec_vh_mass_matrix, MAX_N_COEF*2,rbc->ele[i].vec_vh_mass_inv);
        //     }
        // }
}

void Utilities::matrix_matrix_mult(
        double    **mat,
        double    **matr,
        int      row,
        int      col,
        double    **ans)
{
        int      i, j, k;
    
        for(i = 0; i < row; i++)
        {
            //std::cout<<"i = "<<i<<std::endl;
            for(j = 0; j < col; j++)
            {
                //std::cout<<"j = "<<i<<std::endl;
                ans[i][j] = 0.0;
                for(k = 0; k < col; k++){
                   // std::cout<<"k = "<<i<<std::endl;
                    ans[i][j] += mat[i][k]*matr[k][j];
                    //std::cout<<"mat["<<i<<"]["<<k<<"] = "<<mat[i][k]<<std::endl;
                    //std::cout<<"mat["<<k<<"]["<<j<<"] = "<<mat[k][j]<<std::endl;
                }
            }
        }
}

/* comput Ax */
void Utilities::matrix_vec_mult(
        double    **mat,
        double    *vec,
        int       row,
        int       col,
        double    *ans)
{
        int      i, j;

        for(i = 0; i < row; i++)
        {
            ans[i] = 0.0;
            for(j = 0; j < col; j++)
            {
                ans[i] += mat[i][j]*vec[j];
            }
        }
}

void Utilities::soln_at_pt(
         double   *dg_poly,
         double   *crds,
         double   *cent,
         double   sqrt_area,
         double    *soln)
{
         int      i;
         double    val;

         for(i = 0; i < N_EQN; i++)
             soln[i] = 0.0;

         for(i = 0; i < MAX_N_COEF; i++)
         {
             val = vh_val_loc_div_free_basis(crds,cent,sqrt_area,i);
             soln[0] += dg_poly[i]*val;
         }
}

double Utilities::inter_integr_7_quad_flux_LDG(
        double    *vt[], // vertices of the tri after affine mapping. 
        double    *cent,
        double    area,
        double    **affine_trans, 
        double    **inv_affine_tran,
        double    *q[],       //input soln: flux q of LDG in local coords.
        double    crds[][3],  //output. quadrature pts where vec{q} is evaluated. 
        double    fluxx[][4], //output: [][# eq]
        double    fluxy[][4],
	double    fluxz[][4])
{
        int       i, j;
        double    *pcrds[3], one_3rd, sqrt_area;
        double    emid0[3], emid1[3], emid2[3], *pemid[3];
        double    g_vh[3];
        // double    soln_at_pt[4];

        // one_3rd = 1.0/3.0; 

        // for(j = 0; j < 3; j++) cent[j] = one_3rd*(vt[0][j] + vt[1][j] + vt[2][j]);

        for(i = 0; i < 3; i++)
        {
            // crds[i][0] = vt[i][0];
            // crds[i][1] = vt[i][1];
            memcpy(crds[i], vt[i], sizeof(double)*3); 
        }

        for(i = 0; i < 3; i++)
        {
            emid0[i] = 0.5*(vt[0][i] + vt[1][i]);
            emid1[i] = 0.5*(vt[1][i] + vt[2][i]);
            emid2[i] = 0.5*(vt[2][i] + vt[0][i]);
        }

        pemid[0] = emid0;
        pemid[1] = emid1;
        pemid[2] = emid2;
        // crds[3][0] = emid0[0];
        // crds[3][1] = emid0[1];
        memcpy(crds[3], emid0, sizeof(double)*3);
        // crds[4][0] = emid1[0];
        // crds[4][1] = emid1[1];
        memcpy(crds[4], emid1, sizeof(double)*3);
        // crds[5][0] = emid2[0];
        // crds[5][1] = emid2[1];
        memcpy(crds[5], emid2, sizeof(double)*3);
        // crds[6][0] = cent[0];
        // crds[6][1] = cent[1];
        memcpy(crds[6], cent, sizeof(double)*3);

        // area = triangle_area_3d(vt[0], vt[1], vt[2]);
        sqrt_area = sqrt(area);

        // Compute flux q at the point
        for(i = 0; i < 7; i++)
        {
            soln_at_pt(q[0], crds[i], cent, sqrt_area, fluxx[i]);
            soln_at_pt(q[1], crds[i], cent, sqrt_area, fluxy[i]);
        }
}


/*
 * edge_integr_q_hat_flux(): Computes \int_{edge} q_hat\cdot co-norm_{edge} v_h ds
 * */
/*
 * Main references: 
 * 1. Unified Analysis of Discontinuous Galerkin Methods. 
 *         Arnold et al. SIAM J. Numer. Anal. 39(5), 2002  
 * 2. High Order Discontinuous Galerkin Methods for Elliptic Problems on Surfaces.
 *         Antonietti et al. SIAM J. Numer. Anal. 53(2), 2015
 * 3. Analysis of the Discontinuous Galerkin Method for Elliptic Problems on Surfaces.
 *         Dedner et al. IMA J. of Numer. Anal. 33, 2013
 * 4. The Compact Discontinuous Galerkin (CDG) method for Elliptic Problems. 
 *         J. Peraire and P.-O. Persson. SIAM J. Sci. Comput. 30(4), 2008 
 *         \hat{q} = \avg{q} -C11[[u]] - C12[[q]]
 *  The first implementation uses definition in Ref. 1.
 *  The second implementation uses definition in Ref. 2.
 *
 * */
void Utilities::edge_integr_q_hat_flux(
        std::list<ELEMENT>::iterator &it_ele,
        std::list<ELEMENT>::iterator &it_adj_ele,
        double                  *prev_soln,        // prev rk-stepping u in element local crds
        double                  *prev_nb_soln,     // prev rk-stepping u in neighbor element local crds
        double                  *q_prev_soln[],    // prev rk-stepping q in element local crds
        double                  *q_prev_nb_soln[], // prev rk-stepping q in neighbor element local crds
        int                     side,
        int                     indx,
        double                  *e_q_int, // \int_{edge} q_hat\cdot n v_h ds
        // double                  u_hat[][10], // on edge u_hat, [# eq.][# quadra pt]
        double                  dt,
        int                     rk_iter,
        double                  *affine_pt[], // vertices of the tri after affine mapping
        double                  *aff_cent,
        double                  **adj_aff_time_inv_aff, 
                                 //neighbors' affine_tran times inverse of affine_tran of current ele (it_ele).
        double                  **aff_time_inv_adj_aff, 
                                 //affine_tran of current ele. times inverse of neighbors' affine_tran 
        double                  *adj_aff_cent, // center of neighbor element after affine trans. in
                                               // its own local coords.
        double                  **inv_affine_tran,
        double                  **inv_adj_affine_tran,
	RBC                     *rbc)
{
        int       i, j, k, nbside, dim = 3, positive_id; 
        double    **affine_trans = it_ele->affine_trans;
        double    **adj_ele_affine_trans = it_adj_ele->affine_trans;
        double    conorm[4], adj_ele_conorm[4], 
                  dtemp, dtemp_jump[10], dtemp_avg[10], wei, dtemp_hat[10]; 
        double    aff_conorm[4], aff_adj_ele_conorm[4], avg_conorm[4];
        double    side_vec[4], aff_side_vec[4], B[4], soln[3][4]; 
        double    length, sqrt_area, nb_sqrt_area, penalty_a, penalty_b[3]; 
                  //  penalty_a, penalty_b[3] are C11 and C12.
        double    co_n_pos[4]; // co-normal of the re-defined positive side of edge.
        double    pcrds[3][4], qcrds[4], qcrds_in_ngbr[4], soln_pt_in_ngbr[4];
        double    glb_qcrds[4], glb_qcrds_in_ngbr[4]; 
        CVector   *vt[3];
        double    tmpf[4], vh_pt, q_soln[3][4], q_nb_soln[3][4], //q_soln[dim][# eq]; 
                  q_nb_soln_pt_in_cur[4][4]; //  q_nb_soln_pt_in_cur[# eq][dim]

        double    soln_u[4], nb_soln_u[4], nb_soln_u_in_cur[4], nb_soln_pt_in_cur[4]; 
        double    u_jump[5][4]; //[# eq.][]
        // double    u_hat[10][4]; //[# quadra pt][# eq.] 
        double    q_hat[10][4][3]; // [# quadra pt][# eq.][vec component]
        // double    fluxx[4], fluxy[4], fluxz[4]; 
        // double    glb_fluxx[4], glb_fluxy[4], glb_fluxz[4]; 
        // double    nbfluxx[4], nbfluxy[4], nbfluxz[4]; 
        // double    glb_nbfluxx[4], glb_nbfluxy[4], glb_nbfluxz[4]; 
        double    flux[10][4], debug_vec[4];  //flux[#quadr][#eq]

        int       debug = NO; 
        int       COMPUT_aff_conorm_ON_FLY = NO; 

        /**
        if((it_ele->Id == 765 && indx == 0 && side == 0) 
          )
        {
            debug = YES; 
            printf("\n--------------------------------------\n"); 
            printf("-----------  Ele[%d] entered edge_integr_q_hat_flux(),"
                   " side = %d\n", it_ele->Id, side); 
        }
        **/

        for(nbside = 0; nbside < 3; nbside++)
        {
            if(it_adj_ele->adj_ele[nbside] == it_ele->Id)
                // if (it_ele->Id == 0){
                //     std::cout<<nbside<<std::endl;
                //     std::cout<<"lalala_it_adj_ele->"<<it_adj_ele->Id<<std::endl;
                //     std::cout<<it_adj_ele->adj_ele[nbside]<<std::endl;
                // }
                break;
        }
        if(nbside >= 3 || it_ele->Id == it_adj_ele->Id)
        {
            printf("ERROR: edge_integr_q_hat_flux(),"
                " adj_ele, ele do not share common edge\n");
            printf("Or Ids of eles = %d, %d\n", it_ele->Id, it_adj_ele->Id);
            //clean_up(ERROR); //20191006
        }
              
        length = it_ele->length_side[side];
        sqrt_area = sqrt(it_ele->area0);
        /*Compute weight C11 in Ref. 1 or \beta_e_h in Ref. 2*/
        // penalty_a = 1.0;  //Ref. 1.
        penalty_a = 1.0/length;  // Ref. 2
        //penalty_a = 0.0;  // alternating flux. This seems to lead instability for diffusion test!!!

        // if(debugging("Surf_Cahn_Hill"))
        //     penalty_a = 0.0; // alternating flux
   
        nb_sqrt_area = sqrt(it_adj_ele->area0);

        if(YES == COMPUT_aff_conorm_ON_FLY){
            /* Transform co-normals of it_ele and adj. ele to local crds. of the it_ele. */
            memcpy(conorm, it_ele->side_norm[side].x, sizeof(double)*3 );
            conorm[3] = 0.0;
            matrix_vec_mult(affine_trans, conorm, 4, 4, aff_conorm);

            memcpy(adj_ele_conorm, it_adj_ele->side_norm[nbside].x, sizeof(double)*3 );
            adj_ele_conorm[3] = 0.0;
            matrix_vec_mult(affine_trans, adj_ele_conorm, 4, 4, aff_adj_ele_conorm);
        }
        else
        {
            memcpy(aff_conorm, it_ele->aff_conorm[side], sizeof(double)*4); 
            memcpy(aff_adj_ele_conorm,it_ele->aff_adj_ele_conorm[side], sizeof(double)*4);
        }
 
        /*begin DEBUG*/
        if(it_ele->Id == 0)//YES == debug)
        {
        //     std::cout<<"it_ele->"<<it_ele->Id<<std::endl;
        //     std::cout<<"side = "<<side<<std::endl;
        //     std::cout<<"it_adj_ele->"<<it_adj_ele->Id<<std::endl;
        //     std::cout<<"nbside = "<<nbside<<std::endl;
        //     printf("Glb    conormal[%10.9g, %10.9g, %10.9g]\n", it_ele->side_norm[side].x[0],
        //        it_ele->side_norm[side].x[1], it_ele->side_norm[side].x[2]); 
        //     printf("Glb nb conormal[%10.9g, %10.9g, %10.9g]\n", it_adj_ele->side_norm[nbside].x[0],
        //        it_adj_ele->side_norm[nbside].x[1], it_adj_ele->side_norm[nbside].x[2]); 

        //     printf("Loc    conormal[%g, %g, %g]\n", aff_conorm[0], 
        //            aff_conorm[1], aff_conorm[2]); 
        //     printf("Loc nb conormal[%g, %g, %g]\n", aff_adj_ele_conorm[0], 
        //            aff_adj_ele_conorm[1], aff_adj_ele_conorm[2]); 
        //     printf("side length %10.9g, %10.9g\n\n",it_ele->length_side[side],
        //            it_adj_ele->length_side[nbside]); 
        //            std::cout<<"adj_ele_side_lengths: "<<it_adj_ele->length_side[0]<<" "<<it_adj_ele->length_side[1]<<" "<<it_adj_ele->length_side[2]<<std::endl;
        }
        /*end DEBUG*/

        for(j = 0; j < 3; j++)
            avg_conorm[j] = 0.5*(aff_conorm[j] - aff_adj_ele_conorm[j]); 

        dtemp = Mag3d(avg_conorm);
        for(j = 0; j < 3; j++) avg_conorm[j] = avg_conorm[j]/dtemp;

        /*Compute weight \beta in Ref. 2 or C12 in Ref. 1*/
        /*
        for(j = 0; j < 3; j++) penalty_b[j] = -0.5*avg_conorm[j]; 
        */
        positive_id = LDG_C12_vector(aff_conorm, aff_adj_ele_conorm, 
                       it_ele->Id, it_adj_ele->Id, 3, penalty_b);
        // for(j = 0; j < 3; j++) co_n_pos[j] = 2.0*penalty_b[j];
        /* This is to use alternating flux(when C11 = 0). 
         * \hat{q}\dot n =  q^+\dot n^+ - C11[[u]].
         */
        if(positive_id == it_ele->Id){
            for(j = 0; j < 3; j++) co_n_pos[j] = 2.0*penalty_b[j];
        }
        else{
            for(j = 0; j < 3; j++) co_n_pos[j] = -2.0*penalty_b[j];
        }

        /** Purely for debugging **/
        // for(j = 0; j < 3; j++) aff_conorm[j] = avg_conorm[j];
        // for(j = 0; j < 3; j++) aff_adj_ele_conorm[j] = -avg_conorm[j];
        /** END::: Purely for debugging **/

        if(YES == debug)
        {
            /*double debug_C12[4] = {0.0, 0.0, 0.0, 0.0}; 
            printf("Local-coords C12 vector[%g, %g, %g]\n", 
                   penalty_b[0], penalty_b[1], penalty_b[2]); 
            memcpy(B, penalty_b, sizeof(double)*3); B[3] = 0.0; 
            matrix_vec_mult(inv_affine_tran, B, 4, 4, debug_C12);
            print_vector("Global-coords C12 vector:", 4, debug_C12, "%10.9g "); */
        }

        /* test code */
        /**
        if(it_ele->Id == 0 && indx == 0)
        {
            printf("On ele[%d], side[%d]\n", it_ele->Id, side); 
            print_vector("aff_conorm", 3, aff_conorm, "%g "); 
            for(j = 0; j < 3; j++)
                aff_adj_ele_conorm[j] *= -1.0; 
            print_vector("-aff_adj_conorm", 3, aff_adj_ele_conorm, "%g "); 
            print_vector("avg_conorm", 3, avg_conorm, "%g "); 
            if(0 == side)
            {
                 memcpy(side_vec, it_ele->side_vector[side].x, sizeof(double)*3 ); 
                 side_vec[3] = 0.0;
                 matrix_vec_mult(affine_trans, side_vec, 4, 4, aff_side_vec);
                 print_vector("aff_side_vec", 3, aff_side_vec, "%g ");
                 for(j=0;j<3;j++) B[j] = rbc->node[it_ele->vertex[0]].x[j];
                 B[3] = 1.0;
                 matrix_vec_mult(affine_trans, B, 4, 4, soln[0]);
                 for(j=0;j<3;j++) B[j] = rbc->node[it_ele->vertex[1]].x[j];
                 matrix_vec_mult(affine_trans, B, 4, 4, soln[1]);
                 for(j=0;j<3;j++) B[j] = rbc->node[it_ele->vertex[2]].x[j];
                 matrix_vec_mult(affine_trans, B, 4, 4, soln[2]);
                 print_vector("vt0", 4, soln[0], "%g ");
                 print_vector("vt1", 4, soln[1], "%g ");
                 print_vector("vt2", 4, soln[2], "%g ");
              
            }
            printf("\n"); 
        }
        **/
        /* END:: test code */

        for(i = 0; i < N_EQN; i++) e_q_int[i] = 0.0;

        // BEGIN:: use saved flux value. 
        // std::cout<<"it_ele->qhat_f_flag[side] = "<<it_ele->qhat_f_flag[side]<<std::endl;
        if(true == it_ele->qhat_f_flag[side]){
            for(k = 0; k < Gauss_N; k++){
                for(i = 0; i < dim; i++)
                    qcrds[i] = (affine_pt[(side+1)%3][i] + affine_pt[side][i])/2.0 +
                           (affine_pt[(side+1)%3][i] - affine_pt[side][i])/2.0*q[k];    
                vh_pt = vh_val_loc_div_free_basis(qcrds, aff_cent, sqrt_area, indx); 
 
                // for(i = 0; i < N_EQN; i++)
                //     flux[k][i] = it_ele->q_hat_flux_store[side][k][i]; 

                // flux[#quadr][#eq]
                for(i = 0; i < N_EQN; i++) {
                    // tmpf[i] = flux[k][i]*vh_pt*qw[k];
                    if (it_ele->Id == 5){
                        // std::cout<<"q_hat_flux "<<it_ele->q_hat_flux_store[side][k][i]<<" "<<vh_pt<<" "<<qw[k]<<std::endl;
                        // std::cout<<"tmpf["<<i<<"] = "<<it_ele->q_hat_flux_store[side][k][i]*vh_pt*qw[k]<<std::endl;
                    }
                    tmpf[i] = it_ele->q_hat_flux_store[side][k][i]*vh_pt*qw[k];
                    e_q_int[i] += tmpf[i];
                }
            }
            if(1 == Gauss_N) {
                for(i = 0; i < N_EQN+1; i++) e_q_int[i] *= length;
            }
            else {
                for(i = 0; i < N_EQN; i++) e_q_int[i] *= (0.5*length);
            }
            return; 
        } 
        // END:: use saved flux value. 
        
        for(k = 0; k < Gauss_N; k++)
        {
            for(i = 0; i < dim; i++){ // dim = 3
                qcrds[i] = (affine_pt[(side+1)%3][i] + affine_pt[side][i])/2.0 +
                           (affine_pt[(side+1)%3][i] - affine_pt[side][i])/2.0*q[k];
            }

            vh_pt = vh_val_loc_div_free_basis(qcrds, aff_cent, sqrt_area, indx);

            soln_at_pt(prev_soln, qcrds, aff_cent, sqrt_area, soln_u);

            soln_at_pt(q_prev_soln[0], qcrds, aff_cent, sqrt_area, q_soln[0]);
            soln_at_pt(q_prev_soln[1], qcrds, aff_cent, sqrt_area, q_soln[1]);
            // std::cout<<q_soln[0][0]<<" "<< q_soln[1][0]<<std::endl;
            //std::cout<<" "<<std::endl;
            for(i = 0; i < N_EQN; i++) q_soln[2][i] = 0.0; 
           
            qcrds[3] = 1.0; 
            matrix_vec_mult(adj_aff_time_inv_aff, qcrds, 4, 4, qcrds_in_ngbr); 
            matrix_vec_mult(inv_affine_tran, qcrds, 4, 4, glb_qcrds); 

            /*begin debug*/
            if(YES == debug)
            {
                printf("Gauss_quadr[%d] affine_crds[%g %g %g %g]\n",
                            k, qcrds[0], qcrds[1], qcrds[2],  qcrds[3]); 
                printf("Gauss_quadr[%d] in glb_crds[%g %g %g %g]\n",
                            k, glb_qcrds[0], glb_qcrds[1], glb_qcrds[2], glb_qcrds[3]); 
                printf("Gauss_quadr[%d] in neigh_loc_crds[%g %g %g %g]\n",
                          k, qcrds_in_ngbr[0], qcrds_in_ngbr[1], qcrds_in_ngbr[2],
                          qcrds_in_ngbr[3]); 
                printf("Ele local uh = %14.12g\n", soln_u[0]); 
                printf("Ele local q vector = (%14.12g, %14.12g)\n", q_soln[0][0], q_soln[1][0]); 
                // memcpy(B, qcrds, sizeof(double)*2); 
                // B[2] = soln_u[0]; B[3] = 1.0; 
                // matrix_vec_mult(inv_affine_tran, B, 4, 4, debug_vec);
                // print_vector("uh in global-coords:", 4, debug_vec, "%10.9g ");  
            }
            /*END::: begin debug*/

            soln_at_pt(prev_nb_soln, qcrds_in_ngbr, adj_aff_cent, nb_sqrt_area, nb_soln_u);

            soln_at_pt(q_prev_nb_soln[0], qcrds_in_ngbr, adj_aff_cent, nb_sqrt_area, q_nb_soln[0]);
            soln_at_pt(q_prev_nb_soln[1], qcrds_in_ngbr, adj_aff_cent, nb_sqrt_area, q_nb_soln[1]);

            /*** 
            // code for advection
            if(debugging("linear_adv0"))
            {
                flux_from_soln_glb_pt(glb_qcrds, NULL, glb_fluxx, glb_fluxy, glb_fluxz);
                for(i = 0; i < N_EQN; i++)
                {
                    tmpflux[0] = glb_fluxx[i];
                    tmpflux[1] = glb_fluxy[i];
                    tmpflux[2] = glb_fluxz[i];
                    tmpflux[3] = 0.0;
                    matrix_vec_mult(affine_trans, tmpflux, 4, 4, tmpflux2);
                    fluxx[i] = tmpflux2[0];
                    fluxy[i] = tmpflux2[1];
                    fluxz[i] = tmpflux2[2];
                }
            }
            ***/

            /*begin debug*/
            /*
            if(YES == debug)
            {
                printf("           glb_fluxx[%g, %g, %g]\n", glb_fluxx[0], glb_fluxy[0], glb_fluxz[0]); 
                printf("   affined glb_fluxx[%g, %g, %g]\n", fluxx[0], fluxy[0], fluxz[0]); 
            }
            */
            /*END::: begin debug*/

            /*transform q_flux soln (vector) point from local coords defined on adj element
             * to local coords defined on current element*/
            soln_pt_in_ngbr[0] = q_nb_soln[0][0]; 
            soln_pt_in_ngbr[1] = q_nb_soln[1][0]; 
            soln_pt_in_ngbr[2] = 0.0; 
            soln_pt_in_ngbr[3] = 0.0; 
            matrix_vec_mult(aff_time_inv_adj_aff, soln_pt_in_ngbr, 4, 4, q_nb_soln_pt_in_cur[0]); 

            if(YES == debug)
            {
                printf("Adj. Ele local-crds:     q vector = (%14.12g, %14.12g)\n", 
                  q_nb_soln[0][0], q_nb_soln[1][0]);
                printf("Adj. Ele in-it_ele crds: q vector = (%14.12g, %14.12g, %14.12g, %g)\n", 
                  q_nb_soln_pt_in_cur[0][0], q_nb_soln_pt_in_cur[0][1], 
                             q_nb_soln_pt_in_cur[0][2], q_nb_soln_pt_in_cur[0][3]);
            }

            /**
            soln_pt_in_ngbr[0] = 0.0;
            soln_pt_in_ngbr[1] = q_nb_soln[1][0]; 
            soln_pt_in_ngbr[2] = 0.0; 
            soln_pt_in_ngbr[3] = 0.0; 
            matrix_vec_mult(aff_time_inv_adj_aff, soln_pt_in_ngbr, 4, 4, q_nb_soln_pt_in_cur[1]); 
            **/

            /*Transform soln point from local coords defined on adj element
             * to coords defined on current element (it_ele)*/
            /* 
             * ??????????????????????
             * nb_soln_u[0] is the height of the graph of the soln defined 
             * on local coords of adj. element.
             * When use it, do we need to tranform it into coords. of it_ele??????
             * 
             */
            memcpy(soln_pt_in_ngbr, qcrds_in_ngbr, sizeof(double)*2);
            soln_pt_in_ngbr[2] = nb_soln_u[0];
            soln_pt_in_ngbr[3] = 1.0;
            matrix_vec_mult(aff_time_inv_adj_aff, soln_pt_in_ngbr, 4, 4, nb_soln_pt_in_cur);
            nb_soln_u_in_cur[0] = nb_soln_pt_in_cur[2];

            /*begin debug*/
            if(YES == debug)
            {
                qcrds_in_ngbr[3] = 1.0;
                // matrix_vec_mult(inv_adj_affine_tran,  qcrds_in_ngbr, 4, 4, glb_qcrds_in_ngbr);
                // printf("Gauss_quadr[%d] in nb affine_crds[%g %g %g %g]\n",
                //        k, qcrds_in_ngbr[0], qcrds_in_ngbr[1], qcrds_in_ngbr[2], qcrds_in_ngbr[3]); 
                // printf("Gauss_quadr[%d] in nb glb_crds[%g %g %g %g]\n",
                //        k, glb_qcrds_in_ngbr[0], glb_qcrds_in_ngbr[1], glb_qcrds_in_ngbr[2],
                //             glb_qcrds_in_ngbr[3]); 
                printf("nb_soln in adj. ele local coords = %14.12g, in it_ele local coords = %14.12g\n",
                                nb_soln_u[0], nb_soln_u_in_cur[0]);

                double debug_vec_adj_uh[4];
                memcpy(B, qcrds_in_ngbr, sizeof(double)*2); 
                B[2] = nb_soln_u[0]; B[3] = 1.0; 
                matrix_vec_mult(inv_adj_affine_tran, B, 4, 4, debug_vec_adj_uh);
                // print_vector("Adj. uh in global-coords:", 4, debug_vec_adj_uh, "%10.9g ");  

                double debug_vec_uh[4], tmp_u_jump[4]; 
                memcpy(B, qcrds, sizeof(double)*2); 
                B[2] = soln_u[0]; B[3] = 1.0; 
                matrix_vec_mult(inv_affine_tran, B, 4, 4, debug_vec_uh);
                // print_vector("uh in global-coords:", 4, debug_vec_uh, "%10.9g ");  

                tmp_u_jump[0] = penalty_a*(debug_vec_uh[2]*it_ele->side_norm[side].x[0] +
                                           debug_vec_adj_uh[2]*it_adj_ele->side_norm[nbside].x[0]); 
                tmp_u_jump[1] = penalty_a*(debug_vec_uh[2]*it_ele->side_norm[side].x[1] +
                                           debug_vec_adj_uh[2]*it_adj_ele->side_norm[nbside].x[1]); 
                tmp_u_jump[2] = penalty_a*(debug_vec_uh[2]*it_ele->side_norm[side].x[2] +
                                           debug_vec_adj_uh[2]*it_adj_ele->side_norm[nbside].x[2]); 
                // print_vector("Test in Global-coords tmp_u_jump:", 4, tmp_u_jump, "%10.9g ");
            }
            /*END::: begin debug*/

            /***
                // tmp_alpha = PI; 
                tmp_alpha = 0.0; 
                for(i = 0; i < N_EQN; i++)
                {
                    q_flux[k][i] = 0.5*(fluxx[i]*avg_conorm[0]+fluxy[i]*avg_conorm[1] + 
                                      fluxz[i]*avg_conorm[2] +
                                      nbfluxx[i]*avg_conorm[0]+nbfluxy[i]*avg_conorm[1] + 
                                      nbfluxz[i]*avg_conorm[2] -
                                 tmp_alpha*(nb_soln_u_in_cur[i]-q_soln[i]));
                }
            ***/

            /*Ref. 1: Unified Analysis of Discontinuous Galerkin Methods.*/
            /*Need to compute fluxes: \hat q*/
            /*Compute u jump*/
            /**
            for(i = 0; i < N_EQN; i++)
            {
                // u_jump[i][0] = penalty_a*(nb_soln_u_in_cur[i]*aff_adj_ele_conorm[0] + 
                //                           soln_u[i]*aff_conorm[0]);
                // u_jump[i][1] = penalty_a*(nb_soln_u_in_cur[i]*aff_adj_ele_conorm[1] + 
                //                           soln_u[i]*aff_conorm[1]);
                // u_jump[i][2] = penalty_a*(nb_soln_u_in_cur[i]*aff_adj_ele_conorm[2] + 
                //                           soln_u[i]*0.0); 
                u_jump[i][0] = penalty_a*(nb_soln_u[i]*aff_adj_ele_conorm[0] + 
                                          soln_u[i]*aff_conorm[0]);
                u_jump[i][1] = penalty_a*(nb_soln_u[i]*aff_adj_ele_conorm[1] + 
                                          soln_u[i]*aff_conorm[1]);
                u_jump[i][2] = penalty_a*(nb_soln_u[i]*aff_adj_ele_conorm[2] + 
                                          soln_u[i]*0.0); 
            }

            if(YES == debug)
            {
                printf("---u jump in local-coords is[%g, %g, %g]\n", u_jump[0][0],
                    u_jump[0][1], u_jump[0][2]); 
                memcpy(B, u_jump[0], sizeof(double)*3); B[3] = 0.0;
                matrix_vec_mult(inv_affine_tran, B, 4, 4, debug_vec);
                print_vector("Global-coords u_jump:", 4, debug_vec, "%10.9g ");
            }

            // q_soln & q_nb_soln_pt_in_cur needs to be changed for multiple eqs. 
            for(i = 0; i < N_EQN; i++)
            {
                q_hat[k][i][0] = 0.5*(q_nb_soln_pt_in_cur[i][0] + q_soln[0][i]); 
                q_hat[k][i][1] = 0.5*(q_nb_soln_pt_in_cur[i][1] + q_soln[1][i]); 
                q_hat[k][i][2] = 0.5*(q_nb_soln_pt_in_cur[i][2] + 0.0); 
            }
            if(YES == debug)
            {
                printf("\n q_avg (%g %g %g)\n ", 
                      q_hat[k][0][0], q_hat[k][0][1], q_hat[k][0][2]); 
                memcpy(B, q_hat[k][0], sizeof(double)*3); B[3] = 0.0;
                matrix_vec_mult(inv_affine_tran, B, 4, 4, debug_vec);
                print_vector("Global-coords q_avg:", 4, debug_vec, "%10.9g ");
            }

            //aff_adj_ele_conorm, aff_conorm
            for(i = 0; i < N_EQN; i++)
            {
                dtemp[i] = q_nb_soln_pt_in_cur[i][0]*aff_adj_ele_conorm[0] + 
                           q_nb_soln_pt_in_cur[i][1]*aff_adj_ele_conorm[1] + 
                           q_nb_soln_pt_in_cur[i][2]*aff_adj_ele_conorm[2] + 
                           q_soln[0][i]*aff_conorm[0] + q_soln[1][i]*aff_conorm[1];
            }
            if(YES == debug)
            {
                printf("\n q_jump (%10.9g)\n ", dtemp[0]); 
            }

            for(i = 0; i < N_EQN; i++)
            {
                q_hat[k][i][0] += (-dtemp[i]*penalty_b[0] - penalty_a*u_jump[i][0]);
                q_hat[k][i][1] += (-dtemp[i]*penalty_b[1] - penalty_a*u_jump[i][1]);
                q_hat[k][i][2] += (-dtemp[i]*penalty_b[2] - penalty_a*u_jump[i][2]);
            }
            if(YES == debug)
            {
                printf("\n Local-coords q_flux (%g %g %g)\n ", 
                         q_hat[k][0][0], q_hat[k][0][1], q_hat[k][0][2]); 
                memcpy(B, q_hat[k][0], sizeof(double)*3); B[3] = 0.0;
                matrix_vec_mult(inv_affine_tran, B, 4, 4, debug_vec);
                print_vector("Global-coords q_flux", 4, debug_vec, "%10.9g ");
            }
                
            // for(i = 0; i < N_EQN; i++)
            //     u_hat[i][k] = 0.5*(nb_soln_u_in_cur[i] + soln_u[i]) + 
            //                        Dot3d(penalty_b, u_jump[i]);

            for(i = 0; i < N_EQN; i++)
            {
                flux[k][i] = Dot3d(q_hat[k][i], avg_conorm); 
                tmpf[i] = flux[k][i]*vh_pt*qw[k];
                e_q_int[i] += tmpf[i];
            }

            if(YES == debug)
            {
                printf("Ele[%d]-edge(%d) q[%d] q_hat dot n = %10.9g\n", 
                       it_ele->Id, side, k, flux[k][0]); 
            }
            **/
            /*******END::: Ref. 1*******/
            /***************************/

            /*Ref. 2: High Order Discontinuous Galerkin Methods for
             *  Elliptic Problems on Surfaces. Antonietti.*/
            /*Need to compute fluxes: \hat q*/
            /*Compute u jump*/
            for(i = 0; i < N_EQN; i++)
                u_jump[i][0] = (soln_u[i] - nb_soln_u[i]);

            if(YES == debug){
                printf("---u jump in it_ele local-coords is[%14.12g]\n", u_jump[0][0]);
            }

            //aff_adj_ele_conorm, aff_conorm
            // Both q_soln & q_nb_soln_pt_in_cur are defined in it_ele coords.
            // q_soln & q_nb_soln_pt_in_cur needs to be changed for multiple eqs. 
            // q_soln[dim][# eq]; q_nb_soln_pt_in_cur[# eq][dim]
            // first compute q average.
            for(i = 0; i < N_EQN; i++)
            {
                dtemp_avg[i] = -(q_nb_soln_pt_in_cur[i][0]*aff_adj_ele_conorm[0]  
                               + q_nb_soln_pt_in_cur[i][1]*aff_adj_ele_conorm[1]  
                               + q_nb_soln_pt_in_cur[i][2]*aff_adj_ele_conorm[2]) + 
                               (q_soln[0][i]*aff_conorm[0] + q_soln[1][i]*aff_conorm[1]);
                dtemp_avg[i] *= 0.5; 
            }
            if(YES == debug) {
                double temp_val; 
                temp_val = q_soln[0][0]*aff_conorm[0] + q_soln[1][0]*aff_conorm[1];
                temp_val += -q_nb_soln[0][0]*adj_ele_conorm[0] - q_nb_soln[1][0]*adj_ele_conorm[1];
                printf("\n it_ele crds: q_avg (%14.12g), test-val %14.12g\n ",
                    dtemp_avg[0], temp_val); 
            }

            //Now compute q jump
            for(i = 0; i < N_EQN; i++)
            {
                dtemp_jump[i] = q_nb_soln_pt_in_cur[i][0]*aff_adj_ele_conorm[0] + 
                                q_nb_soln_pt_in_cur[i][1]*aff_adj_ele_conorm[1] + 
                                q_nb_soln_pt_in_cur[i][2]*aff_adj_ele_conorm[2] + 
                           q_soln[0][i]*aff_conorm[0] + q_soln[1][i]*aff_conorm[1];
            }
            /**
            if(YES == debug)
            {
                printf("\n q_jump (%10.9g) dot_q_nb(%13.12g) dot_q_n(%13.12g)\n ", 
                         dtemp_jump[0], Dot3d(q_nb_soln_pt_in_cur[0], aff_adj_ele_conorm ), 
                          q_soln[0][i]*aff_conorm[0] + q_soln[1][i]*aff_conorm[1]); 
            }
            **/

            wei = Dot3d(penalty_b,co_n_pos);
            for(i = 0; i < N_EQN; i++){
                dtemp_hat[i] = (dtemp_avg[i] + dtemp_jump[i]*wei - penalty_a*u_jump[i][0]);
                // std::cout<<"dtemp_avg "<<dtemp_avg[i]<<std::endl;;
                // std::cout<<"dtemp_jump & wei "<<dtemp_jump[i]<<" , "<<wei<<std::endl;;
                // std::cout<<"penalty_a & u_jump "<<penalty_a<<" , "<<u_jump[i][0]<<std::endl;;
            }

            if(YES == debug) {
                printf("\n it_ele local-coords q_flux (%14.12g), C12 = %g\n ", 
                    dtemp_hat[0], wei); 
            }

            for(i = 0; i < N_EQN; i++)
                for(j = 0; j < 3; j++)
                    q_hat[k][i][j] = dtemp_hat[i]*aff_conorm[j];

            /**
            // TMP: try alternating flux
            if(positive_id == it_ele->Id){
                for(i = 0; i < N_EQN; i++)
                    for(j = 0; j < 3; j++)
                        q_hat[k][i][j] = q_soln[j][i];
            }
            else{
                for(i = 0; i < N_EQN; i++)
                    for(j = 0; j < 3; j++)
                        q_hat[k][i][j] = q_nb_soln_pt_in_cur[i][j];
            }
            // END::: TMP: try alternating flux
            **/

            // flux[#quadr][#eq]
            for(i = 0; i < N_EQN; i++)
            {
                flux[k][i] = Dot3d(q_hat[k][i], aff_conorm); 
                // std::cout<<"q_hat "<<q_hat[k][i][0]<<" "<<q_hat[k][i][1]<<" "<<q_hat[k][i][2]<<std::endl;
                // std::cout<<"aff_conorm "<<aff_conorm[0]<<" "<<aff_conorm[1]<<" "<<aff_conorm[2]<<std::endl;
                // flux[k][i] = Dot3d(q_hat[k][i], avg_conorm); 
                tmpf[i] = flux[k][i]*vh_pt*qw[k];
                e_q_int[i] += tmpf[i];
            }

            if(YES == debug)
            {
                printf("Ele[%d]-edge(%d) q[%d]"
                  " q_hat dot n = %13.12g q_hat[%13.12g, %13.12g, %13.12g]\n", 
                       it_ele->Id, side, k, flux[k][0], 
                       q_hat[k][0][0], q_hat[k][0][1], q_hat[k][0][2]); 
            }
            /*******END::: Ref. 2*******/

        }//END::: for(k = 0; k < Gauss_N; k++)

        for(i = 0; i < N_EQN; i++) e_q_int[i] *= (0.5*length);

        if(YES == debug)
            printf("Ele[%d]-edge(%d) q_hat dot n = %13.12g\n",
                     it_ele->Id, side, e_q_int[0]); 

        /**
        if((it_ele->Id == 765 && side == 0) ||
           (it_ele->Id == 766 && side == 1)) {
            printf("Ele[%d]-edge(%d) indx %d; q_hat dot n = %13.12g\n",
                     it_ele->Id, side, indx, e_q_int[0]); 
        }
        **/

        if(false == it_ele->qhat_f_flag[side]){
            it_ele->qhat_f_flag[side]= true; 
            it_adj_ele->qhat_f_flag[nbside]= true; 
 
            for(k = 0; k < Gauss_N; k++){
                for(i = 0; i < N_EQN; i++)
                    it_ele->q_hat_flux_store[side][k][i] = flux[k][i]; 
            }
            for(k = 0; k < Gauss_N; k++){
                for(i = 0; i < N_EQN; i++)
                    it_adj_ele->q_hat_flux_store[nbside][Gauss_N-1-k][i] = -flux[k][i]; 
            }
        }
}

/* edge_integr_new_u_hat_flux() computes: 
 * \int_e \hat{u_h} co-norm_{e} \cdot \vec{r} ds,
 * \vec{r} is the vector test function.  
 * For diffusion problem, it pairs with edge_integr_q_hat_flux()
 * to give alternating flux. 
 */
/*
 * Main references: 
 * 1. Unified Analysis of Discontinuous Galerkin Methods. 
 *         Arnold et al. SIAM J. Numer. Anal. 39(5), 2002  
 * 2. High Order Discontinuous Galerkin Methods for Elliptic Problems on Surfaces.
 *         Antonietti et al. SIAM J. Numer. Anal. 53(2), 2015
 * 3. Analysis of the Discontinuous Galerkin Method for Elliptic Problems on Surfaces.
 *         Dedner et al. IMA J. of Numer. Anal. 33, 2013
 * 4. The Compact Discontinuous Galerkin (CDG) method for Elliptic Problems. 
 *         J. Peraire and P.-O. Persson. SIAM J. Sci. Comput. 30(4), 2008 
 *         \hat{u} = \avg{u} + C12 \cdot [[u]].
 *  The first implementation uses definition in Ref. 1. 
 *  The second implementation uses definition in Ref. 2.        
 * */
void  Utilities::edge_integr_new_u_hat_flux(
        std::list<ELEMENT>::iterator &it_ele,
        std::list<ELEMENT>::iterator &it_adj_ele,
        double                  *cur_soln, // in element local crds
        double                  *cur_nb_soln, // in neighbor element local crds
        int                     side,
        int                     indx,     
        double                  *e_u_int, // \int_{edge} u_hat n \cdot (v_h) ds
        double                  dt,
        int                     rk_iter,
        double                  *affine_pt[], // vertices of the tri after affine mapping
        double                  *aff_cent,
        double                  **adj_aff_time_inv_aff, 
                                 //neighbors' affine_tran times inverse of affine_tran of current
        double                  **aff_time_inv_adj_aff, 
                                 //affine_tran of current times inverse of neighbors' affine_tran 
        double                  *adj_aff_cent, // center of neighbor element after affine trans.
        double                  **inv_affine_tran,
        double                  **inv_adj_affine_tran,
	RBC                     *rbc,
        bool                    *it_f_flag,
        bool                    *it_adj_f_flag,
        double                  it_hat_flux_store[][4][2],
        double                  it_adj_hat_flux_store[][4][2])
{
        int       i, j, k, nbside, dim = 3, positive_id; 
        double    **affine_trans = it_ele->affine_trans;
        double    **adj_ele_affine_trans = it_adj_ele->affine_trans;
        double    conorm[4], adj_ele_conorm[4], dtemp; 
        double    aff_conorm[4], aff_adj_ele_conorm[4], avg_conorm[4] = {0.0, 0.0, 0.0, 0.0};
        double    side_vec[4], aff_side_vec[4], B[4], soln[3][4]; 
        double    length, sqrt_area, nb_sqrt_area, penalty_a, penalty_b[3]; 
                  //  penalty_a, penalty_b[3] are C11 and C12.
        double    co_n_pos[4]; // co-normal of the re-defined positive side of edge. 
        double    pcrds[3][4], qcrds[4], qcrds_in_ngbr[4], soln_pt_in_ngbr[4];
        double    glb_qcrds[4], glb_qcrds_in_ngbr[4]; 
        CVector   *vt[3];
        double    tmpf[4], vec_vh_pt[3] = {0.0, 0.0, 0.0}; 
        // double    q_soln[3][4], q_nb_soln[3][4], q_nb_soln_pt_in_cur[4][4];
                  // q_nb_soln_pt_in_cur[# eq][dim]
        double    soln_u[4], nb_soln_u[4], nb_soln_u_in_cur[4], nb_soln_pt_in_cur[4]; 
        double    u_jump[4][4]; //[# eq.][]
        double    u_hat[4][10]; //[# eq.][# quadra pt] 
        double    flux[10][4], debug_vec[4]; 

        int       debug = NO; 
        int       COMPUT_aff_conorm_ON_FLY = NO; 

        /**
        if((it_ele->Id == 0 && indx == 0 && side == 0)||
           (it_ele->Id == 2 && indx == 0 && side == 1))
        {
            debug = YES; 
            printf("\n\n&&&&&&&&&&&&&&&&\n"
                   " ele[%d] entered edge_integr_new_u_hat_flux(),"
                   " side = %d, index = %d\n", it_ele->Id, side, indx); 
            printf("&&&&&&&&&&&&&&&&&&&&&&\n"); 
        }
        **/

        /*
        // if(it_ele->Id == 0 && indx == 0)
        if(it_ele->Id == 0)
        {
            debug = YES; 
            printf("\n\n********ele[%d] entered edge_integr_new_u_hat_flux(),"
                   " side = %d, index = %d\n", it_ele->Id, side, indx); 
        }
        */

        for(nbside = 0; nbside < 3; nbside++)
        {
            if(it_adj_ele->adj_ele[nbside] == it_ele->Id)
                break;
        }
        if(nbside >= 3 || it_ele->Id == it_adj_ele->Id)
        {
            printf("ERROR: edge_integr_new_u_hat_flux(),"
                   " adj_ele, ele do not share common edge\n");
            printf("Or Ids of eles = %d, %d\n", it_ele->Id, it_adj_ele->Id);
            //clean_up(ERROR); //20191006
        }
              
        length = it_ele->length_side[side];
        sqrt_area = sqrt(it_ele->area0);
        // penalty_a = 10.0/length; 
        // penalty_a = 1.0; 
        penalty_a = 0.0; 
   
        nb_sqrt_area = sqrt(it_adj_ele->area0);

        if(YES == COMPUT_aff_conorm_ON_FLY){
            memcpy(conorm, it_ele->side_norm[side].x, sizeof(double)*3 );
            conorm[3] = 0.0;
            matrix_vec_mult(affine_trans, conorm, 4, 4, aff_conorm);

            memcpy(adj_ele_conorm, it_adj_ele->side_norm[nbside].x, sizeof(double)*3 );
            adj_ele_conorm[3] = 0.0;
            matrix_vec_mult(affine_trans, adj_ele_conorm, 4, 4, aff_adj_ele_conorm);
        }
        else
        {
            memcpy(aff_conorm, it_ele->aff_conorm[side], sizeof(double)*4);
            memcpy(aff_adj_ele_conorm,it_ele->aff_adj_ele_conorm[side], sizeof(double)*4);
        }

        /*begin DEBUG*/
        /*
        if(YES == debug)
        {
            printf("Glb    conormal[%g, %g, %g]\n", it_ele->side_norm[side].x[0],
               it_ele->side_norm[side].x[1], it_ele->side_norm[side].x[2]); 
            printf("Glb nb conormal[%g, %g, %g]\n", it_adj_ele->side_norm[nbside].x[0],
               it_adj_ele->side_norm[nbside].x[1], it_adj_ele->side_norm[nbside].x[2]); 

            printf("Loc    conormal[%g, %g, %g]\n", aff_conorm[0], 
                   aff_conorm[1], aff_conorm[2]); 
            printf("Loc nb conormal[%g, %g, %g]\n", aff_adj_ele_conorm[0], 
                   aff_adj_ele_conorm[1], aff_adj_ele_conorm[2]); 
            printf("side length %g, %g\n\n",it_ele->length_side[side],
                   it_adj_ele->length_side[nbside]); 
        }
        */
        /*end DEBUG*/

        for(j = 0; j < 3; j++)
            avg_conorm[j] = 0.5*(aff_conorm[j] - aff_adj_ele_conorm[j]);    
        dtemp = Mag3d(avg_conorm);
        for(j = 0; j < 3; j++)
            avg_conorm[j] = avg_conorm[j]/dtemp;

        /*
        for(j = 0; j < 3; j++) penalty_b[j] = -0.5*avg_conorm[j]; 
        */
        positive_id = LDG_C12_vector(aff_conorm, aff_adj_ele_conorm, 
                       it_ele->Id, it_adj_ele->Id, 3, penalty_b);
        /*
         * This is to use alternating flux. \hat{u} = u^-.
         */
        if(positive_id == it_ele->Id) {
            for(j = 0; j < 3; j++) co_n_pos[j] = 2.0*penalty_b[j]; 
        }
        else {
            for(j = 0; j < 3; j++) co_n_pos[j] = -2.0*penalty_b[j]; 
        }

        /** Purely for debugging **/
        // for(j = 0; j < 3; j++) aff_conorm[j] = avg_conorm[j];
        // for(j = 0; j < 3; j++) aff_adj_ele_conorm[j] = -avg_conorm[j];
        /** END::: Purely for debugging **/

        if(YES == debug)
        {
            /*double debug_C12[4] = {0.0, 0.0, 0.0, 0.0};
            printf("Local-coords C12 vector[%g, %g, %g], co-nor-pos [%g, %g, %g]\n",
                   penalty_b[0], penalty_b[1], penalty_b[2],
                   co_n_pos[0], co_n_pos[1], co_n_pos[2]);
            print_vector("aff_conorm in local:", 4, aff_conorm, "%10.9g ");
            memcpy(B, penalty_b, sizeof(double)*3); B[3] = 0.0;
            matrix_vec_mult(inv_affine_tran, B, 4, 4, debug_C12);
            print_vector("Global-coords C12 vector:", 4, debug_C12, "%10.9g ");
            memcpy(B, co_n_pos, sizeof(double)*3); B[3] = 0.0;
            matrix_vec_mult(inv_affine_tran, B, 4, 4, debug_C12);
            print_vector("Global-coords co_n_pos vector:", 4, debug_C12, "%10.9g ");
            */
        }

        for(i = 0; i < N_EQN; i++) e_u_int[i] = 0.0;

        // BEGIN:: use saved flux value.
        // if(true == it_ele->uhat_f_flag[side]){
        if(true == it_f_flag[side]){
            for(k = 0; k < Gauss_N; k++){
                for(i = 0; i < dim; i++)
                    qcrds[i] = (affine_pt[(side+1)%3][i] + affine_pt[side][i])/2.0 +
                           (affine_pt[(side+1)%3][i] - affine_pt[side][i])/2.0*q[k];
         
                vec_vh_val_loc_div_free_basis(qcrds, aff_cent, sqrt_area, indx, vec_vh_pt);

                for(i = 0; i < N_EQN; i++) {
                    flux[k][i] = Dot3d(vec_vh_pt, aff_conorm);
                    // tmpf[i] = flux[k][i]*it_ele->u_hat_flux_store[side][k][i]*qw[k];
                    tmpf[i] = flux[k][i]*it_hat_flux_store[side][k][i]*qw[k];
                    e_u_int[i] += tmpf[i];
                }
            }
            if(1 == Gauss_N) {
                for(i = 0; i < N_EQN+1; i++) e_u_int[i] *= length;
            }
            else {
                for(i = 0; i < N_EQN; i++) e_u_int[i] *= (length*0.5); 
            }
            return; 
        }
        // END:: use saved flux value.

        for(k = 0; k < Gauss_N; k++)
        {
            for(i = 0; i < dim; i++)
                qcrds[i] = (affine_pt[(side+1)%3][i] + affine_pt[side][i])/2.0 +
                           (affine_pt[(side+1)%3][i] - affine_pt[side][i])/2.0*q[k];

            vec_vh_val_loc_div_free_basis(qcrds, aff_cent, sqrt_area, indx, vec_vh_pt);

            soln_at_pt(cur_soln, qcrds, aff_cent, sqrt_area, soln_u);

            qcrds[3] = 1.0; 
            matrix_vec_mult(adj_aff_time_inv_aff, qcrds, 4, 4, qcrds_in_ngbr); 
            matrix_vec_mult(inv_affine_tran, qcrds, 4, 4, glb_qcrds); 

            /*begin debug*/
            /*
            if(YES == debug)
            {
                printf("Gauss_quadr[%d] affine_crds[%g %g %g]\n",
                            k, qcrds[0], qcrds[1], qcrds[2]); 
                printf("Gauss_quadr[%d] in glb_crds[%g %g %g]\n",
                            k, glb_qcrds[0], glb_qcrds[1], glb_qcrds[2]); 
            }
            */
            /*END::: begin debug*/

            soln_at_pt(cur_nb_soln, qcrds_in_ngbr, adj_aff_cent, nb_sqrt_area, nb_soln_u);

            /*transform soln point from local coords defined on adj element
             * to coords defined on current element (it_ele)*/
            /*
             * ??????????????????????
             * nb_soln_u[0] is the height of the graph of the soln 
             * defined on local coords of adj. element.
             * When use it, do we need to tranform it into coords. of it_ele??????
             *
             */
            memcpy(soln_pt_in_ngbr, qcrds_in_ngbr, sizeof(double)*2);
            soln_pt_in_ngbr[2] = nb_soln_u[0];
            soln_pt_in_ngbr[3] = 1.0;
            matrix_vec_mult(aff_time_inv_adj_aff, soln_pt_in_ngbr, 4, 4, nb_soln_pt_in_cur);
            nb_soln_u_in_cur[0] = nb_soln_pt_in_cur[2];

            /*begin debug*/
            /*
            if(YES == debug)
            {
                qcrds_in_ngbr[3] = 1.0;
                matrix_vec_mult(inv_adj_affine_tran,  qcrds_in_ngbr, 4, 4, glb_qcrds_in_ngbr);
                printf("Gauss_quadr[%d] in neigh affine_crds[%g %g %g]\n",
                            k, qcrds_in_ngbr[0], qcrds_in_ngbr[1], qcrds_in_ngbr[2]); 
                printf("Gauss_quadr[%d] in neigh glb_crds[%g %g %g]\n",
                            k, glb_qcrds_in_ngbr[0], glb_qcrds_in_ngbr[1], glb_qcrds_in_ngbr[2]); 
            }
            */
            /*END::: begin debug*/

            /*Ref. 1: Unified Analysis of Discontinuous Galerkin Methods.*/
            /*Need to compute fluxes: \hat u and \hat q*/
            /*Compute u jump*/
            /**
            for(i = 0; i < N_EQN; i++)
            {
                // u_jump[i][0] = penalty_a*(nb_soln_u_in_cur[i]*aff_adj_ele_conorm[0] + 
                //                           soln_u[i]*aff_conorm[0]);
                // u_jump[i][1] = penalty_a*(nb_soln_u_in_cur[i]*aff_adj_ele_conorm[1] + 
                //                           soln_u[i]*aff_conorm[1]);
                // u_jump[i][2] = penalty_a*(nb_soln_u_in_cur[i]*aff_adj_ele_conorm[2] + 
                //                           soln_u[i]*0.0); 
                u_jump[i][0] = penalty_a*(nb_soln_u[i]*aff_adj_ele_conorm[0] +
                                          soln_u[i]*aff_conorm[0]);
                u_jump[i][1] = penalty_a*(nb_soln_u[i]*aff_adj_ele_conorm[1] +
                                          soln_u[i]*aff_conorm[1]);
                u_jump[i][2] = penalty_a*(nb_soln_u[i]*aff_adj_ele_conorm[2] +
                                          soln_u[i]*0.0);
            }

            for(i = 0; i < N_EQN; i++)
            {
                // u_hat[k][i] = 0.5*(nb_soln_u_in_cur[i] + soln_u[i]) 
                //               - Dot3d(penalty_b, u_jump[i]);
                u_hat[k][i] = 0.5*(nb_soln_u[i] + soln_u[i]) 
                              - Dot3d(penalty_b, u_jump[i]);
            } 
            if(YES == debug)
            {
                printf("---u jump in local-coords is[%g, %g, %g]\n", u_jump[0][0],
                    u_jump[0][1], u_jump[0][2]);
                memcpy(B, u_jump[0], sizeof(double)*3); B[3] = 0.0;
                matrix_vec_mult(inv_affine_tran, B, 4, 4, debug_vec);
                print_vector("Global-coords u_jump:", 4, debug_vec, "%10.9g ");
                printf("----u_hat in local-coords = %10.9g\n", u_hat[k][0]); 
                printf("----u_avg in local-coords = %10.9g\n",
                       0.5*(nb_soln_u_in_cur[0] + soln_u[0])); 
                printf("----u(%g, %g) in both local-coords\n",
                       nb_soln_u[0], soln_u[0]); 
            }
            for(i = 0; i < N_EQN; i++)
            {
                flux[k][i] = Dot3d(vec_vh_pt, avg_conorm); 
                tmpf[i] = flux[k][i]*u_hat[k][i]*qw[k];
                e_u_int[i] += tmpf[i];
            }

            if(YES == debug)
            {
                printf("Edge u_hat_flux[at %d] = %10.9g (). n dot r = %10.9g\n", 
                      k, tmpf[0], Dot3d(vec_vh_pt, avg_conorm)); 
            }
            **/
            /*******END::: Ref. 1*******/
            /***************************/
            /*Ref. 2: High Order Discontinuous Galerkin Methods for 
             *  Elliptic Problems on Surfaces..*/
            /*Compute u jump*/
            for(i = 0; i < N_EQN; i++)
                u_jump[i][0] = (soln_u[i] - nb_soln_u[i]); 

            //for diffusion problem. It alternates direction with \hat{q}.
            for(i = 0; i < N_EQN; i++)
            {
                    // u_hat[k][i] = 0.5*(nb_soln_u_in_cur[i] + soln_u[i]) 
                    //               - u_jump[i][0]*Dot3d(penalty_b,co_n_pos);
                u_hat[k][i] = 0.5*(nb_soln_u[i] + soln_u[i]) 
                              - u_jump[i][0]*Dot3d(penalty_b,co_n_pos);
            }

            /*begin debug*/
            if(YES == debug)
            {
                printf("k(%d)----u jump in local-coords is[%13.12g], C12*n = %13.12g\n", 
                    k, u_jump[0][0], Dot3d(penalty_b,co_n_pos));
                printf("k(%d)----u_hat in local-coords = %14.12g\n", k, u_hat[k][0]); 
                // printf("----u_avg(%g, %g) in K^+ local-coords = %10.9g\n",
                //        nb_soln_u_in_cur[0], soln_u[0], 
                //        0.5*(nb_soln_u_in_cur[0] + soln_u[0])); 
                printf("k(%d)----u(%14.12g) nb_u(%14.12g) in both local-coords\n",
                       k, soln_u[0], nb_soln_u[0]); 
            }
            /*END::: begin debug*/

            for(i = 0; i < N_EQN; i++)
            {
                flux[k][i] = Dot3d(vec_vh_pt, aff_conorm); 
                // flux[k][i] = Dot3d(vec_vh_pt, avg_conorm); 
                tmpf[i] = flux[k][i]*u_hat[k][i]*qw[k];
                e_u_int[i] += tmpf[i];
            }

            if(YES == debug)
            {
                printf("Edge u_hat_flux[at %d] = %13.12g (). n dot r = %13.12g, indx=%d\n", 
                      k, tmpf[0], Dot3d(vec_vh_pt, aff_conorm), indx); 
            }
            /*******END::: Ref. 2*******/

        }//END::: for(k = 0; k < Gauss_N; k++)

        for(i = 0; i < N_EQN; i++) e_u_int[i] *= (0.5*length);

        if(YES == debug)
            printf("Edge side[%d] u_hat_flux = %13.12g, indx=%d\n", side, e_u_int[0], indx); 

        /**
        if((it_ele->Id == 765 && side == 0) ||
           (it_ele->Id == 766 && side == 1)) {
            for(k = 0; k < Gauss_N; k++)
            printf("Ele[%d]-edge(%d) indx %d; quadrature[%d] u_hat = %13.12g\n",
                     it_ele->Id, side, indx, k, u_hat[k][0]);
        }
        **/

        // if(false == it_ele->uhat_f_flag[side]){
        if(false == it_f_flag[side]){
            // it_ele->uhat_f_flag[side] = true;
            // it_adj_ele->uhat_f_flag[nbside] = true;
            it_f_flag[side] = true;
            it_adj_f_flag[nbside] = true;

            for(k = 0; k < Gauss_N; k++){
                for(i = 0; i < N_EQN; i++)
                    // it_ele->u_hat_flux_store[side][k][i] = u_hat[k][i];
                    it_hat_flux_store[side][k][i] = u_hat[k][i];
            }
            for(k = 0; k < Gauss_N; k++){
                for(i = 0; i < N_EQN; i++)
                    // it_adj_ele->u_hat_flux_store[nbside][Gauss_N-1-k][i] = u_hat[k][i];
                    it_adj_hat_flux_store[nbside][Gauss_N-1-k][i] = u_hat[k][i];
            }
        }
}

void Utilities::inverse_matrix(
        double      **mat,
        int         size,
        double      **inv)
{              
                                 
        int        i, j, ii, jj, k;
        static int size_ = 0, *indx, improve = 3;
        static long double *col = NULL, **a, *colcp, *r;
        long double d, sdp;
        long double **iden;

        if(NULL == inv) return;

        return inverse_matrix_gj(mat, size, inv);

}

// inverse by Gauss-Jordan
void Utilities:: inverse_matrix_gj(
        double      **mat,
        int        size,
        double      **inv)
{
        static int size_ = 0;
        static long double **b = NULL, **a;
        int i, j;
        long double **tmpm, **iden;

        if(size == 1)
        {
            inv[0][0] = 1.0/mat[0][0];
            return;
        }

        if(size != size_)
        {
            if(b != NULL)
            {
                //free(a); free(b);
            }
            size_ = size;
            FT_matrix(&a, (size+1), (size+1), sizeof(long double));
            FT_matrix(&b, (size+1), (size+1), sizeof(long double));
        }
        for(j = 0; j < size; j++)
        {
            for(i = 0; i < size; i++)
                a[i+1][j+1] = mat[i][j];
        }

        for(i = 1; i <= size; i++)
        {
            b[i][1] = 0.0;
            b[i][2] = 0.0;
        }
        b[1][1] = 1.0;
        b[2][2] = 1.0;

        gaussj(a,size,b,2);

        for(j = 0; j < size; j++)
        {
            for(i = 0; i < size; i++)
                inv[i][j] = a[i+1][j+1];
        }

}

// gauss-Jordan elimination with pivoting
void Utilities::gaussj(
        long double **a,
        int         n,
        long double **b,
        int         m)
{
        static int size_ = 0,  *indxc = NULL, *indxr, *ipiv;
        int   i, icol, irow, j, k , l, ll;
        long double big, dum, pivinv, temp;

        if(n != size_)
        {
            if(indxc != NULL)
            {
                //free(indxc); free(indxr); free(ipiv);
            }
            FT_vector(&indxc, (n+1),  sizeof(int));
            FT_vector(&indxr, (n+1),  sizeof(int));
            FT_vector(&ipiv, (n+1),  sizeof(int));
            size_ = n;
        }

        for(i = 1; i <= n; i++) ipiv[i] = 0;
        for(i = 1; i <= n; i++)
        {
            big = 0.0;
            for(j = 1; j <= n; j++)
            {
                if(ipiv[j] != 1)
                {
                    for(k =1; k<=n; k++)
                    {
                        if(ipiv[k] == 0)
                        {
                            if(fabsl(a[j][k]) >= big)
                            {
                                big = fabsl(a[j][k]);
                                irow = j;
                                icol = k;
                            }
                        }
                    }
                }
            }
            ++(ipiv[icol]);

            if(irow != icol)
            {
                for(l = 1; l <= n; l++) SWAP(a[irow][l], a[icol][l]);
                for(l = 1; l <= m; l++) SWAP(b[irow][l], b[icol][l]);
            }
            indxr[i] = irow;
            indxc[i] = icol;
            if(a[icol][icol] == 0.0)
            {
                printf("ERROR() gaussj, main ele = 0.0 at icol %d\n", icol);
                //clean_up(ERROR); //20191006
            }
            pivinv = 1.0/a[icol][icol];
            a[icol][icol] = 1.0;
            for(l = 1; l <= n; l++) a[icol][l] *= pivinv;
            for(l = 1; l <= m; l++) b[icol][l] *= pivinv;
            for(ll= 1; ll <= n; ll++)
            {
                if(ll != icol)
                {
                    dum = a[ll][icol];
                    a[ll][icol] = 0.0;
                    for(l = 1; l <= n; l++) a[ll][l] -= a[icol][l]*dum;
                    for(l = 1; l <= m; l++) b[ll][l] -= b[icol][l]*dum;
                }
            }
        }

        for(l =n; l >=1; l--)
        {
            if(indxr[l] != indxc[l])
            {
                for(k =1; k <= n; k++)
                    SWAP(a[k][indxr[l]], a[k][indxc[l]]);
            }
        }
}

double Utilities::triangle_area_3d(
       double   *p1,
       double   *p2,
       double   *p3)
{
       double a;
       double alpha;
       double area;
       double b;
       double base;
       double c;
       double dot;
       double height;

       dot =
        (p2[0] - p1[0]) * (p3[0] - p1[0]) +
        (p2[1] - p1[1]) * (p3[1] - p1[1]) +
        (p2[2] - p1[2]) * (p3[2] - p1[2]);

       base = enorm0_3d(p1, p2);
       if (base == 0.0)
       {
           height = 0.0;
       }
       else
       {
           alpha = dot / ( base * base );
           a = p3[0] - p1[0] - alpha * (p2[0] - p1[0]);
           b = p3[1] - p1[1] - alpha * (p2[1] - p1[1]);
           c = p3[2] - p1[2] - alpha * (p2[2] - p1[2]);

           height = sqrt ( a * a + b * b + c * c );
       }
       area = 0.5 * base * height;
       return area;
}

void Utilities::comp_mass_matrix(
        int      n_coeff,
        double   *pt0,
        double   *pt1,
        double   *pt2,
        double   *cent,
        double   sqrt_area,
        double   **mass_m)
{
        int      i, j, k;
        // double   *pcrds[3], midpt[3][2], tmp[2];
        long double det, area;
        double   crds[20][2];
        long  double tmpans[20];
        static double w[16] ={0.144315607677787,0.095091634267285,0.095091634267285,0.095091634267285,
                             0.103217370534718, 0.103217370534718,0.103217370534718,
                             0.032458497623198,0.032458497623198,0.032458497623198,
                             0.027230314174435,0.027230314174435,0.027230314174435,
                             0.027230314174435,0.027230314174435,0.027230314174435};
        double   tmp_area;

        tri_quadrature_16_pts(pt0, pt1, pt2, crds);

        det = (long double)(pt1[0]-pt0[0])*(pt2[1]-pt0[1]) -
              (long double)(pt2[0]-pt0[0])*(pt1[1]-pt0[1]);
        area = det*0.5;

        if(MAX_N_COEF == 1)
        {
            mass_m[0][0] = det*0.5;
        }
        else if(MAX_N_COEF == 3 || MAX_N_COEF == 6 || MAX_N_COEF == 10 ||  MAX_N_COEF == 15)
        {
            for(k = 0; k < MAX_N_COEF; k++)
            {
                for(i = 0; i < MAX_N_COEF; i++)
                {
                    for(j = 0; j < 16; j++)
                        tmpans[j] = normalized_val(crds, cent, sqrt_area, j, i)*
                                normalized_val(crds, cent, sqrt_area, j, k);
                    mass_m[k][i] = 0.0;
                    for(j = 0; j < 16; j++)
                        mass_m[k][i] += tmpans[j]*w[j];
                    mass_m[k][i] *= area;
                }
            }
        }
        else
        {
            printf("ERROR: implement comp_mass_matrix for MAX_N_COEF = %d\n",
                MAX_N_COEF);
            //clean_up(ERROR); //20191006
        }
}


void Utilities::comp_mass_matrix_of_vec_vh(
        int      n_coeff,
        double   *pt0,
        double   *pt1,
        double   *pt2,
        double   *cent,
        double   sqrt_area,
        double   **mass_m)
{
        int      i, j, k;
        long double det, area;
        double   crds[20][2];
        long  double tmpans[20];
        static double w[16] ={0.144315607677787,0.095091634267285,0.095091634267285,
                              0.095091634267285,
                             0.103217370534718, 0.103217370534718,0.103217370534718,
                             0.032458497623198,0.032458497623198,0.032458497623198,
                             0.027230314174435,0.027230314174435,0.027230314174435,
                             0.027230314174435,0.027230314174435,0.027230314174435};
        double   tmp_area;
        double   vec1[3], vec2[3]; 

        tri_quadrature_16_pts(pt0, pt1, pt2, crds);

        det = (long double)(pt1[0]-pt0[0])*(pt2[1]-pt0[1]) -
              (long double)(pt2[0]-pt0[0])*(pt1[1]-pt0[1]);
        area = det*0.5;

        if(MAX_N_COEF == 1)
        {
            mass_m[0][0] = det*0.5; mass_m[0][1] = 0.0;
            mass_m[1][0] = 0.0;     mass_m[1][1] = det*0.5;
        }
        else if(MAX_N_COEF == 3 || MAX_N_COEF == 6 || 
                MAX_N_COEF == 10 ||  MAX_N_COEF == 15)
        {
            for(k = 0; k < MAX_N_COEF*2; k++)
            {
                for(i = 0; i < MAX_N_COEF*2; i++)
                {
                    for(j = 0; j < 16; j++)
                    {
                        // tmpans[j] = normalized_val(crds, cent, sqrt_area, j, i)*
                        //         normalized_val(crds, cent, sqrt_area, j, k);
                        normalized_vec_vh_val(crds, cent, sqrt_area, j, i, vec1); 
                        normalized_vec_vh_val(crds, cent, sqrt_area, j, k, vec2); 
                        tmpans[j] = Dot2d(vec1, vec2);
                    }
                    mass_m[k][i] = 0.0;
                    for(j = 0; j < 16; j++)
                        mass_m[k][i] += tmpans[j]*w[j];
                    mass_m[k][i] *= area;
                }
            }
        }
        else
        {
            printf("ERROR: implement comp_mass_matrix_of_vec_vh for MAX_N_COEF = %d\n",
                MAX_N_COEF);
            //clean_up(ERROR); //20191006
        }
}

double Utilities::vh_val_loc_div_free_basis(
        double *crds,
        double *cent,
        double sqrt_area,
        int    indx)
{
        double ans;
        switch(indx)
        {
        case 0:
            ans = 1.0;
        break;
        case 1:
            ans = (crds[0]-cent[0])/sqrt_area;
        break;
        case 2:
            ans = (crds[1]-cent[1])/sqrt_area;
        break;
        case 3:
            ans = sqr((crds[0]-cent[0])/sqrt_area);
        break;
        case 4:
            ans = ((crds[0]-cent[0])/sqrt_area)*((crds[1]-cent[1])/sqrt_area);
        break;
        case 5:
            ans = sqr((crds[1]-cent[1])/sqrt_area);
        break;
        case 6:
            ans = cub((crds[0]-cent[0])/sqrt_area);
        break;
        case 7:
            ans = sqr((crds[0]-cent[0])/sqrt_area)*((crds[1]-cent[1])/sqrt_area);
        break;
        case 8:
            ans = ((crds[0]-cent[0])/sqrt_area)*sqr((crds[1]-cent[1])/sqrt_area);
        break;
        case 9:
            ans = cub((crds[1]-cent[1])/sqrt_area);
        break;
        default:
            printf("ERROR vh_val_loc_div_free_basis, implement 2D degree %d\n", indx);
            //clean_up(ERROR); //20191006
        }
        return ans;
}

int Utilities::LDG_C12_vector(
	double *n_pos,
	double *n_neg,
	int    id_pos,
	int    id_neg,
        int    dim, 
	double *C12)
{
        double S1, S2; 
        int    i, positive; 

	if(id_pos > id_neg)
        {
            S1 = 1.0;
            S2 = 0.0;
            positive = id_pos;
        }
        else
        {
            S1 = 0.0;
            S2 = 1.0;
            positive = id_neg;
        }

        for(i = 0; i < dim; i++)
            C12[i] = 0.5*(S1*n_pos[i] + S2*n_neg[i]);
        return positive;
}

void Utilities::vec_vh_val_loc_div_free_basis(
        double  *crds,
        double  *cent,
        double  sqrt_area,
        int     indx,
        double  *ans)
{
        switch(indx)
        {
        case 0:
            ans[0] = 1.0; ans[1] = 0.0;
        break;
        case 1:
            ans[0] = 0.0; ans[1] = 1.0;
        break;
        case 2:
            ans[0] = (crds[0]-cent[0])/sqrt_area; ans[1] = 0.0;
        break;
        case 3:
            ans[1] = (crds[0]-cent[0])/sqrt_area; ans[0] = 0.0;
        break;

        case 4:
            ans[0] = (crds[1]-cent[1])/sqrt_area; ans[1] = 0.0;
        break;
        case 5:
            ans[1] = (crds[1]-cent[1])/sqrt_area; ans[0] = 0.0;
        break;

        case 6:
            ans[0] = sqr((crds[0]-cent[0])/sqrt_area); ans[1] = 0.0;
        break;
        case 7:
            ans[1] = sqr((crds[0]-cent[0])/sqrt_area); ans[0] = 0.0;
        break;

        case 8:
            ans[0] = ((crds[0]-cent[0])/sqrt_area)*((crds[1]-cent[1])/sqrt_area); ans[1] = 0.0;
        break;
        case 9:
            ans[1] = ((crds[0]-cent[0])/sqrt_area)*((crds[1]-cent[1])/sqrt_area); ans[0] = 0.0;
        break;

        case 10:
            ans[0] = sqr((crds[1]-cent[1])/sqrt_area); ans[1] = 0.0;
        break;
        case 11:
            ans[1] = sqr((crds[1]-cent[1])/sqrt_area); ans[0] = 0.0;
        break;

        /**
        case 6:
            ans[0] = cub((crds[0]-cent[0])/sqrt_area);
        break;
        case 7:
            ans[0] = sqr((crds[0]-cent[0])/sqrt_area)*((crds[1]-cent[1])/sqrt_area);
        break;
        case 8:
            ans[0] = ((crds[0]-cent[0])/sqrt_area)*sqr((crds[1]-cent[1])/sqrt_area);
        break;
        case 9:
            ans[0] = cub((crds[1]-cent[1])/sqrt_area);
        break;
        **/
        default:
            printf("ERROR vec_vh_val_loc_div_free_basis, implement 2D degree %d\n", indx);
            //clean_up(ERROR); //20191006
        }
}

double Utilities::enorm0_3d(
       double   *p0,
       double   *p1)
{
       double value;

       value = sqrt(
         sqr(p1[0] - p0[0])  +
         sqr(p1[1] - p0[1])  +
         sqr(p1[2] - p0[2]));

       return value;
}

void Utilities::tri_quadrature_16_pts(
        double       *pcrds0,
        double       *pcrds1,
        double       *pcrds2,
        double       crds[][2])
{
        static double xg[16] = {0.33333333333333333,0.081414823414554,0.459292588292723,
                               0.459292588292723,0.65886138449648,0.17056930775176,
                               0.17056930775176,0.898905543365938,0.050547228317031,
                               0.050547228317031,0.008394777409958001,0.263112829634638,
                               0.728492392955404,0.263112829634638,0.728492392955404,
                               0.008394777409958001};
        static double yg[16] = {0.33333333333333333,0.459292588292723,0.459292588292723,
                               0.081414823414554,0.17056930775176,0.17056930775176 ,
                               0.65886138449648,0.050547228317031,0.050547228317031,
                               0.898905543365938,0.263112829634638,0.728492392955404,
                               0.008394777409958001,0.008394777409958001,0.263112829634638,
                               0.728492392955404};

        int i;
        for(i = 0; i < 16; i++)
        {
            crds[i][0] = pcrds0[0] + (pcrds1[0]-pcrds0[0])*xg[i] + (pcrds2[0]-pcrds0[0])*yg[i];
            crds[i][1] = pcrds0[1] + (pcrds1[1]-pcrds0[1])*xg[i] + (pcrds2[1]-pcrds0[1])*yg[i];
        }
}

double Utilities::normalized_val(
        double   crds[][2],
        double   *cent,
        double   sqrt_area,
        int      pos,
        int      indx)
{
        long double tmpx, tmpy;

        switch(indx)
        {
        case 0:
            return 1.0;
        break;
        case 1:
            tmpx = crds[pos][0] - cent[0];
            return tmpx/sqrt_area;
        break;
        case 2:
            tmpy = crds[pos][1] - cent[1];
            return tmpy/sqrt_area;
        break;
        case 3:
            tmpx = crds[pos][0] - cent[0];
            return sqr(tmpx/sqrt_area);
        break;
        case 4:
            tmpx = (crds[pos][0] - cent[0])/sqrt_area;
            tmpy = (crds[pos][1] - cent[1])/sqrt_area;
            return tmpx*tmpy;
        break;
        case 5:
            tmpy = crds[pos][1] - cent[1];
            return sqr(tmpy/sqrt_area);
        break;
        case 6:
            tmpx = crds[pos][0] - cent[0];
            return cub(tmpx/sqrt_area);
        break;
        case 7:
            tmpx = (crds[pos][0] - cent[0])/sqrt_area;
            tmpy = (crds[pos][1] - cent[1])/sqrt_area;
            return (sqr(tmpx)*tmpy);
        break;
        case 8:
            tmpx = (crds[pos][0] - cent[0])/sqrt_area;
            tmpy = (crds[pos][1] - cent[1])/sqrt_area;
            return (tmpx*sqr(tmpy));
        break;
        case 9:
            tmpy = (crds[pos][1] - cent[1])/sqrt_area;
            return cub(tmpy);
        break;
        case 10:
             tmpx = (crds[pos][0] - cent[0])/sqrt_area;
             return sqr(tmpx)*sqr(tmpx);
        break;
        case 11:
             tmpx = (crds[pos][0] - cent[0])/sqrt_area;
             tmpy = (crds[pos][1] - cent[1])/sqrt_area;
             return cub(tmpx)*tmpy;
        break;
        case 12:
             tmpx = (crds[pos][0] - cent[0])/sqrt_area;
             tmpy = (crds[pos][1] - cent[1])/sqrt_area;
             return sqr(tmpx)*sqr(tmpy);
        break;
        case 13:
             tmpx = (crds[pos][0] - cent[0])/sqrt_area;
             tmpy = (crds[pos][1] - cent[1])/sqrt_area;
             return (tmpx)*cub(tmpy);
        break;
        case 14:
             tmpy = (crds[pos][1] - cent[1])/sqrt_area;
             return sqr(tmpy)*sqr(tmpy);
        break;
        }

        printf("ERROR: normalized_val()\n");
        //clean_up(ERROR); //20191006
}

void Utilities::normalized_vec_vh_val(
        double  crds[][2],
        double  *cent,
        double  sqrt_area,
        int     pos,
        int     indx,
        double  *ans)
{
        switch(indx)
        {
        case 0:
            ans[0] = 1.0; ans[1] = 0.0;
        break;
        case 1:
            ans[0] = 0.0; ans[1] = 1.0;
        break;
        case 2:
            ans[0] = (crds[pos][0]-cent[0])/sqrt_area; ans[1] = 0.0;
        break;
        case 3:
            ans[1] = (crds[pos][0]-cent[0])/sqrt_area; ans[0] = 0.0;
        break;

        case 4:
            ans[0] = (crds[pos][1]-cent[1])/sqrt_area; ans[1] = 0.0;
        break;
        case 5:
            ans[1] = (crds[pos][1]-cent[1])/sqrt_area; ans[0] = 0.0;
        break;

        case 6:
            ans[0] = sqr((crds[pos][0]-cent[0])/sqrt_area); ans[1] = 0.0;
        break;
        case 7:
            ans[1] = sqr((crds[pos][0]-cent[0])/sqrt_area); ans[0] = 0.0;
        break;

        case 8:
            ans[0] = ((crds[pos][0]-cent[0])/sqrt_area)*((crds[pos][1]-cent[1])/sqrt_area); ans[1] = 0.0;
        break;
        case 9:
            ans[1] = ((crds[pos][0]-cent[0])/sqrt_area)*((crds[pos][1]-cent[1])/sqrt_area); ans[0] = 0.0;
        break;

        case 10:
            ans[0] = sqr((crds[pos][1]-cent[1])/sqrt_area); ans[1] = 0.0;
        break;
        case 11:
            ans[1] = sqr((crds[pos][1]-cent[1])/sqrt_area); ans[0] = 0.0;
        break;
        default:
            printf("ERROR normalized_vec_vh_val, implement 2D degree %d\n", indx);
            //clean_up(ERROR); //20191006
        }

}



#if !defined(min)
#define	    min(a,b)     (((a) < (b)) ? (a) : (b))
#endif /* !defined(min) */
typedef  void	    *POINTER;	/* Pointer to Unknown Data Types. */

void Utilities::LDG_Surf_diffusion(
        //POINTER *RBC_pointer, //will throw in the address of the rbc struct
        RBC *rbc,
        RBC *n_rbc,
        double sqrt_D,
        double k_0,
        double k_1,
        double k_2,
        double k_3,
        double k_4,
        double k_ss,
        double beta,
        double gamma,
        double q1,
        double h,
        // double sink_coef,
        // double source_coef,
        double fr_time,
        double time_step_factor,
        //int        *iperm,
        //double     *dh,
        //Wave       *wv,
        //Wave       *nwv,
        //Front      *fr,
        //Front      *nfr,
        double     dt,
        int        INCL_ADV
        )
{
        //RBC           *rbc = (RBC *)fr->RBC_pointer[0];
        //POINTER       *RBC_pointers = fr->RBC_pointer;
        //RBC           *n_rbc = (RBC *)nfr->RBC_pointer[0];
        //POINTER       *n_RBC_pointers = nfr->RBC_pointer;
        int           Num_RBCs = 1, i_NumRBC, i;
        std::list<ELEMENT>::iterator it_ele, n_it_ele;
        std::list<ELEMENT>::iterator it_adj_ele[3], n_it_adj_ele[3];
        // list<ELEMENT>::iterator ***sv_it_adj_ele, ***sv_n_it_adj_ele;
        int           rk_iter, side;
        ELEMENT       adj_ele[3]; 
        static int    first = YES;
        double        max_dt; 
        /*Gas_param     **prmslst;*/ //20191006
        int           num_params; 
        double        sqrt_diff_ceof; // this is used as diffusion coefficient in
                                  // current implementation.
        double        time, min_edge_len = HUGE_VAL;  

        // DEBUG_ENTER(LDG_Surf_diffusion)
    //std::cout<<"ERROR 1"<<std::endl;
        // printf("\n\n\nEntered LDG_Surf_diffusion(), dt = %g, step %d\n", dt, fr->step); 
        // fflush(stdout);
        // printf("*********************************\n\n"); 

        if(YES == first)
        {
            first = NO;
            if(1 == Gauss_N)
            {
                q[0] = 0.0; qw[0] = 2.0; 
            }
            else if(2 == Gauss_N)
            {
                q[0] = -1.0/sqrt(3.0); q[1] = 1.0/sqrt(3.0);
                qw[0] = 1.0; qw[1] = 1.0;
            }
            else if(3 == Gauss_N)
            {
                q[0] = -sqrt(0.6); q[1] = 0.0; q[2] = sqrt(0.6);
                qw[0] = 5.0/9.0; qw[1] = 8.0/9.0; qw[2] = 5.0/9.0;
            }
            else if(Gauss_N == 4)
            {
                q[0] = -0.86113631159405257522; q[1] = -0.33998104358485626480;
                q[2] =  0.33998104358485626480; q[3] = 0.86113631159405257522;
                qw[0] =  0.34785484513745385737; qw[1] = 0.65214515486254614263;
                qw[2] =  0.65214515486254614263; qw[3] = 0.34785484513745385737;
            }
            else
            {
                q[0] = 0.0;
                qw[0] = 2.0;
                Gauss_N = 1; 
            }
        }
//std::cout<<"ERROR 2"<<std::endl;
        //num_params = return_params_list(&prmslst);
        // for(i = 0; i < num_params; i++) fprint_Gas_param(stdout,prmslst[i]);
        sqrt_diff_ceof = sqrt_D;
        //sqrt_diff_ceof = sqrt(prmslst[0]->shear_visc); 

        // printf("WRNING: LDG_Surf_diffusion(), sqrt_diff_ceof = %g\n", sqrt_diff_ceof);
        // clean_up(0); 

        //newdt = HUGE_VAL; //20191006, not sure why newdt is never defined in the source code.
//std::cout<<"ERROR 3"<<std::endl;
        for(rk_iter = 0; rk_iter < RK_STEP; rk_iter++)
        {
            if(1 == RK_STEP) {
                time = fr_time;
                //time = fr->time; 
            }
            else if(2 == RK_STEP) {
                if(rk_iter == 0)
                    time = fr_time;
                    // time = fr->time; 
                else
                    time = fr_time + dt;
                    // time = fr->time + dt; 
            }
            else {
                if(YES ==  INCL_ADV) {
                     printf("ERROR: LDG_Surf_diffusion(). need to set time\n");
                     //clean_up(ERROR); //20191006
                }
            }
//std::cout<<"ERROR 4"<<std::endl;
            for(i_NumRBC = 0; i_NumRBC < Num_RBCs; i_NumRBC++) {
                //n_rbc = (RBC*)n_RBC_pointers;
                //n_rbc = (RBC*)n_RBC_pointers[i_NumRBC];
                min_edge_len = HUGE_VAL;
//std::cout<<"ERROR 5"<<std::endl;
                // printf("\n\n---STEP %d------rk_iter = %d, rbc = %p, new rbc = %p\n\n",
                //                 fr->step, rk_iter, rbc, n_rbc);

                it_ele = rbc->Curr_ele->begin(); 
                n_it_ele = n_rbc->Curr_ele->begin(); 
               // std::cout<<"ERROR 6"<<std::endl;
                for (int k = 0;(it_ele != rbc->Curr_ele->end()); it_ele++, k++)
                {
                    if (k >= rbc->ENum){
                        break;
                    }
                    it_ele->uhat_f_flag[0] = it_ele->uhat_f_flag[1] 
                                = it_ele->uhat_f_flag[2] = false; 
                    it_ele->qhat_f_flag[0] = it_ele->qhat_f_flag[1] 
                                = it_ele->qhat_f_flag[2] = false; 
                }
//std::cout<<"ERROR 7"<<std::endl;
                it_ele = rbc->Curr_ele->begin(); 
                n_it_ele = n_rbc->Curr_ele->begin(); 
//std::cout<<"ERROR 8"<<std::endl;
                for(i=0;(it_ele != rbc->Curr_ele->end());i++, it_ele++, n_it_ele++)
                {
                    if (i >= rbc->ENum){
                        break;
                    }
                    if (it_ele->vertex[0] >= (INT_MAX-100) || it_ele->vertex[1] >= (INT_MAX-100) || it_ele->vertex[2] >= (INT_MAX-100)){
                        continue;
                    }
                    if (it_ele->vertex[0] <= (-INT_MAX+100) || it_ele->vertex[1] <= (-INT_MAX+100) || it_ele->vertex[2] <= (-INT_MAX+100)){
                        continue;
                    }
                    //std::cout<<"ERROR 8.1"<<std::endl;
                    if(RK_STEP == (rk_iter + 1))
                    {
                        for(side = 0; side < 3; side++)
                        {
                            if(it_ele->length_side[side] < min_edge_len)
                                min_edge_len = it_ele->length_side[side]; 
                        }
                    }
                   //std::cout<<"ERROR 8.2"<<std::endl;
                    for(side = 0; side < 3; side++){
                        it_adj_ele[side] = rbc->sv_it_adj_ele[it_ele->Id][side];
                        // if (it_ele->Id %200 == 0){
                        //     std::cout<<"double-checking if sv_it_adj_ele is functioning correctly (Id: "<< it_ele->Id<<"): "<<rbc->sv_it_adj_ele[it_ele->Id][side]->Id<<std::endl;
                        // }
                    }
 
                   //std::cout<<"ERROR 8.3"<<std::endl; 
                    for(side = 0; side < 3; side++){
                        n_it_adj_ele[side] = n_rbc->sv_it_adj_ele[n_it_ele->Id][side];
                        // if (it_ele->Id % 200 == 0){
                        //     std::cout<<"double-checking if sv_it_adj_ele is functioning correctly (for n_it_adj, Id: "<< it_ele->Id<<") : "<<rbc->sv_it_adj_ele[n_it_ele->Id][side]->Id<<std::endl;
                        // }
                        //std::cout<<"n_it_adj_ele"<<n_it_adj_ele[side]<<std::endl;
                        }

                   // std::cout<<"ERROR 8.4"<<std::endl;
                    for(side = 0; side < 3; side++) {
                        if(n_it_adj_ele[side]->Id != it_adj_ele[side]->Id) {
                          //  std::cout<<"made it in the loop?"<<std::endl;
                            printf("ERROR: LDG_Surf_diffusion().\n");
                            printf("side %d. New and old ele list not"
                                   " consistent(1-3) %d, %d\n",
                                    side, n_it_adj_ele[side]->Id, it_adj_ele[side]->Id);
                           // clean_up(ERROR); //20191006
                        }
                    }

                    //std::cout<<"ERROR 8.5"<<std::endl;
                    // adv_fw(it_ele, n_it_ele, it_adj_ele, rbc, dt, rk_iter);
                    // std::cout<<it_ele->affine_pt[1][0]<<std::endl;
                    // std::cout<<it_ele->affine_pt[1][1]<<std::endl;
                    // std::cout<<it_ele->affine_pt[1][2]<<std::endl;
                    // adv_fw_LDG_update_u(it_ele, n_it_ele, it_adj_ele, n_it_adj_ele, 
                    //                     rbc, dt, time, sqrt_diff_ceof, sink_coef, source_coef, rk_iter, INCL_ADV);

                    adv_fw_LDG_update_u(it_ele, n_it_ele, it_adj_ele, n_it_adj_ele, 
                                        rbc, dt, time, sqrt_diff_ceof, k_0, k_1, k_2,
                                        k_3, k_4, k_ss, beta, gamma, q1, h, 
                                        rk_iter, INCL_ADV);
                    
                    //Should we update b, the inhibitor of u, here as well? But since the evolution of b takes the form of an ODE,
                    //how should it be handled? FOR NOW, let's update it in the adv_fw_LDG_update_u function.
                    
                    
                    //std::cout<<"ERROR 8.6"<<std::endl;
                    if((rk_iter + 1) == RK_STEP)
                    {
                        double tmp_max_dt; 
                        /*if(YES == INCL_ADV)
                           tmp_max_dt = time_step_on_tri(nfr, n_it_ele);*/
                        // max_dt = /*Time_step_factor(fr)*/time_step_factor*0.5*
                        //          /*sqr(min_edge_len)*/(min_edge_len*min_edge_len)/(sqrt_D*sqrt_D)/*prmslst[0]->shear_visc*/;
                        max_dt = time_step_factor*0.5*
                                 (min_edge_len*min_edge_len)/(sqrt_D*sqrt_D);
                        /*if(YES == INCL_ADV)
                            max_dt = min(tmp_max_dt, max_dt); */
                        //newdt = min(newdt, max_dt); //20191006
                    }

                    
                } // for(i=0;(it_ele != rbc->Curr_ele->end());i++, it_ele++, n_it_ele++)
//std::cout<<"ERROR 9"<<std::endl;
                /*:: Compute q by using updated u ::*/
                it_ele = rbc->Curr_ele->begin();
                n_it_ele = n_rbc->Curr_ele->begin();
               // std::cout<<"ERROR 10"<<std::endl;
                for(i=0;(it_ele != rbc->Curr_ele->end());i++, it_ele++, n_it_ele++)
                {
                    if (i >= rbc->ENum){
                        break;
                    }
                    if (it_ele->vertex[0] >= (INT_MAX-100) || it_ele->vertex[1] >= (INT_MAX-100) || it_ele->vertex[2] >= (INT_MAX-100)){
                        continue;
                    }
                    if (it_ele->vertex[0] <= (-INT_MAX+100) || it_ele->vertex[1] <= (-INT_MAX+100) || it_ele->vertex[2] <= (-INT_MAX+100)){
                        continue;
                    }
                   
                    for(side = 0; side < 3; side++)
                        it_adj_ele[side] = rbc->sv_it_adj_ele[it_ele->Id][side];
                    for(side = 0; side < 3; side++) 
                        n_it_adj_ele[side] = n_rbc->sv_it_adj_ele[n_it_ele->Id][side];

                    for(side = 0; side < 3; side++) {
                        if(n_it_adj_ele[side]->Id != it_adj_ele[side]->Id) {
                            printf("ERROR: LDG_Surf_diffusion(). 2\n");
                            printf("side %d. New and old ele list not"
                                   " consistent(1-3) %d, %d\n",
                                    side, n_it_adj_ele[side]->Id, it_adj_ele[side]->Id);
                            //clean_up(ERROR); //20191006
                        }
                    }

                    adv_fw_LDG_update_q(it_ele, n_it_ele, it_adj_ele, n_it_adj_ele, 
                                rbc, dt, sqrt_diff_ceof, rk_iter);

                    
                } // 2nd: for(i=0;(it_ele != rbc->Curr_ele->end());i++, it_ele++, n_it_ele++)
              //  std::cout<<"ERROR 11"<<std::endl;
            } // for(i_NumRBC = 0; i_NumRBC < Num_RBCs; i_NumRBC++)

            /*if((rk_iter + 1) == RK_STEP)
            {
                Solid_max_dt(nfr) = newdt; 
            }*/ //20191006, not sure what this is for.
        } // for(rk_iter = 0; rk_iter < RK_STEP; rk_iter++)

        time = fr_time + dt;//fr->time + dt; //20191006
        

        // DEBUG_LEAVE(LDG_Surf_diffusion)
}



// void Utilities::Allocate_memory(RBC *rbc) 
// {
//         int   i, j, k; 

//         /////// allocate memory for the first cell 
// 	/* RBC *rbc=(RBC*)front->RBC_pointer[0]; */ //20191008

//        /* for(i = 0; i < MAXD; i++)
//         {
//             for(j = 0; j < 2; j++)
//             {
//                // rbc->Buf_ele_num[i][j] = 0;
//                // rbc->Buf_ele[i][j] = NULL;
//                // rbc->Buf_ele_index[i][j] = NULL;
//             }
//         }*/

//         FT_vector(&rbc->nodeloc, rbc->NNum, sizeof(CVector));
//         FT_vector(&rbc->ele, rbc->ENum, sizeof(ELEMENT));
//         FT_vector(&rbc->node, rbc->NNum, sizeof(CVector));
//         FT_vector(&rbc->normal_node_fitting, rbc->NNum, sizeof(CVector));
//         FT_vector(&rbc->nw, rbc->NNum, sizeof(double));
// 	FT_vector(&rbc->la[0],rbc->NNum,sizeof(CVector));
//         FT_vector(&rbc->la[1],rbc->NNum,sizeof(CVector));
//         FT_vector(&rbc->node_nb,rbc->NNum*8*2,sizeof(int));

//         // added by S. Xu
//         FT_vector(&rbc->node0, rbc->NNum, sizeof(CVector));
//         FT_vector(&rbc->nodep, rbc->NNum, sizeof(CVector));
//         FT_vector(&rbc->nodeVel, rbc->NNum, sizeof(CVector));
//         FT_vector(&rbc->fn, rbc->NNum, sizeof(CVector));
//         /*FT_vector(&rbc->fund_geom, rbc->NNum, sizeof(FUND_GEOM));*/
//         //END:: added by S. Xu

// 	using namespace std;
// 	//scalar(&rbc->Curr_ele_index,sizeof(list<int>));
//         //scalar(&rbc->Curr_ele,sizeof(list<ELEMENT>));

// 	//rbc->Curr_ele_index = new list<int>;
// 	//rbc->Curr_ele = new list<ELEMENT>;
//        /* rbc->cur_node_index = new list<int>;
//         rbc->cur_node_priority = new list<int>;
//         rbc->curr_link_index = new list<int>;

        
//         rbc->Mean_curv =new double [rbc->NNum];
//         rbc->Gauss_curv =new double [rbc->NNum];

//         rbc->geom_g = new double*[rbc->NNum];
//         for(i=0;i<rbc->NNum;i++)
//             rbc->geom_g[i] = new double[3];

//         rbc->inv_geom_uv = new double*[rbc->NNum];
//         for(i=0;i<rbc->NNum;i++)
//             rbc->inv_geom_uv[i] = new double[4];
//         rbc->der_inv_geom_uv = new double*[rbc->NNum];
//         for(i=0;i<rbc->NNum;i++)
//             rbc->der_inv_geom_uv[i] = new double[4];

//         rbc->lambda_inex = new double[rbc->NNum];

//         // added by S. Xu
//         rbc->fMarker = new double[rbc->NNum];
//         rbc->chem = new double[rbc->NNum];*/
//         //END:: added by S. Xu
//         /**************************/

//         /**
//         rbc->nodeloc = (CVector *) malloc(sizeof(CVector)*rbc->NNum);
//         rbc->ele = (ELEMENT *) malloc(sizeof(ELEMENT)*rbc->ENum);
//         rbc->node = (CVector *) malloc(sizeof(CVector)*rbc->NNum);
//         rbc->nw = (double *) malloc(sizeof(double)*rbc->NNum);
//         **/

//         rbc->sv_it_adj_ele = Allocate_2D_matrix<list<ELEMENT>::iterator>( rbc->ENum, 3);
// }

double Utilities::current_area(RBC *rbc)
{
    std::list<ELEMENT> *Cur_ele = rbc->Curr_ele;
	std::list<ELEMENT>::iterator it=Cur_ele->begin();
        int      j;
        CVector  vt[MAXD];
        double   dtemp, one_third, total_area = 0.0;

        one_third = 1.0/3.0;
        if (Cur_ele->begin() == Cur_ele->end()){
            std::cout<<"Cur_ele has its begin == end...."<<std::endl;
        }
        int COUNTING = 0;
        int k = 0;
        while(it!=Cur_ele->end())
        {
            if (k >= rbc->ENum){
                break;
            }
            if (it->vertex[0] >= (INT_MAX-100) || it->vertex[1] >= (INT_MAX-100) || it->vertex[2] >= (INT_MAX-100)){
                it++;
                continue;
            }
            if (it->vertex[0] <= (-INT_MAX+100) || it->vertex[1] <= (-INT_MAX+100) || it->vertex[2] <= (-INT_MAX+100)){
                it++;
                continue;
            }
            COUNTING += 1;
            //std::cout<<"Cause of seg fault 1"<<std::endl;
            for(j = 0; j < 3; j++){
               // std::cout<<it->vertex[j]<<std::endl;
                vt[j] = rbc->node[it->vertex[j]];
            }
            //std::cout<<"Cause of seg fault 2"<<std::endl;
            it->out_norm=Cross(vt[1]-vt[0],vt[2]-vt[0]);
            //std::cout<<"Cause of seg fault 3"<<std::endl;
            dtemp=Modul(it->out_norm);
            //std::cout<<"Cause of seg fault 4"<<std::endl;
            it->area=0.5*dtemp;
            it->area0 = it->area;
            //std::cout<<"Cause of seg fault 5"<<std::endl;
            //std::cout<<(it->area)<<std::endl;

            //total_area += it->area;
            //std::cout<<"Cause of seg fault 6"<<std::endl;
            it->center[0] = one_third*(vt[0].x[0] + vt[1].x[0] + vt[2].x[0]);
            it->center[1] = one_third*(vt[0].x[1] + vt[1].x[1] + vt[2].x[1]);
            it->center[2] = one_third*(vt[0].x[2] + vt[1].x[2] + vt[2].x[2]);
            //std::cout<<"Cause of seg fault 7"<<std::endl;
            if (it == Cur_ele->begin()){
               // std::cout<<it->area<<std::endl;
               // std::cout<<it->out_norm.x[0]<<std::endl;
               // std::cout<<it->out_norm.x[1]<<std::endl;
               // std::cout<<it->out_norm.x[2]<<std::endl;
               //ALTHOUGH there is something weird about this function,
               //such that, even though the computation does not really normalize
               //the out_norm vector, the actually saved version is the
               //normalized vectors.
            }
            it++;
            k+=1;
        }
        std::cout<<"COUNTING = "<<COUNTING<<std::endl;

#if defined(__MPI__)
        pp_global_sum(&total_area, 1); 
#endif /* #if defined(__MPI__) */
        return total_area;
};

void Utilities::build_ele_adj_by_iter(RBC* rbc, CoordInfoVecs& coordInfoVecs)
{
	int    i, side;
	std::list<ELEMENT>::iterator  it_ele; 
    std::list<ELEMENT>::iterator it_ele2;
        std::list<ELEMENT>::iterator it_adj_ele[3]; 
         
        // rbc->sv_it_adj_ele = Allocate_2D_matrix<list<ELEMENT>::iterator>(
        //                          (rbc)->ENum, 3);

        // std::cout<<"rbc->sv_it_adj_ele = " << rbc->sv_it_adj_ele<<std::endl; 

        it_ele = rbc->Curr_ele->begin();
        
        for(i=0;(it_ele != rbc->Curr_ele->end());i++, it_ele++) {
            if (i >= rbc->ENum){
                break;
            }
            if (it_ele->vertex[0] >= (INT_MAX-100) || it_ele->vertex[1] >= (INT_MAX-100) || it_ele->vertex[2] >= (INT_MAX-100)){
                continue;
            }
            if (it_ele->vertex[0] <= (-INT_MAX+100) || it_ele->vertex[1] <= (-INT_MAX+100) || it_ele->vertex[2] <= (-INT_MAX+100)){
                continue;
            }
            for(side = 0; side < 3; side++) {
                it_ele2 = rbc->Curr_ele->begin();
                it_adj_ele[side] = rbc->Curr_ele->begin();
                //std::cout<<it_adj_ele[side]->Id<<std::endl;
                // if (side == 0){
                //     int max_search = coordInfoVecs.triangles2Triangles_1[i];
                //     //std::cout<<"max_search = "<<max_search<<std::endl;
                //     for (int j = 0; it_ele2 != rbc->Curr_ele->end();j++, it_ele2++){
                //         it_adj_ele[side] = it_ele2;
                //         if (j == max_search){
                //            // std::cout<<"j = "<<j<<std::endl;
                //             break;
                //         }
                //     }
                // }
                // if (side == 1){
                //     int max_search = coordInfoVecs.triangles2Triangles_2[i];
                //     //std::cout<<"max_search = "<<max_search<<std::endl;
                //     for (int j = 0; it_ele2 != rbc->Curr_ele->end();j++, it_ele2++){
                //         it_adj_ele[side] = it_ele2;
                //         if (j == max_search){
                //       //      std::cout<<"j = "<<j<<std::endl;
                //             break;
                //         }
                //     }
                // }
                // if (side == 2){
                //     int max_search = coordInfoVecs.triangles2Triangles_3[i];
                //     //std::cout<<"max_search = "<<max_search<<std::endl;
                //     for (int j = 0; it_ele2 != rbc->Curr_ele->end();j++, it_ele2++){
                //         it_adj_ele[side] = it_ele2;
                //         if (j == max_search){
                //         //    std::cout<<"j = "<<j<<std::endl;
                //             break;
                //         }
                //     }
                // }
                it_adj_ele[side] = find_if(it_adj_ele[side], rbc->Curr_ele->end(),
                             std::bind2nd( loc_EleNumber(), it_ele->adj_ele[side]) );    
                if(it_adj_ele[side] == rbc->Curr_ele->end())
                {
                    printf("ERROR: build_ele_adj_by_iter()."
                             " adj ele not found in list\n");
                    //clean_up(ERROR);
                }

                rbc->sv_it_adj_ele[it_ele->Id][side] = it_adj_ele[side];
                //std::cout<<"sv_it_adj_ele["<<it_ele->Id<<"]["<<side<<"] = "<<rbc->sv_it_adj_ele[it_ele->Id][side]->Id<<std::endl;
                // std::cout<<"it_adj_ele["<<side<<"] = "<<it_adj_ele[side]->Id<<std::endl;
            }
        }
};

void Utilities::comput_aff_conormals(RBC* rbc)
{
	std::list<ELEMENT>::iterator it_ele;
        std::list<ELEMENT>::iterator it_adj_ele[3];
        int       i, side, nbside; 
        double    **affine_trans; 
	double    conorm[4], adj_ele_conorm[4];
        double    aff_conorm[4], aff_adj_ele_conorm[4]; 

	it_ele = rbc->Curr_ele->begin();
	for(int i = 0;(it_ele != rbc->Curr_ele->end()); i++, it_ele++){
        if (i >= rbc->ENum){
            break;
        }
        if (it_ele->vertex[0] >= (INT_MAX-100) || it_ele->vertex[1] >= (INT_MAX-100) || it_ele->vertex[2] >= (INT_MAX-100)){
            continue;
        }
        if (it_ele->vertex[0] <= (-INT_MAX+100) || it_ele->vertex[1] <= (-INT_MAX+100) || it_ele->vertex[2] <= (-INT_MAX+100)){
            continue;
        }
        //std::cout<<it_ele->Id<<std::endl;
	    for(side = 0; side < 3; side++){
                it_adj_ele[side] = rbc->sv_it_adj_ele[it_ele->Id][side];
                //std::cout<<it_adj_ele[side]->Id<<std::endl;
        }
            affine_trans = it_ele->affine_trans;
            for(side = 0; side < 3; side++)
            {
                for(nbside = 0; nbside < 3; nbside++)
                {
                    if(it_adj_ele[side]->adj_ele[nbside] == it_ele->Id)
                        break;
                }
                
                if(nbside >= 3 || it_ele->Id == it_adj_ele[side]->Id)
                {
                    printf("ERROR: comput_aff_conormals(),"
                        " adj_ele, ele do not share common edge\n");
                    printf("Or Ids of eles = %d, %d\n", it_ele->Id, it_adj_ele[side]->Id);
                    //clean_up(ERROR);
                }
                //std::cout<<"it_ele->side_norm[side].x"<<it_ele->side_norm[side].x[0]<<" "<<it_ele->side_norm[side].x[1]<<" "<<it_ele->side_norm[side].x[2]<<std::endl;
                memcpy(conorm, it_ele->side_norm[side].x, sizeof(double)*3 );
                conorm[3] = 0.0;
                // std::cout<<"conorm = "<<conorm[0]<<" "<<conorm[1]<<" "<<conorm[2]<<std::endl;
                
                matrix_vec_mult(affine_trans, conorm, 4, 4, aff_conorm);
                
                memcpy(it_ele->aff_conorm[side], aff_conorm, sizeof(double)*4);
                if ( i == 0){
                    // std::cout<<"affine_trans = "<<affine_trans[0][0]<<" "<<affine_trans[0][1]<<" "<<affine_trans[0][2]<<" "<<affine_trans[0][3]<<std::endl;
                    // std::cout<<"affine_trans = "<<affine_trans[1][0]<<" "<<affine_trans[1][1]<<" "<<affine_trans[1][2]<<" "<<affine_trans[1][3]<<std::endl;
                    // std::cout<<"affine_trans = "<<affine_trans[2][0]<<" "<<affine_trans[2][1]<<" "<<affine_trans[2][2]<<" "<<affine_trans[2][3]<<std::endl;
                    // std::cout<<"affine_trans = "<<affine_trans[3][0]<<" "<<affine_trans[3][1]<<" "<<affine_trans[3][2]<<" "<<affine_trans[3][3]<<std::endl;
                    // std::cout<<"it_ele->aff_conorm[] ="<<it_ele->aff_conorm[0][0]<<" "<<it_ele->aff_conorm[0][1]<<" "<<it_ele->aff_conorm[0][2]<<" "<<it_ele->aff_conorm[0][3]<<std::endl;
                    // std::cout<<"it_ele->aff_conorm[] ="<<it_ele->aff_conorm[1][0]<<" "<<it_ele->aff_conorm[1][1]<<" "<<it_ele->aff_conorm[1][2]<<" "<<it_ele->aff_conorm[1][3]<<std::endl;
                    // std::cout<<"it_ele->aff_conorm[] ="<<it_ele->aff_conorm[2][0]<<" "<<it_ele->aff_conorm[2][1]<<" "<<it_ele->aff_conorm[2][2]<<" "<<it_ele->aff_conorm[2][3]<<std::endl;
                    // std::cout<<"it_ele->aff_conorm[<<"side"<<] ="<<rbc->ele[i].aff_conorm[0][0]<<" "<<rbc->ele[i].aff_conorm[0][1]<<" "<<rbc->ele[i].aff_conorm[0][2]<<" "<<rbc->ele[i].aff_conorm[0][3]<<std::endl;
                    // std::cout<<"it_ele->aff_conorm[<<"side"<<] ="<<rbc->ele[i].aff_conorm[1][0]<<" "<<rbc->ele[i].aff_conorm[1][1]<<" "<<rbc->ele[i].aff_conorm[1][2]<<" "<<rbc->ele[i].aff_conorm[1][3]<<std::endl;
                    // std::cout<<"it_ele->aff_conorm[<<"side"<<] ="<<rbc->ele[i].aff_conorm[2][0]<<" "<<rbc->ele[i].aff_conorm[2][1]<<" "<<rbc->ele[i].aff_conorm[2][2]<<" "<<rbc->ele[i].aff_conorm[2][3]<<std::endl;
                    std::cout<<" "<<std::endl;
                }
                memcpy(adj_ele_conorm, it_adj_ele[side]->side_norm[nbside].x, 
                  sizeof(double)*3);
                  
                adj_ele_conorm[3] = 0.0;
                
                matrix_vec_mult(affine_trans, adj_ele_conorm, 4, 4, aff_adj_ele_conorm);
                
                memcpy(it_ele->aff_adj_ele_conorm[side], aff_adj_ele_conorm, sizeof(double)*4);
                
            }
            
        }
        
};
void Utilities::comput_aff_cent_vtx(RBC* rbc)
{
	std::list<ELEMENT>::iterator it_ele;
        int       i, j, side; 
        double    **affine_trans; 
        CVector   *vt[3];
        double    B[4]; 

	it_ele = rbc->Curr_ele->begin();
	for(int k = 0 ;(it_ele != rbc->Curr_ele->end()); it_ele++, k++)
        {
            if (k >= rbc->ENum){
                break;
            }
            if (it_ele->vertex[0] >= (INT_MAX-100) || it_ele->vertex[1] >= (INT_MAX-100) || it_ele->vertex[2] >= (INT_MAX-100)){
                continue;
            }
            if (it_ele->vertex[0] <= (-INT_MAX+100) || it_ele->vertex[1] <= (-INT_MAX+100) || it_ele->vertex[2] <= (-INT_MAX+100)){
                continue;
            }
            for(j=0;j<3;j++)
                vt[j] = &(rbc->node[it_ele->vertex[j]]);
            for(j=0;j<3;j++)
            {
                memcpy(B, vt[j]->x, sizeof(double)*3);
                B[3] = 1.0;
                matrix_vec_mult(it_ele->affine_trans, B, 4, 4, it_ele->affine_pt[j]);
            }
            memcpy(B, it_ele->center, sizeof(double)*3);
            B[3] = 1.0;
            matrix_vec_mult(it_ele->affine_trans, B, 4, 4, it_ele->aff_cent);
        }
};

void Utilities::accurate_sphere_diffusion_initializer(
        RBC                      *rbc,
        std::list<ELEMENT>::iterator  &it_ele, 
        double                   *dg_rho_soln,
        double                   *dg_q[])
{
        double       crds[16][2], v_conU[3][8], glb_crds[16][4], coords[3];
        static double w16[16] ={0.144315607677787,0.095091634267285,0.095091634267285,
                             0.095091634267285,
                             0.103217370534718, 0.103217370534718,0.103217370534718,
                             0.032458497623198,0.032458497623198,0.032458497623198,
                             0.027230314174435,0.027230314174435,0.027230314174435,
                             0.027230314174435,0.027230314174435,0.027230314174435};
        static double w3[3] = {0.33333333333333, 0.33333333333333, 0.33333333333333};
        int           i, j, k, l, indx;
        long double   q[9], qw[9], Bn[9];
        double        rhs[10][64], mulrhs[10][64], dens[10], tmp, q_rhs[2][10][64];
        double        **mass_inv, **mass_matrix, B[4], affine_pt[3][4], *affine_pt_p[3], 
                      aff_cent[4], sqrt_area;
        double        **vec_vh_mass_inv, **vec_vh_mass_matrix;

        static double **inv_affine_tran = NULL;
        CVector       *vt[3];
        static double    **coefV = NULL, **coefA;
        static double    **tmp_outR, **tmp_outQ; 
        int           N_nodes = 16; // same as # of quadrature pts.
        double        coefW[64], gamma[64], coefS[64], rhsb[64], len, epsW, surf[64]; 
        int           debug = NO; 

        //Gas_param     **prmslst;
        int           num_params;
        static double sqrt_diff_ceof = 0.0; // this is used as diffusion coefficient in
                                  // current implementation.

        /*
        if(38 == it_ele->Id)
        {
            debug = YES;
            printf("\n\n\n ******rbc = %p, Ele[%d] entered"
             " accurate_sphere_diffusion_initializer()\n", rbc, it_ele->Id); 
        }
        */

        if(NULL == inv_affine_tran)
            FT_matrix(&inv_affine_tran, 4, 4,sizeof(double));

        if(NULL == coefV)
        {
            coefV = new double*[N_nodes];
            for(i = 0; i < N_nodes; i++)
                coefV[i] = new double[15];

            coefA = new double*[N_nodes];
            for(i = 0; i < N_nodes; i++)
                coefA[i] = new double[15];

            tmp_outR = new double*[N_nodes];
            for(i = 0; i < N_nodes; i++)
                tmp_outR[i] = new double[15];

            tmp_outQ = new double*[N_nodes];
            for(i = 0; i < N_nodes; i++)
                tmp_outQ[i] = new double[N_nodes];

            //num_params = return_params_list(&prmslst);
            // for(i = 0; i < num_params; i++)
            //     fprint_Gas_param(stdout,prmslst[i]);
            sqrt_diff_ceof = 0.0;//sqrt(prmslst[0]->shear_visc);//06252020: set sqrt_diff_coef to be 0 at the moment as I have no clue where to locate values for prmslst.
        }

        mass_matrix = it_ele->mass_matrix;
        mass_inv = it_ele->mass_inv;

        inverse_matrix(it_ele->affine_trans, 4, inv_affine_tran);
        for(j=0;j<3;j++)
            vt[j] = &(rbc->node[it_ele->vertex[j]]);

        for(j=0;j<3;j++)
        {
            B[3] = 1.0;
            memcpy(B, vt[j]->x, sizeof(double)*3);
            matrix_vec_mult(it_ele->affine_trans, B, 4, 4, affine_pt[j]);
            affine_pt_p[j] = affine_pt[j];
        }

        for(j=0;j<3;j++) B[j] = it_ele->center[j];
        B[3] = 1.0;
        matrix_vec_mult(it_ele->affine_trans, B, 4, 4, aff_cent);

        if(YES == debug)
        {
            // printf("Orignal tri points and center\n");
            // print_vector("origna vt0", 3, vt[0]->x, "%g ");
            // print_vector("origna vt1", 3, vt[1]->x, "%g ");
            // print_vector("origna vt2", 3, vt[2]->x, "%g ");
            // print_vector("origna cent", 3, it_ele->center, "%g ");
            // printf("Side Length = %g, %g, %g\n",
            //     it_ele->length_side[0], it_ele->length_side[1], it_ele->length_side[2]); 
            // printf("Affined tri points and center\n");
            // print_vector("affine vt0", 3, affine_pt_p[0], "%g ");
            // print_vector("affine vt1", 3, affine_pt_p[1], "%g ");
            // print_vector("affine vt2", 3, affine_pt_p[2], "%g ");
            // print_vector("affine cent", 3, aff_cent, "%g ");
        }

        sqrt_area = sqrt(it_ele->area0);
        tri_quadrature_16_pts(affine_pt_p[0], affine_pt_p[1], affine_pt_p[2], crds);

        /**
        for(i = 0; i < 3; i++)
        {
            for(j=0;j<2;j++) crds[i][j] = affine_pt_p[i][j];
        }
        **/        
 
        for(i = 0; i < 16; i++)
        // for(i = 0; i < 3; i++)
        {
            B[3] = 1.0; B[2] = 0.0; 
            memcpy(B, crds[i], sizeof(double)*2);
            matrix_vec_mult(inv_affine_tran, B, 4, 4, glb_crds[i]);
        }

        // if(YES == debug)
        // {
        //     for(i = 0; i < 16; i++)
        //         printf("Global quadr[%d] = %g, %g, %g, %g\n", i, 
        //           glb_crds[i][0], glb_crds[i][1], glb_crds[i][2], glb_crds[i][3]);
        // }

        for(i = 0; i < N_EQN; i++)
        {
            for(j = 0; j < MAX_N_COEF; j++)
                rhs[i][j] = 0.0;
        }

        for(j = 0; j < N_EQN; j++)
        {
            for(k = 0; k < MAX_N_COEF; k++)
            {
                rhs[j][k] = 0.0;
                for(i = 0; i < 16; i++)
                // for(i = 0; i < 3; i++)
                {
                    rhs[j][k] += w16[i]*accurate_diffusion_sphere(rbc->Center.x, 
                                 glb_crds[i],rbc->avg_radi, it_ele, 0.0)
                                 *vh_val_loc_div_free_basis(crds[i],aff_cent,sqrt_area,k);
                    /*
                    rhs[j][k] += w3[i]*accurate_diffusion_sphere(rbc->Center.x, 
                                   glb_crds[i],rbc->avg_radi,it_ele, 0.0)
                                 *vh_val_loc_div_free_basis(crds[i],aff_cent,sqrt_area,k);
                    */
                }
            }
        }

        for(j = 0; j < N_EQN; j++)
        {
            for(k = 0; k < MAX_N_COEF; k++)
            {
                rhs[j][k] *= it_ele->area0;
            }
        }

        for(i = 0; i < N_EQN; i++)
            matrix_vec_mult(mass_inv, rhs[i], MAX_N_COEF, MAX_N_COEF, mulrhs[i]);

        memcpy(dg_rho_soln, mulrhs[0], sizeof(double)*MAX_N_COEF); 
        // for(indx = 0; indx < MAX_N_COEF; indx++)
        //     dg_rho_soln[indx] = mulrhs[0][indx];

        // std::cout<<"dg_rho_soln = ["<<dg_rho_soln[0]<<" "<<dg_rho_soln[1]<<" "<<dg_rho_soln[2]<<std::endl;

        // The code set initial dt = 0, which let q be computed.
        return; 
};

// void Utilities::Init_cell_intfc_vel(
// 	RBC* rbc
//     //Wave     *wv,
//       //  Front    *fr
//       )
// {
// 	int           i, j, k, i_NumRBC, Num_RBCs = 1;//fr->Num_RBCs;
// 	//RBC           *rbc = (RBC *)fr->RBC_pointer[0];
//        // POINTER       *RBC_pointers = fr->RBC_pointer;
//         double        coords[3], r, phi, theta; 
//         CVector       *vt[3]; 
//         double        ref_tri_vt[3][4] = {{0.0, 0.0, 0.0, 1.0}, 
//                                           {1.0, 0.0, 0.0, 1.0}, {0.0, 1.0, 0.0, 1.0}}; 
//         double        **A, B[4], soln[3][4], **affine_tran, **inv_affine_tran, **jacobi_x_over_x_hat; 
//         double        jac0[3], jac1[3], jac2[3], det, det2; 
//         double        *tri_norm, tmp_vec[4], tang_x[4], tang_y[4], len, tan_tmp[4]; 
//         double        **rot_matr_rows, **translate;
//         double        Eu_crds[3][3] = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
//         double        dg_rho_soln[10], dg_q[2][10], *dg_q_p[2]; 

//         std::list<ELEMENT>::iterator ele_it;

//         FT_matrix(&A, 3, 3,sizeof(double));
//         FT_matrix(&jacobi_x_over_x_hat, 3, 3,sizeof(double));
//         FT_matrix(&affine_tran, 4, 4,sizeof(double));
//         FT_matrix(&inv_affine_tran, 4, 4,sizeof(double));
//         FT_matrix(&rot_matr_rows, 4, 4,sizeof(double));
//         FT_matrix(&translate, 4, 4,sizeof(double));

//         for(j = 0; j < 4; j++)
//             for(k = 0; k < 4; k++)
//                 translate[j][k] = 0.0;

// 	for(i_NumRBC = 0; i_NumRBC < Num_RBCs; i_NumRBC++)
//         {
//             //rbc = (RBC*)RBC_pointers[i_NumRBC];

//             if (1 > 0)//if(debugging("Surf_diff"))
//             {
//                 for(i = 0; i < rbc->NNum; i++)
//                 {
//                     coords[0] = rbc->node[i].x[0] - rbc->Center.x[0]; 
//                     coords[1] = rbc->node[i].x[1] - rbc->Center.x[1]; 
//                     coords[2] = rbc->node[i].x[2] - rbc->Center.x[2]; 
//                     rbc->nodeVel[i].x[0] = 0.0;
//                     rbc->nodeVel[i].x[1] = 0.0;
//                     rbc->nodeVel[i].x[2] = 0.0; 

//                     r = sqrt(sqr(coords[0]) + sqr(coords[1]) + sqr(coords[2])); 
//                     // theta = atan(coords[1]/coords[0]);
//                     theta = atan2(coords[1], coords[0]);
         
//                     // if(fabs(coords[0])< 1.0e-12)
//                     /**
//                     if(i == 35)
//                     {
//                         printf("ERROR: Init_cell_intfc_vel. cords[0] = %g, crds[1] = %g, theta %g r = %g\n", 
//                                coords[0], coords[1], theta, r);
//                         printf("coords = %13.12g, %13.12g, %13.12g\n", coords[0], coords[1], coords[2]); 
//                         rbc->node[i].Printf("node crds");
//                         rbc->Center.Printf("icent crds");
//                         // clean_up(ERROR); 
//                     }
//                     **/

//                     phi = acos(coords[2]/r); 
//                     rbc->chem[i] = sqr(sin(theta))*sqr(sin(theta))*sin(theta)*cos(5.0*phi) + 
//                                    sqr(sin(theta))*sqr(sin(theta))*cos(theta)*cos(4.0*phi);
//                 } // for(i = 0; i < rbc->NNum; i++)

//                 dg_q_p[0] = dg_q[0]; 
//                 dg_q_p[1] = dg_q[1]; 
//                 for(i = 0; i < 2; i++)
//                 {
//                     for(k = 0; k < MAX_N_COEF; k++)
//                         dg_q[i][k] = 0.0;
//                 }

//                 for(ele_it = rbc->Curr_ele->begin();(ele_it!=rbc->Curr_ele->end());ele_it++)
//                 {
//                     accurate_sphere_diffusion_initializer(rbc, ele_it, dg_rho_soln, dg_q_p); 
//                     // for(i = 0; i < MAX_N_COEF; i++)
//                     //     ele_it->dg_rho[0][i] = dg_rho_soln[i]; 

//                     memcpy(ele_it->dg_rho[0], dg_rho_soln, sizeof(double)*MAX_N_COEF);
//                     memcpy(ele_it->dg_q[0][0], dg_q_p[0], sizeof(double)*MAX_N_COEF); 
//                     memcpy(ele_it->dg_q[0][1], dg_q_p[1], sizeof(double)*MAX_N_COEF); 

//                     if(0 == ele_it->Id)
//                     {
//                         //printf("\n\n\nIn Init_cell_intfc_vel(). Ele[%d] init condiiton\n", 
//                                ele_it->Id); 
//                         //print_vector("u  soln", MAX_N_COEF, ele_it->dg_rho[0], "%g "); 
//                         //print_vector("q0 soln", MAX_N_COEF, ele_it->dg_q[0][0], "%g "); 
//                         //print_vector("q1 soln", MAX_N_COEF, ele_it->dg_q[0][1], "%g "); 
//                     }
//                 }
//             }//END:: if(debugging("Surf_diff"))

//         }//END:: for(i_NumRBC = 0; i_NumRBC < Num_RBCs; i_NumRBC++)

//         // free(A); free(affine_tran); free(inv_affine_tran); free(jacobi_x_over_x_hat); 
//         // free(rot_matr_rows); free(translate); 
// };


// double Utilities::derivative_P4_surf(
// 	double   *surf,
// 	double   *x,
//         int      indx)
// {
//         double ans = 0.0; 
//         if(0 == indx) //x derivative
//         {
//             ans += surf[1];
//             ans += surf[3] * x[0]; 
//             ans += surf[4] * x[1]; 
//             ans += surf[6] * 3.0*sqr(x[0])/6.0; 
//             ans += surf[7] * 2.0*x[0]*x[1]*0.5;
//             ans += surf[8] * sqr(x[1])*0.5; 
//             ans += surf[10] * 4.0*cub(x[0])/24.0;
//             ans += surf[11] * 3.0*sqr(x[0])*x[1]/6.0;
//             ans += surf[12] * 2.0*x[0]*sqr(x[1])*0.25;
//             ans += surf[13] * sqr(x[1])*x[1]/6.0;
//         }
//         else
//         {
//             ans += surf[2];
//             ans += surf[4] * x[0]; 
//             ans += surf[5] * x[1];
//             ans += surf[7] * sqr(x[0])*0.5;
//             ans += surf[8] * x[0]*2.0*x[1]*0.5; 
//             ans += surf[9] * 3.0*sqr(x[1])/6.0;
//             ans += surf[11] * sqr(x[0])*x[0]/6.0;
//             ans += surf[12] * sqr(x[0])*2.0*x[1]*0.25;
//             ans += surf[13] * x[0]*3.0*sqr(x[1])/6.0;
//             ans += surf[14] * 4.0*cub(x[1])/24.0;
//         }

//         return ans; 
// };

double Utilities::accurate_diffusion_sphere(
	double                   *glb_cell_cent, 
	double                   *glb_crds,
        double                   avg_radi,
        std::list<ELEMENT>::iterator  &it_ele,
	double                   time)
{
        double        coords[3], B[4]; 
        double        phi, r, theta, chem, tmp; 
        double        **affine_trans = it_ele->affine_trans; 
        static double **inv_affine_tran = NULL;
        int           debug = NO;
        
        /**
        if(1 == it_ele->Id)
        {
            debug = YES;
            printf("Ele[%d] enter accurate_diffusion_sphere at time %g\n", 
                     it_ele->Id, time); 
            printf("Glb center (%g, %g, %g); pt (%g, %g, %g) avg_radi = %g\n",
                 glb_cell_cent[0], glb_cell_cent[1], glb_cell_cent[2], 
                 glb_crds[0], glb_crds[1], glb_crds[2], avg_radi); 
        }
        **/

        if(NULL == inv_affine_tran)
            FT_matrix(&inv_affine_tran, 4, 4,sizeof(double));

        coords[0] = glb_crds[0] - glb_cell_cent[0];
        coords[1] = glb_crds[1] - glb_cell_cent[1];
        coords[2] = glb_crds[2] - glb_cell_cent[2];

        r = sqrt(sqr(coords[0]) + sqr(coords[1]) + sqr(coords[2]));
        /** // for cos(phi) case.
        theta = atan2(coords[1], coords[0]);
        phi = acos(coords[2]/r);
        **/

        phi = atan2(coords[1], coords[0]);
        theta = acos(coords[2]/r);

        chem = sqr(sin(theta))*sqr(sin(theta))*sin(theta)*cos(5.0*phi) +
               sqr(sin(theta))*sqr(sin(theta))*cos(theta)*cos(4.0*phi);
        tmp = exp(-30.0*time/(sqr(r))); 

        /**
        chem = cos(phi); 
        tmp = exp(-2.0*time/(sqr(avg_radi))); 
        **/

        chem *= tmp;
        // return 10.0; 

        // if(YES == debug)
        // {
        //     printf("computed radi = %14.12g, theta = %g, phi = %g\n", r, theta, phi); 
        // }

        return chem; 
};
void Utilities::LDG_Error_of_surf_diffusion_test(
	RBC      *rbc, 
	double   time, 
	double   *diff_error)
{
	std::list<ELEMENT>::iterator it_ele;
        int       i, j, indx; 
        long      ele_max_L8;
        double    coords[3], crds[16][2], glb_crds[16][4]; 
        static double w16[16] ={0.144315607677787,0.095091634267285,0.095091634267285,
                             0.095091634267285,
                             0.103217370534718, 0.103217370534718,0.103217370534718,
                             0.032458497623198,0.032458497623198,0.032458497623198,
                             0.027230314174435,0.027230314174435,0.027230314174435,
                             0.027230314174435,0.027230314174435,0.027230314174435};
        static double w3[3] = {0.33333333333333, 0.33333333333333, 0.33333333333333};

        CVector       *vt[3];
        double        B[4]; 
        double        **mass_inv, **mass_matrix, affine_pt[3][4], *affine_pt_p[3],
                      aff_cent[4], sqrt_area;
	static double **inv_affine_tran = NULL;
        double        num_soln, acc_soln, error_L1 = 0, error_Linf = 0, temp; 
        double        error_L2 = 0.0, temp_sqr; 
        double        total_u = 0.0, cell_total_u, total_area = 0.0, total_area_diff =0.0; 

        if(NULL == inv_affine_tran)
            FT_matrix(&inv_affine_tran, 4, 4,sizeof(double));

	it_ele = rbc->Curr_ele->begin();
	for(int k = 0;(it_ele != rbc->Curr_ele->end());it_ele++, k++)
        {
            if (k >= rbc->ENum){
                break;
            }
            mass_matrix = it_ele->mass_matrix;
            mass_inv = it_ele->mass_inv;

            cell_total_u = 0.0;
            for(indx = 0; indx < MAX_N_COEF; indx++)
                cell_total_u += it_ele->dg_rho[0][indx]*mass_matrix[0][indx]; 
            total_u += cell_total_u;
            total_area += mass_matrix[0][0];
            total_area_diff += fabs(mass_matrix[0][0]-it_ele->area0); 

            inverse_matrix(it_ele->affine_trans, 4, inv_affine_tran);

            for(j=0;j<3;j++)
                vt[j] = &(rbc->node[it_ele->vertex[j]]);

            for(j=0;j<3;j++)
            {
                B[3] = 1.0;
                memcpy(B, vt[j]->x, sizeof(double)*3);
                matrix_vec_mult(it_ele->affine_trans, B, 4, 4, affine_pt[j]);
                affine_pt_p[j] = affine_pt[j];
            }

            for(j=0;j<3;j++) B[j] = it_ele->center[j];
            B[3] = 1.0;
            matrix_vec_mult(it_ele->affine_trans, B, 4, 4, aff_cent);

            sqrt_area = sqrt(it_ele->area0);
            tri_quadrature_16_pts(affine_pt_p[0], affine_pt_p[1], affine_pt_p[2], crds);

            for(i = 0; i < 3; i++)
            {
                for(j=0;j<2;j++)
                    crds[i][j] = affine_pt_p[i][j];
            }

            // for(i = 0; i < 16; i++)
            for(i = 0; i < 3; i++)
            {
                B[3] = 1.0; B[2] = 0.0;
                memcpy(B, crds[i], sizeof(double)*2);
                matrix_vec_mult(inv_affine_tran, B, 4, 4, glb_crds[i]);
            }

            temp = 0.0; temp_sqr = 0.0; 
            // for(i = 0; i < 16; i++)
            for(i = 0; i < 3; i++)
            {
                soln_at_pt(it_ele->dg_rho[0], crds[i], aff_cent, sqrt_area, &num_soln);
                acc_soln = accurate_diffusion_sphere(rbc->Center.x, glb_crds[i], 
                                   rbc->avg_radi, it_ele, time); 
                                   std::cout<<"@ elem "<<it_ele->Id<<", acc_soln = "<<acc_soln<<" "<<"num_soln = "<<num_soln<<std::endl;
                // temp += fabs(acc_soln-num_soln)*w16[i];
                temp += fabs(acc_soln-num_soln)*w3[i];
                temp_sqr += sqr(acc_soln-num_soln)*w3[i];
                if(fabs(acc_soln-num_soln) > error_Linf)
                {
                    error_Linf = fabs(acc_soln-num_soln); 
                    ele_max_L8 = it_ele->Id; 
                }
            }
            error_L1 += temp*it_ele->area0;
            error_L2 += temp_sqr*it_ele->area0;
        }
        // std::cout<<"total_u = "<<total_u<<std::endl;
        std::cout<<"error_L1 = "<<error_L1<<std::endl;
        std::cout<<"error_L2 = "<<error_L2<<std::endl;
        // pp_global_sum(&error_L1, 1);
        // pp_global_max(&error_Linf, 1);
        // pp_global_max(&error_L2, 1);
        // pp_global_lmax(&ele_max_L8, 1);
        // pp_global_sum(&total_u, 1);
        // pp_global_sum(&total_area, 1);
        // pp_global_sum(&total_area_diff, 1);

        // printf("t(%g), L1_e %14.12g; L2_e %14.12g; L8_e %14.12g on ele[%d], total_u %14.12g, t_area %14.12g\n",
        //          time, error_L1, error_L2, error_Linf, ele_max_L8, total_u, total_area);

};

void Utilities::LDG_Error_of_surf_diffusion_test_on_ele(
	RBC      *rbc, 
	std::list<ELEMENT>::iterator it_ele,
	double   time)
{
        int       i, j, indx; 
        long      ele_max_L8;
        double    coords[3], crds[16][2], glb_crds[16][4]; 
        static double w16[16] ={0.144315607677787,0.095091634267285,0.095091634267285,
                             0.095091634267285,
                             0.103217370534718, 0.103217370534718,0.103217370534718,
                             0.032458497623198,0.032458497623198,0.032458497623198,
                             0.027230314174435,0.027230314174435,0.027230314174435,
                             0.027230314174435,0.027230314174435,0.027230314174435};
        static double w3[3] = {0.33333333333333, 0.33333333333333, 0.33333333333333};

        CVector       *vt[3];
        double        B[4]; 
        double        **mass_inv, **mass_matrix, affine_pt[3][4], *affine_pt_p[3],
                      aff_cent[4], sqrt_area;
	static double **inv_affine_tran = NULL;
        double        num_soln, acc_soln, error_L1 = 0, error_Linf = 0, temp; 
        double        total_u = 0.0, cell_total_u, total_area = 0.0, total_area_diff =0.0; 

        if(NULL == inv_affine_tran)
            FT_matrix(&inv_affine_tran, 4, 4,sizeof(double));

        mass_matrix = it_ele->mass_matrix;
        mass_inv = it_ele->mass_inv;

        cell_total_u = 0.0;
        for(indx = 0; indx < MAX_N_COEF; indx++)
            cell_total_u += it_ele->dg_rho[0][indx]*mass_matrix[0][indx]; 
        total_u += cell_total_u;
        total_area += mass_matrix[0][0];
        total_area_diff += fabs(mass_matrix[0][0]-it_ele->area0); 

        inverse_matrix(it_ele->affine_trans, 4, inv_affine_tran);

        for(j=0;j<3;j++)
            vt[j] = &(rbc->node[it_ele->vertex[j]]);

        for(j=0;j<3;j++)
        {
            B[3] = 1.0;
            memcpy(B, vt[j]->x, sizeof(double)*3);
            matrix_vec_mult(it_ele->affine_trans, B, 4, 4, affine_pt[j]);
            affine_pt_p[j] = affine_pt[j];
        }

        for(j=0;j<3;j++) B[j] = it_ele->center[j];
        B[3] = 1.0;
        matrix_vec_mult(it_ele->affine_trans, B, 4, 4, aff_cent);

        sqrt_area = sqrt(it_ele->area0);
        tri_quadrature_16_pts(affine_pt_p[0], affine_pt_p[1], affine_pt_p[2], crds);

        for(i = 0; i < 3; i++)
        {
            for(j=0;j<2;j++)
                crds[i][j] = affine_pt_p[i][j];
        }

        // for(i = 0; i < 16; i++)
        for(i = 0; i < 3; i++)
        {
            B[3] = 1.0; B[2] = 0.0;
            memcpy(B, crds[i], sizeof(double)*2);
            matrix_vec_mult(inv_affine_tran, B, 4, 4, glb_crds[i]);
        }

        temp = 0.0; 
        // for(i = 0; i < 16; i++)
        for(i = 0; i < 3; i++)
        {
            soln_at_pt(it_ele->dg_rho[0], crds[i], aff_cent, sqrt_area, &num_soln);
            acc_soln = accurate_diffusion_sphere(rbc->Center.x, glb_crds[i], 
                               rbc->avg_radi, it_ele, time); 
                               std::cout<<"num_soln = "<<num_soln<<" , acc_soln = "<<acc_soln<<std::endl;
            // printf("At pt[%g, %g], acc_soln = %12.10g, num_soln = %12.10g\n",
            //            crds[i][0], crds[i][1], acc_soln, num_soln);
        }

};
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Utilities::LDG_Surface_Diffusion_Initialize(
    CoordInfoVecs& coordInfoVecs,
    GeneralParams& generalParams,
    HostSetInfoVecs& hostSetInfoVecs,
    AuxVecs& auxVecs,
	RBC      *rbc,
    RBC      *n_rbc,
    uint mem_prealloc
    ){
    transferDtoH(generalParams, coordInfoVecs, hostSetInfoVecs);
    for (int i = 0; i < coordInfoVecs.num_triangles; i++){
        triangles2Triangles_host_vecs(i, hostSetInfoVecs,coordInfoVecs,generalParams, auxVecs);
    }
    transferHtoD(generalParams,coordInfoVecs,hostSetInfoVecs);

    rbc->NNum = generalParams.maxNodeCount;
    rbc->ENum = coordInfoVecs.num_triangles;    

    FT_vector(&rbc->nodeloc, mem_prealloc*rbc->NNum, sizeof(CVector));
    FT_vector(&rbc->ele, mem_prealloc*rbc->ENum, sizeof(ELEMENT));
    FT_vector(&rbc->node, mem_prealloc*rbc->NNum, sizeof(CVector));
    FT_vector(&rbc->normal_node_fitting, mem_prealloc*rbc->NNum, sizeof(CVector));
    FT_vector(&rbc->nw, mem_prealloc*rbc->NNum, sizeof(double));
    FT_vector(&rbc->la[0],mem_prealloc*rbc->NNum,sizeof(CVector));
    FT_vector(&rbc->la[1],mem_prealloc*rbc->NNum,sizeof(CVector));
    FT_vector(&rbc->node_nb,mem_prealloc*rbc->NNum*8*2,sizeof(int));

    // added by S. Xu
    FT_vector(&rbc->node0, mem_prealloc*rbc->NNum, sizeof(CVector));
    FT_vector(&rbc->nodep, mem_prealloc*rbc->NNum, sizeof(CVector));
    FT_vector(&rbc->nodeVel, mem_prealloc*rbc->NNum, sizeof(CVector));
    FT_vector(&rbc->fn, mem_prealloc*rbc->NNum, sizeof(CVector));

    rbc->chem = new double[mem_prealloc*rbc->NNum];

    rbc->Curr_ele_index = new std::list<int>(mem_prealloc*rbc->ENum);
    if (rbc->Curr_ele_index == NULL){
        std::cout<<"Curr_ele_index is null"<<std::endl;
    }
    rbc->Curr_ele = new std::list<ELEMENT>(mem_prealloc*rbc->ENum);
    if (rbc->Curr_ele == NULL){
        std::cout<<"Curr_ele is null"<<std::endl;
    }
    /*FT_vector(&rbc->fund_geom, rbc->NNum, sizeof(FUND_GEOM));*/
    //END:: added by S. Xu
    std::list<ELEMENT>::iterator it_ele;
    
    it_ele = rbc->Curr_ele->begin();
    for(int k = 0; it_ele != rbc->Curr_ele->end(); it_ele++, k++){
        if (k >= rbc->ENum){
            break;
        }
        it_ele->affine_trans = NULL;
        it_ele->inv_affine_trans = NULL;
        it_ele->mass_matrix = NULL;
        it_ele->mass_inv = NULL;
    }
    
    rbc->sv_it_adj_ele = Allocate_2D_matrix<std::list<ELEMENT>::iterator>( mem_prealloc*rbc->ENum, 3);

    make_rbc_structure(rbc, coordInfoVecs, generalParams);
    std::cout<<"make_rbc_structure DONE!"<<std::endl;
    build_ele_adj_by_iter(rbc, coordInfoVecs);
    std::cout<<"build_ele_adj_by_iter DONE!"<<std::endl;
    Comp_ele_affine_trans(rbc, coordInfoVecs);
    std::cout<<"Comp_ele_affine_trans DONE!"<<std::endl;
    comput_aff_cent_vtx(rbc);
    std::cout<<"aff_cent_vtx DONE!"<<std::endl;
    Comp_ele_mass_matrix(rbc, coordInfoVecs);
    std::cout<<"Comp_ele_mass_matrix DONE!"<<std::endl;
    comput_aff_conormals(rbc);
    std::cout<<"comput_aff_conormals DONE!"<<std::endl;

////////////////////////////////////////////////////////////////////////////////////////
////////////////////////// END OF CREATING THE FIRST CELL //////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
    std::cout<<"START building 2nd RBC as backup"<<std::endl;
    n_rbc->NNum = generalParams.maxNodeCount;
    n_rbc->ENum = coordInfoVecs.num_triangles;    

    FT_vector(&n_rbc->nodeloc, mem_prealloc*n_rbc->NNum, sizeof(CVector));
    FT_vector(&n_rbc->ele, mem_prealloc*n_rbc->ENum, sizeof(ELEMENT));
    FT_vector(&n_rbc->node, mem_prealloc*n_rbc->NNum, sizeof(CVector));
    FT_vector(&n_rbc->normal_node_fitting, mem_prealloc*n_rbc->NNum, sizeof(CVector));
    FT_vector(&n_rbc->nw, mem_prealloc*n_rbc->NNum, sizeof(double));
    FT_vector(&n_rbc->la[0],mem_prealloc*n_rbc->NNum,sizeof(CVector));
    FT_vector(&n_rbc->la[1],mem_prealloc*n_rbc->NNum,sizeof(CVector));
    FT_vector(&n_rbc->node_nb,mem_prealloc*n_rbc->NNum*8*2,sizeof(int));

//	added by S. Xu
    FT_vector(&n_rbc->node0, mem_prealloc*n_rbc->NNum, sizeof(CVector));
    FT_vector(&n_rbc->nodep, mem_prealloc*n_rbc->NNum, sizeof(CVector));
    FT_vector(&n_rbc->nodeVel, mem_prealloc*n_rbc->NNum, sizeof(CVector));
    FT_vector(&n_rbc->fn, mem_prealloc*n_rbc->NNum, sizeof(CVector));
    

    n_rbc->Curr_ele_index = new std::list<int>(mem_prealloc*n_rbc->ENum);
    if (n_rbc->Curr_ele_index == NULL){
        std::cout<<"Curr_ele_index is null"<<std::endl;
    }
    n_rbc->Curr_ele = new std::list<ELEMENT>(mem_prealloc*n_rbc->ENum);
    if (n_rbc->Curr_ele == NULL){
        std::cout<<"Curr_ele is null"<<std::endl;
    }
    /*FT_vector(&rbc->fund_geom, rbc->NNum, sizeof(FUND_GEOM));*/
    //END:: added by S. Xu
    std::list<ELEMENT>::iterator n_it_ele;
    n_it_ele = n_rbc->Curr_ele->begin();
    for(int k = 0; n_it_ele != n_rbc->Curr_ele->end(); n_it_ele++, k++){
        if (k >= n_rbc->ENum){
            break;
        }
        n_it_ele->affine_trans = NULL;
        n_it_ele->inv_affine_trans = NULL;
        n_it_ele->mass_matrix = NULL;
        n_it_ele->mass_inv = NULL;
    }
    
    n_rbc->sv_it_adj_ele = Allocate_2D_matrix<std::list<ELEMENT>::iterator>( mem_prealloc*n_rbc->ENum, 3);

    make_rbc_structure(n_rbc, coordInfoVecs, generalParams);
    std::cout<<"2nd make_rbc_structure DONE!"<<std::endl;
    build_ele_adj_by_iter(n_rbc, coordInfoVecs);
    std::cout<<"2nd build_ele_adj_by_iter DONE!"<<std::endl;
    Comp_ele_affine_trans(n_rbc, coordInfoVecs);
    std::cout<<"2nd Comp_ele_affine_trans DONE!"<<std::endl;
    comput_aff_cent_vtx(n_rbc);
    std::cout<<"2nd aff_cent_vtx DONE!"<<std::endl;
    Comp_ele_mass_matrix(n_rbc, coordInfoVecs);
    std::cout<<"2nd Comp_ele_mass_matrix DONE!"<<std::endl;
    comput_aff_conormals(n_rbc);
    std::cout<<"2nd comput_aff_conormals DONE!"<<std::endl;
    

// 	//////////////////////////////////////////////////////////////////////////////////////
// 	/////////////////////// END OF CREATEING THE SECOND CELL /////////////////////////////
// 	//////////////////////////////////////////////////////////////////////////////////////
    it_ele = rbc->Curr_ele->begin();
    double current_soln = 0.0;
    double centroid_z, centroid_x, centroid_y;
    int initial_elem_count = 0;
    n_it_ele = n_rbc->Curr_ele->begin();
    it_ele = rbc->Curr_ele->begin();
    int k = 0;
    for (it_ele = rbc->Curr_ele->begin(); (it_ele != rbc->Curr_ele->end()); it_ele++, n_it_ele++, k++){
        if (k >= rbc->ENum){
            break;
        }
        if (it_ele->vertex[0] >= (INT_MAX-100) || it_ele->vertex[1] >= (INT_MAX-100) || it_ele->vertex[2] >= (INT_MAX-100)
            || it_ele->vertex[0] <= (-INT_MAX+100) || it_ele->vertex[1] <= (-INT_MAX+100) || it_ele->vertex[2] <= (-INT_MAX+100)){
            it_ele->dg_rho[0][0] = INT_MAX;
            it_ele->dg_rho[0][1] = INT_MAX;
            it_ele->dg_rho[0][2] = INT_MAX;
            it_ele->dg_rho[1][0] = INT_MAX;
            it_ele->dg_rho[1][1] = INT_MAX;
            it_ele->dg_rho[1][2] = INT_MAX;
            
        }
        if (n_it_ele->vertex[0] >= (INT_MAX-100) || n_it_ele->vertex[1] >= (INT_MAX-100) || n_it_ele->vertex[2] >= (INT_MAX-100)
            || it_ele->vertex[0] <= (-INT_MAX+100) || n_it_ele->vertex[1] <= (-INT_MAX+100) || n_it_ele->vertex[2] <= (-INT_MAX+100)){
            n_it_ele->dg_rho[0][0] = INT_MAX;
            n_it_ele->dg_rho[0][1] = INT_MAX;
            n_it_ele->dg_rho[0][2] = INT_MAX;
            n_it_ele->dg_rho[1][0] = INT_MAX;
            n_it_ele->dg_rho[1][1] = INT_MAX;
            n_it_ele->dg_rho[1][2] = INT_MAX;
            
        }
    }

    std::list<ELEMENT>::iterator it_area;
    double chem_area = 0.0;
    it_area = rbc->Curr_ele->begin();
    for (int i=0; i<rbc->ENum; i++, it_area++){
        if (coordInfoVecs.triangles2Nodes_1[i] >= (INT_MAX-100) || coordInfoVecs.triangles2Nodes_2[i] >= (INT_MAX-100) || coordInfoVecs.triangles2Nodes_3[i] >= (INT_MAX-100)){
            continue;
        }
        if (coordInfoVecs.triangles2Nodes_1[i] <= (-INT_MAX+100) || coordInfoVecs.triangles2Nodes_2[i] <= (-INT_MAX+100) || coordInfoVecs.triangles2Nodes_3[i] <= (-INT_MAX+100)){
            continue;
        }
        centroid_x = (coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_1[i]] +
            coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_2[i]] + 
            coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_3[i]])/3.0;
        centroid_y = (coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_1[i]] +
            coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_2[i]] + 
            coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_3[i]])/3.0;
        centroid_z = (coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_1[i]] +
            coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_2[i]] + 
            coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_3[i]])/3.0;
        if (centroid_x >= -10000.0){//3.85){// && centroid_z <= 4.0 && centroid_x >= -0.36 && centroid_x <= 0.64 && centroid_y >= 0.6 && centroid_y <= 1.0){
            initial_elem_count += 1;
            chem_area += it_area->area;
            if (it_area->Id == 0 || it_area->Id == 30){std::cout<<"area of triangle["<<it_area->Id<<"] ="<<it_area->area<<std::endl;}
        }
        
    }
    std::cout<<"initial_elem_count = "<<initial_elem_count<<std::endl;
    //std::cout<<"initially assigned conc @ the above elem = "<<1.0/initial_elem_count<<std::endl;
    it_ele = rbc->Curr_ele->begin();
    for (int i=0; (it_ele != rbc->Curr_ele->end()); it_ele++, i++){
        if (i >= rbc->ENum){
            break;
        }
        
        if (coordInfoVecs.triangles2Nodes_1[i] >= (INT_MAX-100) || coordInfoVecs.triangles2Nodes_2[i] >= (INT_MAX-100) || coordInfoVecs.triangles2Nodes_3[i] >= (INT_MAX-100)){
            continue;
        }
        if (coordInfoVecs.triangles2Nodes_1[i] <= (-INT_MAX+100) || coordInfoVecs.triangles2Nodes_2[i] <= (-INT_MAX+100) || coordInfoVecs.triangles2Nodes_3[i] <= (-INT_MAX+100)){
            continue;
        }
        centroid_x = (coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_1[i]] +
            coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_2[i]] + 
            coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_3[i]])/3.0;
        centroid_y = (coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_1[i]] +
            coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_2[i]] + 
            coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_3[i]])/3.0;
        centroid_z = (coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_1[i]] +
            coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_2[i]] + 
            coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_3[i]])/3.0;
                        
        if (centroid_x >= -1000.0){//3.85){// && centroid_z <= 4.0 && centroid_x >= -0.36 && centroid_x <= 0.64 && centroid_y >= 0.6 && centroid_y <= 1.0){
            it_ele->dg_rho[0][0] = 0.04;//1.0/initial_elem_count;//1.0/1000;//rbc->ENum;//1.0;//0.0;
            it_ele->dg_rho[0][1] = 0.0;
            it_ele->dg_rho[0][2] = 0.0;
        }
        else{
            it_ele->dg_rho[0][0] = 0.0;
            it_ele->dg_rho[0][1] = 0.0;
            it_ele->dg_rho[0][2] = 0.0;
        }
        it_ele->b[0][0] = 0.01;
        it_ele->b[0][1] = 0.0;
        it_ele->b[0][2] = 0.0;
        // std::cout<<it_ele->dg_rho[0][0]<<" "<<it_ele->dg_rho[0][1]<<" "<<it_ele->dg_rho[0][2]<<std::endl;
        // std::cout<<it_ele->dg_rho[1][0]<<" "<<it_ele->dg_rho[1][1]<<" "<<it_ele->dg_rho[1][2]<<std::endl;
        
    }

    n_it_ele = n_rbc->Curr_ele->begin();
    int ww = 0;
    for (it_ele = rbc->Curr_ele->begin();(it_ele != rbc->Curr_ele->end()); it_ele++, n_it_ele++, ww++){
        if (ww >= coordInfoVecs.num_triangles){
            break;
        }
        // std::cout<<"initial dg["<<ww<<"] = "<<it_ele->dg_rho[0][0]<<" "<<it_ele->dg_rho[0][1]<<" "<<it_ele->dg_rho[0][2]<<" and b["<<ww<<"] = "<<it_ele->b[0][0]<<" "<<it_ele->b[0][1]<<" "<<it_ele->b[0][2]<<std::endl;
    }
    double diff_error[rbc->ENum];
    //  utilities_ptr->LDG_Error_of_surf_diffusion_test(rbc,0.0,diff_error);
}
void Utilities::LDG_Surface_Diffusion_Solve(
    double edgeswap_iteration,
    CoordInfoVecs& coordInfoVecs,
    GeneralParams& generalParams,
    HostSetInfoVecs& hostSetInfoVecs,
    AuxVecs& auxVecs,
	RBC      *rbc,
    RBC      *n_rbc){     
        std::list<ELEMENT>::iterator it_ele;
        std::list<ELEMENT>::iterator n_it_ele;       

        double dg_q[2][10], *dg_q_p[2]; 
        dg_q_p[0] = dg_q[0]; 
        dg_q_p[1] = dg_q[1]; 
        for(int i = 0; i < 2; i++)
        {
            for(int k = 0; k < MAX_N_COEF; k++)
                dg_q[i][k] = 0.0;
        }
        // std::list<ELEMENT>::iterator ele_it;
        // for(ele_it = rbc->Curr_ele->begin();(ele_it!=rbc->Curr_ele->end());ele_it++){
        // 	// utilities_ptr->accurate_sphere_diffusion_initializer(rbc, ele_it, ele_it->dg_rho[0], dg_q_p); 
        // }

        std::cout<<"before solving PDE"<<std::endl;
        current_area(rbc);
        double t = 0.0;
        double time_step;
        double Dconst = 0.005;//0.05;
        std::cout<<"Dconst = "<<Dconst<<std::endl;
        // double sink_coef = 0.1;
        // std::cout<<"sink_coef = "<<sink_coef<<std::endl;
        // double source_coef = 0.0001;
        // std::cout<<"source_coef = "<<source_coef<<std::endl;
        int maxruntime;
        // if (edgeswap_iteration <0){
        if (edgeswap_iteration <2){
            maxruntime = 1.2e4;//12e4;//2e3+1;
        }
        else{
            maxruntime = 1e3;//12e4;//2e3+1;
        }
        std::cout<<"maxruntime = "<<maxruntime<<std::endl;
        double time_step_actual = 5e-4;//2e-5;
        std::cout<<"time step = "<<time_step_actual<<std::endl;
        int RECORD_FREQUENCY = floor(maxruntime/10);
        std::cout<<"RECORD_FREQUENCY = "<<RECORD_FREQUENCY<<std::endl;
        
        it_ele = rbc->Curr_ele->begin();
        std::list<ELEMENT>::iterator it_record;
        std::vector<double> num_soln(rbc->ENum, INT_MAX);
        std::vector<double> b_soln(rbc->ENum, INT_MAX);
        // double k_0 = 20.0;
        // double k_1 = 25.0;
        // double k_2 = 5.0;
        // double k_3 = 5.0;
        // double k_4 = 1.0;
        // double k_ss = 12;//10.75;
        // double beta = 1.0/1.0;///1.45;
        // double gamma = 1.0;
        // double q1 = 10.0;
        // double h = 10.0;
        
        if (edgeswap_iteration < 1){
            std::cout<<"k_0 = "<<coordInfoVecs.k_0<<std::endl;
            std::cout<<"k_1 = "<<coordInfoVecs.k_1<<std::endl;
            std::cout<<"k_2 = "<<coordInfoVecs.k_2<<std::endl;
            std::cout<<"k_3 = "<<coordInfoVecs.k_3<<std::endl;
            std::cout<<"k_4 = "<<coordInfoVecs.k_4<<std::endl;
            std::cout<<"k_ss = "<<coordInfoVecs.k_ss<<std::endl;
            std::cout<<"beta = "<<coordInfoVecs.beta<<std::endl;
            std::cout<<"gamma = "<<coordInfoVecs.gamma<<std::endl;
            std::cout<<"q = "<<coordInfoVecs.q1<<std::endl;
            std::cout<<"h = "<<coordInfoVecs.h<<std::endl;
        }
        
        std::list<ELEMENT>::iterator it_sum;
        double total_chem = 0.0;
        double total_area = 0.0;
        for (int i = 0; i < maxruntime ; i++){
            generalParams.current_total_sim_step += 1;
            // std::cout<<" "<<std::endl;
            // std::cout<<"i = "<<i<<std::endl;
            if (i == 0){
                time_step = 0.0;
            }
            else{
                time_step = time_step_actual;
                // if (i == 1){
                // 	std::cout<<"time step = "<<time_step<<std::endl;
                // }
            }
            // utilities_ptr->LDG_Surf_diffusion(rbc, n_rbc, Dconst, sink_coef, source_coef, t, 1.0, time_step, 0);//This time step, 1e-4, may be too big for stability, check FLcondition.
            LDG_Surf_diffusion(rbc, n_rbc, Dconst, coordInfoVecs.k_0, 
                                coordInfoVecs.k_1, coordInfoVecs.k_2, coordInfoVecs.k_3,
                                coordInfoVecs.k_4, coordInfoVecs.k_ss, coordInfoVecs.beta,
                                coordInfoVecs.gamma, coordInfoVecs.q1,coordInfoVecs.h,
                                t, 1.0, time_step, 0);

            total_chem = 0.0;
            total_area = 0.0;
            int v = 0;
            // #pragma omp parallel
            {
                // std::cout<<"OpenMP thread number (?) = "<<omp_get_num_threads()<<std::endl;
                for (it_sum=rbc->Curr_ele->begin(); it_sum != rbc->Curr_ele->end(); it_sum++, v++){
                    if (v >= rbc->ENum){
                        break;
                    }    
                    if (it_sum->vertex[0] >= (INT_MAX-100) || it_sum->vertex[1] >= (INT_MAX-100) || it_sum->vertex[2] >= (INT_MAX-100)
                    || it_sum->vertex[0] <= (-INT_MAX+100) || it_sum->vertex[1] <= (-INT_MAX+100) || it_sum->vertex[2] <= (-INT_MAX+100)){
                        continue;
                    }
                
                    soln_at_pt(it_sum->dg_rho[0], it_sum->aff_cent, it_sum->aff_cent, sqrt(it_sum->area), &num_soln[v]);
                    soln_at_pt(it_sum->b[0], it_sum->aff_cent, it_sum->aff_cent, sqrt(it_sum->area), &b_soln[v]);
                    total_chem += num_soln[v];
                    total_area += it_sum->area;
                }
            }
            if (i == 0){
                std::cout<<"cell surface area = "<<total_area<<std::endl;
            }
            int vv = 0;
            // #pragma omp parallel
            {
                for (it_sum=rbc->Curr_ele->begin(); it_sum != rbc->Curr_ele->end(); it_sum++, vv++){
                    if (vv >= rbc->ENum){
                        break;
                    }
                    it_sum->b[0][0] = it_sum->b[0][0] + time_step*coordInfoVecs.k_4*(total_chem/total_area - coordInfoVecs.k_ss)*it_sum->b[0][0];
                    it_sum->b[0][1] = 0.0;
                    it_sum->b[0][2] = 0.0;
                }
            }

            t += time_step;
            // current_soln = 0.0;
            n_it_ele = n_rbc->Curr_ele->begin();
            int w = 0;

            // #pragma omp parallel
            {
                for (it_ele = rbc->Curr_ele->begin(); (it_ele != rbc->Curr_ele->end()); it_ele++, n_it_ele++, w++){
                    if (w >= rbc->ENum){
                        break;
                    }

                    if (it_ele->vertex[0] >= (INT_MAX-100) || it_ele->vertex[1] >= (INT_MAX-100) || it_ele->vertex[2] >= (INT_MAX-100)
                    || it_ele->vertex[0] <= (-INT_MAX+100) || it_ele->vertex[1] <= (-INT_MAX+100) || it_ele->vertex[2] <= (-INT_MAX+100)){
                        continue;
                    }
                    //utilities_ptr->soln_at_pt(it_ele->dg_rho[0], rbc->nodes[i], rbc->Center,)
                    if (i == (maxruntime-1)){
                        soln_at_pt(it_ele->dg_rho[0], it_ele->aff_cent, it_ele->aff_cent, sqrt(it_ele->area), &num_soln[w]);
                        //coordInfoVecs.soln_per_triangle[w] = num_soln[i];
                        //std::cout<<w<<" "<<num_soln[w]<<std::endl;
                        
                    }
                    //  current_soln += it_ele->dg_rho[0][0];
                    if ((i%200 == 0 || i%200 == 1 )&& (it_ele->Id == 1 || it_ele->Id == 5)){
                        
                    }
                    
                }
            }
            //if (i == 0 || i == (maxruntime-1)){
            
            if (RECORD_FREQUENCY == 0){
                if (i == (maxruntime-1) || i == (maxruntime)){
                    int v = 0;
                    total_chem = 0.0;
                    double total_b_chem = 0.0;
                    // #pragma omp parallel
                    {
                        for (it_record=rbc->Curr_ele->begin(); it_record != rbc->Curr_ele->end(); it_record++, v++){
                            if (v >= rbc->ENum){
                                break;
                            }
                            if (it_record->vertex[0] >= (INT_MAX-100) || it_record->vertex[1] >= (INT_MAX-100) || it_record->vertex[2] >= (INT_MAX-100)
                            || it_record->vertex[0] <= (-INT_MAX+100) || it_record->vertex[1] <= (-INT_MAX+100) || it_record->vertex[2] <= (-INT_MAX+100)){
                                continue;
                            }
                        
                            soln_at_pt(it_record->dg_rho[0], it_record->aff_cent, it_record->aff_cent, sqrt(it_record->area), &num_soln[v]);
                            soln_at_pt(it_record->b[0], it_record->aff_cent, it_record->aff_cent, sqrt(it_record->area), &b_soln[v]);
                            coordInfoVecs.soln_per_triangle[v] = num_soln[v];
                            coordInfoVecs.b_per_triangle[v] = b_soln[v];
                            total_chem += num_soln[v];
                            total_b_chem += b_soln[v];
                        
                        }
                    }
                    // std::cout<<"current computation time = "<<i<<"*"<<time_step<<std::endl;
                    // std::cout<<"current total chem conc = "<<total_chem<<std::endl;
                    // std::cout<<"current total b conc = "<<total_b_chem<<std::endl;
                    std::cout<<"Total number of sim steps = "<<generalParams.current_total_sim_step<<" current total chem conc = "<<total_chem<<" current total b conc = "<<total_b_chem<<std::endl;
                    //storage->storeVariables();
                    
                }
            }
            else{
                if (i == (maxruntime-1) || i == (maxruntime)){//i%RECORD_FREQUENCY == 0
                    int v = 0;
                    total_chem = 0.0;
                    double total_b_chem = 0.0;
                    // #pragma omp parallel
                    {
                        for (it_record=rbc->Curr_ele->begin(); it_record != rbc->Curr_ele->end(); it_record++, v++){
                            if (v >= rbc->ENum){
                                break;
                            }
                            if (it_record->vertex[0] >= (INT_MAX-100) || it_record->vertex[1] >= (INT_MAX-100) || it_record->vertex[2] >= (INT_MAX-100)
                            || it_record->vertex[0] <= (-INT_MAX+100) || it_record->vertex[1] <= (-INT_MAX+100) || it_record->vertex[2] <= (-INT_MAX+100)){
                                continue;
                            }
                        
                            soln_at_pt(it_record->dg_rho[0], it_record->aff_cent, it_record->aff_cent, sqrt(it_record->area), &num_soln[v]);
                            soln_at_pt(it_record->b[0], it_record->aff_cent, it_record->aff_cent, sqrt(it_record->area), &b_soln[v]);
                            coordInfoVecs.soln_per_triangle[v] = num_soln[v];
                            coordInfoVecs.b_per_triangle[v] = b_soln[v];
                            total_chem += num_soln[v];
                            total_b_chem += b_soln[v];
                        
                        }
                    }
                    
                    std::cout<<"Total number of sim steps = "<<generalParams.current_total_sim_step<<" current total chem conc = "<<total_chem<<" current total b conc = "<<total_b_chem<<std::endl;
                    // storage->storeVariables();
                    
                    
                    // utilities_ptr->LDG_Error_of_surf_diffusion_test(rbc,t,diff_error);
                    //  std::cout<<"current_conc = "<<current_soln<<std::endl;
                }
            }
            if (i == (maxruntime-1)){
                std::cout<<"End of the chem diff simulation"<<std::endl;
            }
        }
        // n_it_ele = n_rbc->Curr_ele->begin();
        // for (it_ele = rbc->Curr_ele->begin();(it_ele != rbc->Curr_ele->end()); it_ele++, n_it_ele++){
            // utilities_ptr->soln_at_pt(it_ele->dg_rho[0], rbc->nodes[i], rbc->Center,)
            // std::cout<<it_ele->dg_rho[0][0]<<" , "<<it_ele->dg_rho[0][1]<<" , "<<it_ele->dg_rho[0][2]<<" , "<<std::endl;
            // std::cout<<"it dg[1] = "<<it_ele->dg_rho[1][0]<<" "<<it_ele->dg_rho[1][1]<<" "<<it_ele->dg_rho[1][2]<<std::endl;
            // std::cout<<"n_it dg[0] = "<<n_it_ele->dg_rho[0][0]<<" "<<n_it_ele->dg_rho[0][1]<<" "<<n_it_ele->dg_rho[0][2]<<std::endl;
            // std::cout<<"n_it dg[1] = "<<n_it_ele->dg_rho[1][0]<<" "<<n_it_ele->dg_rho[1][1]<<" "<<n_it_ele->dg_rho[1][2]<<std::endl;
            // break;
            // std::cout<<it_ele->dg_q[0][0][0]<<" "<<it_ele->dg_q[0][0][1]<<" "<<it_ele->dg_q[0][0][2]<<std::endl;
            // std::cout<<it_ele->dg_q[0][1][0]<<" "<<it_ele->dg_q[0][1][1]<<" "<<it_ele->dg_q[0][1][2]<<std::endl;
            // std::cout<<n_it_ele->dg_q[1][0][0]<<" "<<n_it_ele->dg_q[1][0][1]<<" "<<n_it_ele->dg_q[1][0][2]<<std::endl;
            // std::cout<<n_it_ele->dg_q[1][1][0]<<" "<<n_it_ele->dg_q[1][1][1]<<" "<<n_it_ele->dg_q[1][1][2]<<std::endl;
        // }
        // double diff_error[rbc->ENum];
        // utilities_ptr->LDG_Error_of_surf_diffusion_test(rbc,0.005,diff_error);
    }

    void Utilities::LDG_Surface_Diffusion_Structure_Rebuild(
    CoordInfoVecs& coordInfoVecs,
    GeneralParams& generalParams,
    HostSetInfoVecs& hostSetInfoVecs,
    AuxVecs& auxVecs,
	RBC      *rbc,
    RBC      *n_rbc
    ){
        transferDtoH(generalParams, coordInfoVecs, hostSetInfoVecs);
        for (int i = 0; i < coordInfoVecs.num_triangles; i++){
            triangles2Triangles_host_vecs(i, hostSetInfoVecs,coordInfoVecs,generalParams, auxVecs);
        }
        transferHtoD(generalParams,coordInfoVecs,hostSetInfoVecs);

        rbc->NNum = generalParams.maxNodeCount;
        rbc->ENum = coordInfoVecs.num_triangles;    

    //     FT_vector(&rbc->nodeloc, rbc->NNum, sizeof(CVector));
    //     FT_vector(&rbc->ele, rbc->ENum, sizeof(ELEMENT));
    //     FT_vector(&rbc->node, rbc->NNum, sizeof(CVector));
    //     FT_vector(&rbc->normal_node_fitting, rbc->NNum, sizeof(CVector));
    //     FT_vector(&rbc->nw, rbc->NNum, sizeof(double));
    //     FT_vector(&rbc->la[0],rbc->NNum,sizeof(CVector));
    //     FT_vector(&rbc->la[1],rbc->NNum,sizeof(CVector));
    //     FT_vector(&rbc->node_nb,rbc->NNum*8*2,sizeof(int));

    //     // added by S. Xu
    //     FT_vector(&rbc->node0, rbc->NNum, sizeof(CVector));
    //     FT_vector(&rbc->nodep, rbc->NNum, sizeof(CVector));
    //     FT_vector(&rbc->nodeVel, rbc->NNum, sizeof(CVector));
    //     FT_vector(&rbc->fn, rbc->NNum, sizeof(CVector));
        // rbc->Curr_ele_index= new std::list<int>(rbc->ENum);
        // if (rbc->Curr_ele_index == NULL){
            // std::cout<<"Curr_ele_index is null"<<std::endl;
        // }
        
        /*FT_vector(&rbc->fund_geom, rbc->NNum, sizeof(FUND_GEOM));*/
        //END:: added by S. Xu
        std::list<ELEMENT>::iterator it_ele;
        
        it_ele = rbc->Curr_ele->begin();
        for(int k = 0; it_ele != rbc->Curr_ele->end(); it_ele++, k++){
            if (k >= rbc->ENum){
                break;
            }
            it_ele->affine_trans = NULL;
            it_ele->inv_affine_trans = NULL;
            it_ele->mass_matrix = NULL;
            it_ele->mass_inv = NULL;
        }
        
        // rbc->sv_it_adj_ele = Allocate_2D_matrix<std::list<ELEMENT>::iterator>( rbc->ENum, 3);

        make_rbc_structure(rbc, coordInfoVecs, generalParams);
        // std::cout<<"make_rbc_structure rebuild DONE!"<<std::endl;
        build_ele_adj_by_iter(rbc, coordInfoVecs);
        // std::cout<<"build_ele_adj_by_iter rebuild DONE!"<<std::endl;
        Comp_ele_affine_trans(rbc, coordInfoVecs);
        // std::cout<<"Comp_ele_affine_trans rebuild DONE!"<<std::endl;
        comput_aff_cent_vtx(rbc);
        // std::cout<<"aff_cent_vtx rebuild DONE!"<<std::endl;
        Comp_ele_mass_matrix(rbc, coordInfoVecs);
        // std::cout<<"Comp_ele_mass_matrix rebuild DONE!"<<std::endl;
        comput_aff_conormals(rbc);
        // std::cout<<"comput_aff_conormals rebuild DONE!"<<std::endl;
        std::cout<<"RBC rebuild done!"<<std::endl;

    // ////////////////////////////////////////////////////////////////////////////////////////
    // ////////////////////////// END OF CREATING THE FIRST CELL //////////////////////////////
    // ////////////////////////////////////////////////////////////////////////////////////////
    //     std::cout<<"START building 2nd RBC as backup"<<std::endl;
        n_rbc->NNum = generalParams.maxNodeCount;
        n_rbc->ENum = coordInfoVecs.num_triangles;    

    //     FT_vector(&n_rbc->nodeloc, n_rbc->NNum, sizeof(CVector));
    //     FT_vector(&n_rbc->ele, n_rbc->ENum, sizeof(ELEMENT));
    //     FT_vector(&n_rbc->node, n_rbc->NNum, sizeof(CVector));
    //     FT_vector(&n_rbc->normal_node_fitting, n_rbc->NNum, sizeof(CVector));
    //     FT_vector(&n_rbc->nw, n_rbc->NNum, sizeof(double));
    //     FT_vector(&n_rbc->la[0],n_rbc->NNum,sizeof(CVector));
    //     FT_vector(&n_rbc->la[1],n_rbc->NNum,sizeof(CVector));
    //     FT_vector(&n_rbc->node_nb,n_rbc->NNum*8*2,sizeof(int));

    // //	added by S. Xu
    //     FT_vector(&n_rbc->node0, n_rbc->NNum, sizeof(CVector));
    //     FT_vector(&n_rbc->nodep, n_rbc->NNum, sizeof(CVector));
    //     FT_vector(&n_rbc->nodeVel, n_rbc->NNum, sizeof(CVector));
    //     FT_vector(&n_rbc->fn, n_rbc->NNum, sizeof(CVector));
    //     rbc->chem = new double[rbc->NNum];

        
        // n_rbc->Curr_ele_index = new std::list<int>(n_rbc->ENum);
        // if (n_rbc->Curr_ele_index == NULL){
        //     std::cout<<"Curr_ele_index is null"<<std::endl;
        // }
        
        /*FT_vector(&rbc->fund_geom, rbc->NNum, sizeof(FUND_GEOM));*/
        //END:: added by S. Xu
        std::list<ELEMENT>::iterator n_it_ele;
        n_it_ele = n_rbc->Curr_ele->begin();
        for(int k = 0; n_it_ele != n_rbc->Curr_ele->end(); n_it_ele++, k++){
            if (k >= n_rbc->ENum){
                break;
            }
            n_it_ele->affine_trans = NULL;
            n_it_ele->inv_affine_trans = NULL;
            n_it_ele->mass_matrix = NULL;
            n_it_ele->mass_inv = NULL;
        }
        
        // n_rbc->sv_it_adj_ele = Allocate_2D_matrix<std::list<ELEMENT>::iterator>( n_rbc->ENum, 3);

        make_rbc_structure(n_rbc, coordInfoVecs, generalParams);
        // std::cout<<"2nd make_rbc_structure rebuild DONE!"<<std::endl;
        build_ele_adj_by_iter(n_rbc, coordInfoVecs);
        // std::cout<<"2nd build_ele_adj_by_iter rebuild DONE!"<<std::endl;
        Comp_ele_affine_trans(n_rbc, coordInfoVecs);
        // std::cout<<"2nd Comp_ele_affine_trans rebuild DONE!"<<std::endl;
        comput_aff_cent_vtx(n_rbc);
        // std::cout<<"2nd aff_cent_vtx rebuild DONE!"<<std::endl;
        Comp_ele_mass_matrix(n_rbc, coordInfoVecs);
        // std::cout<<"2nd Comp_ele_mass_matrix rebuild DONE!"<<std::endl;
        comput_aff_conormals(n_rbc);
        // std::cout<<"2nd comput_aff_conormals rebuild DONE!"<<std::endl;
        std::cout<<"N_RBC rebuild done!"<<std::endl;
        

    // 	//////////////////////////////////////////////////////////////////////////////////////
    // 	/////////////////////// END OF CREATEING THE SECOND CELL /////////////////////////////
    // 	//////////////////////////////////////////////////////////////////////////////////////
    }

    void Utilities::LDG_Surface_Diffusion_Structure_Rebuild_postGrowth(
    CoordInfoVecs& coordInfoVecs,
    GeneralParams& generalParams,
    HostSetInfoVecs& hostSetInfoVecs,
    AuxVecs& auxVecs,
	RBC      *rbc,
    RBC      *n_rbc
    ){
        transferDtoH(generalParams, coordInfoVecs, hostSetInfoVecs);
        for (int i = 0; i < coordInfoVecs.num_triangles; i++){
            triangles2Triangles_host_vecs(i, hostSetInfoVecs,coordInfoVecs,generalParams, auxVecs);
        }
        transferHtoD(generalParams,coordInfoVecs,hostSetInfoVecs);

        rbc->NNum = generalParams.maxNodeCount;
        rbc->ENum = coordInfoVecs.num_triangles;    

        // rbc->Curr_ele_index= new std::list<int>(rbc->ENum);
        // if (rbc->Curr_ele_index == NULL){
            // std::cout<<"Curr_ele_index is null"<<std::endl;
        // }
        
        /*FT_vector(&rbc->fund_geom, rbc->NNum, sizeof(FUND_GEOM));*/
        //END:: added by S. Xu
        std::list<ELEMENT>::iterator it_ele;
        
        it_ele = rbc->Curr_ele->begin();
        for(int k = 0; it_ele != rbc->Curr_ele->end(); it_ele++, k++){
            if (k >= rbc->ENum){
                break;
            }
            it_ele->affine_trans = NULL;
            it_ele->inv_affine_trans = NULL;
            it_ele->mass_matrix = NULL;
            it_ele->mass_inv = NULL;
        }
        
        // rbc->sv_it_adj_ele = Allocate_2D_matrix<std::list<ELEMENT>::iterator>( rbc->ENum, 3);

        make_rbc_structure(rbc, coordInfoVecs, generalParams);
        // std::cout<<"make_rbc_structure rebuild DONE!"<<std::endl;
        build_ele_adj_by_iter(rbc, coordInfoVecs);
        // std::cout<<"build_ele_adj_by_iter rebuild DONE!"<<std::endl;
        Comp_ele_affine_trans(rbc, coordInfoVecs);
        // std::cout<<"Comp_ele_affine_trans rebuild DONE!"<<std::endl;
        comput_aff_cent_vtx(rbc);
        // std::cout<<"aff_cent_vtx rebuild DONE!"<<std::endl;
        Comp_ele_mass_matrix(rbc, coordInfoVecs);
        // std::cout<<"Comp_ele_mass_matrix rebuild DONE!"<<std::endl;
        comput_aff_conormals(rbc);
        // std::cout<<"comput_aff_conormals rebuild DONE!"<<std::endl;
        std::cout<<"RBC rebuild done!"<<std::endl;

        std::vector<double> new_concentration(3,0.0) , new_concentration_b(3,0.0);
        double new_area_from_growth = 0.0;
        it_ele = rbc->Curr_ele->begin();
        for (int i=0; (it_ele != rbc->Curr_ele->end()); it_ele++, i++){
            if (i >= rbc->ENum){
                break;
            }
            if (i == coordInfoVecs.num_triangles-4 ||
                 i == coordInfoVecs.num_triangles-3 ||
                 i == coordInfoVecs.num_triangles-2 ||
                 i == coordInfoVecs.num_triangles-1){
                     new_area_from_growth += it_ele->area0;
                 }
        }
        std::cout<<"new area gained = "<<new_area_from_growth<<std::endl;
        
        it_ele = rbc->Curr_ele->begin();
        for (int i=0; (it_ele != rbc->Curr_ele->end()); it_ele++, i++){
            if (i >= rbc->ENum){
                break;
            }
            if (generalParams.triangle_undergoing_growth.size() != 2){
                break;
            }
            if (i == generalParams.triangle_undergoing_growth[0] || i == generalParams.triangle_undergoing_growth[1]){
                new_concentration[0] += it_ele->dg_rho[0][0];
                new_concentration[1] += it_ele->dg_rho[0][1];
                new_concentration[2] += it_ele->dg_rho[0][2];
                new_concentration_b[0] += it_ele->b[0][0];
                new_concentration_b[1] += it_ele->b[0][1];
                new_concentration_b[2] += it_ele->b[0][2];
            }
            
            if (i == coordInfoVecs.num_triangles-4){
                new_concentration[0] = new_concentration[0]/2.0;
                new_concentration[1] = new_concentration[1]/2.0;
                new_concentration[2] = new_concentration[2]/2.0;
                new_concentration_b[0] = new_concentration_b[0]/2.0;
                new_concentration_b[1] = new_concentration_b[1]/2.0;
                new_concentration_b[2] = new_concentration_b[2]/2.0;
            } // This is for NON-interpolating case


            if (i == coordInfoVecs.num_triangles-4 ||
                 i == coordInfoVecs.num_triangles-3 ||
                 i == coordInfoVecs.num_triangles-2 ||
                 i == coordInfoVecs.num_triangles-1){
                     it_ele->dg_rho[0][0] = new_concentration[0];
                    it_ele->dg_rho[0][1] = new_concentration[1];
                    it_ele->dg_rho[0][2] = new_concentration[2];
                    // This is for NON-terpolating case

                    // it_ele->dg_rho[0][0] = new_concentration[0]*(it_ele->area0)/new_area_from_growth;//1.0/initial_elem_count;//1.0/1000;//rbc->ENum;//1.0;//0.0;
                    // it_ele->dg_rho[0][1] = new_concentration[1]*(it_ele->area0)/new_area_from_growth;//0.0;
                    // it_ele->dg_rho[0][2] = new_concentration[2]*(it_ele->area0)/new_area_from_growth;//0.0;
                    // // This is for interpolating case

                    it_ele->b[0][0] = new_concentration_b[0];
                    it_ele->b[0][1] = new_concentration_b[1];//0.0;
                    it_ele->b[0][2] = new_concentration_b[2];//0.0;
                    // This is for NON-interpolating case

                    // it_ele->b[0][0] = new_concentration_b[0]/2.0;
                    // it_ele->b[0][1] = new_concentration_b[1]/2.0;
                    // it_ele->b[0][2] = new_concentration_b[2]/2.0;
                    // // This is for interpolating case
                 }
            // std::cout<<it_ele->dg_rho[0][0]<<" "<<it_ele->dg_rho[0][1]<<" "<<it_ele->dg_rho[0][2]<<std::endl;
            // std::cout<<it_ele->dg_rho[1][0]<<" "<<it_ele->dg_rho[1][1]<<" "<<it_ele->dg_rho[1][2]<<std::endl;
        }

    // ////////////////////////////////////////////////////////////////////////////////////////
    // ////////////////////////// END OF CREATING THE FIRST CELL //////////////////////////////
    // ////////////////////////////////////////////////////////////////////////////////////////
    //     std::cout<<"START building 2nd RBC as backup"<<std::endl;
        n_rbc->NNum = generalParams.maxNodeCount;
        n_rbc->ENum = coordInfoVecs.num_triangles;    
        
        // n_rbc->Curr_ele_index = new std::list<int>(n_rbc->ENum);
        // if (n_rbc->Curr_ele_index == NULL){
            // std::cout<<"Curr_ele_index is null"<<std::endl;
        // }
        
        /*FT_vector(&rbc->fund_geom, rbc->NNum, sizeof(FUND_GEOM));*/
        //END:: added by S. Xu
        std::list<ELEMENT>::iterator n_it_ele;
        n_it_ele = n_rbc->Curr_ele->begin();
        for(int k = 0; n_it_ele != n_rbc->Curr_ele->end(); n_it_ele++, k++){
            if (k >= n_rbc->ENum){
                break;
            }
            n_it_ele->affine_trans = NULL;
            n_it_ele->inv_affine_trans = NULL;
            n_it_ele->mass_matrix = NULL;
            n_it_ele->mass_inv = NULL;
        }
        
        // n_rbc->sv_it_adj_ele = Allocate_2D_matrix<std::list<ELEMENT>::iterator>( n_rbc->ENum, 3);

        make_rbc_structure(n_rbc, coordInfoVecs, generalParams);
        // std::cout<<"2nd make_rbc_structure rebuild DONE!"<<std::endl;
        build_ele_adj_by_iter(n_rbc, coordInfoVecs);
        // std::cout<<"2nd build_ele_adj_by_iter rebuild DONE!"<<std::endl;
        Comp_ele_affine_trans(n_rbc, coordInfoVecs);
        // std::cout<<"2nd Comp_ele_affine_trans rebuild DONE!"<<std::endl;
        comput_aff_cent_vtx(n_rbc);
        // std::cout<<"2nd aff_cent_vtx rebuild DONE!"<<std::endl;
        Comp_ele_mass_matrix(n_rbc, coordInfoVecs);
        // std::cout<<"2nd Comp_ele_mass_matrix rebuild DONE!"<<std::endl;
        comput_aff_conormals(n_rbc);
        // std::cout<<"2nd comput_aff_conormals rebuild DONE!"<<std::endl;
        std::cout<<"N_RBC rebuild done!"<<std::endl;
        
        n_it_ele = n_rbc->Curr_ele->begin();
        for (int i=0; (n_it_ele != n_rbc->Curr_ele->end()); n_it_ele++, i++){
            if (i >= n_rbc->ENum){
                break;
            }
            if (generalParams.triangle_undergoing_growth.size() > 2){
                std::cout<<"More than one pair of triangles undergoing growth! Something went wrong"<<std::endl;
                break;
            }
            if (i == generalParams.triangle_undergoing_growth[0] || i == generalParams.triangle_undergoing_growth[1]){
                new_concentration[0] += n_it_ele->dg_rho[0][0];
                new_concentration[1] += n_it_ele->dg_rho[0][1];
                new_concentration[2] += n_it_ele->dg_rho[0][2];
                new_concentration_b[0] += n_it_ele->b[0][0];
                new_concentration_b[1] += n_it_ele->b[0][1];
                new_concentration_b[2] += n_it_ele->b[0][2];
            }
            if (i == coordInfoVecs.num_triangles-4){
                new_concentration[0] = new_concentration[0]/2.0;
                new_concentration[1] = new_concentration[1]/2.0;
                new_concentration[2] = new_concentration[2]/2.0;
                new_concentration_b[0] = new_concentration_b[0]/2.0;
                new_concentration_b[1] = new_concentration_b[1]/2.0;
                new_concentration_b[2] = new_concentration_b[2]/2.0;
            } //This is for the non-interpolating case.

            if (i == coordInfoVecs.num_triangles-4 ||
                 i == coordInfoVecs.num_triangles-3 ||
                 i == coordInfoVecs.num_triangles-2 ||
                 i == coordInfoVecs.num_triangles-1){
                    n_it_ele->dg_rho[0][0] = new_concentration[0];
                    n_it_ele->dg_rho[0][1] = new_concentration[1];
                    n_it_ele->dg_rho[0][2] = new_concentration[2];
                     // This is for non-interpolating case

                    // n_it_ele->dg_rho[0][0] = new_concentration[0]*(it_ele->area0)/new_area_from_growth;//1.0/initial_elem_count;//1.0/1000;//rbc->ENum;//1.0;//0.0;
                    // n_it_ele->dg_rho[0][1] = new_concentration[1]*(it_ele->area0)/new_area_from_growth;//0.0;
                    // n_it_ele->dg_rho[0][2] = new_concentration[2]*(it_ele->area0)/new_area_from_growth;//0.0;
                    // // This is for interpolating case
                    
                    n_it_ele->b[0][0] = new_concentration_b[0];
                    n_it_ele->b[0][1] = new_concentration_b[1];//0.0;
                    n_it_ele->b[0][2] = new_concentration_b[2];//0.0;
                    // This is for non-interpolating case

                    // n_it_ele->b[0][0] = new_concentration_b[0]/2.0;
                    // n_it_ele->b[0][1] = new_concentration_b[1]/2.0;
                    // n_it_ele->b[0][2] = new_concentration_b[2]/2.0;
                    // // This is for interpolating case
                 }
            // std::cout<<it_ele->dg_rho[0][0]<<" "<<it_ele->dg_rho[0][1]<<" "<<it_ele->dg_rho[0][2]<<std::endl;
            // std::cout<<it_ele->dg_rho[1][0]<<" "<<it_ele->dg_rho[1][1]<<" "<<it_ele->dg_rho[1][2]<<std::endl;
        }

    // 	//////////////////////////////////////////////////////////////////////////////////////
    // 	/////////////////////// END OF CREATEING THE SECOND CELL /////////////////////////////
    // 	//////////////////////////////////////////////////////////////////////////////////////
    }