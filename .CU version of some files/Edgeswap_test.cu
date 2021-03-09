#include "System.h"
#include <random>
#include "Edgeswap_test.h"
#include <math.h>


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

Edgeswap::Edgeswap(CoordInfoVecs& coordInfoVecs) {
            
    int nnode = coordInfoVecs.nodeLocX.size();
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
    nndata = nndata_temp;
};




//The goal is to perform the swap without explicitly creating a copy of the whole data structure.
//This is achieved by extracting the full info of a smaller affected system.
int Edgeswap::edge_swap_device_vecs(
    int iedge, 
    GeneralParams& generalParams,
    CoordInfoVecs& coordInfoVecs,
    LinearSpringInfoVecs& linearSpringInfoVecs,
    BendingTriangleInfoVecs& bendingTriangleInfoVecs,
    AreaTriangleInfoVecs& areaTriangleInfoVecs) {

    int alpha = 0;
        
    int HEAD,TAIL;
    int H0, T0,H1,H2,T1,T2;
    int edge_start, edge_end;
    int a1, b1, c1, a2, b2, c2;
    double temp_bend = 0.0;
    double linear_spring_constant;
    double bend_spring_constant;
        
    if ( coordInfoVecs.edges2Triangles_1[iedge] != coordInfoVecs.edges2Triangles_2[iedge]){
        H0 = coordInfoVecs.edges2Triangles_1[iedge];//index of the 1st triangle to i-th edge
        T0 = coordInfoVecs.edges2Triangles_2[iedge];//index of the 2nd triangle to i-th edge
        edge_start = coordInfoVecs.edges2Nodes_1[iedge];//index of the 1st node of i-th edge
        edge_end = coordInfoVecs.edges2Nodes_2[iedge];//index of the 2nd node of i-th edge

        a1 = coordInfoVecs.triangles2Edges_1[H0];//index of the 1st node of triangle H0
        b1 = coordInfoVecs.triangles2Edges_2[H0];//index of the 2nd node of triangle H0
        c1 = coordInfoVecs.triangles2Edges_3[H0];//index of the 3rd node of triangle H0
        
        a2 = coordInfoVecs.triangles2Edges_1[T0];//index of the 1st node of triangle T0
        b2 = coordInfoVecs.triangles2Edges_2[T0];
        c2 = coordInfoVecs.triangles2Edges_3[T0];
        
        //Now we identify the edge indices associated with the small subsystem.
        //This gives us the indices for H1, H2, T1, T2 (see the figure below).
        if (a1 != iedge && coordInfoVecs.edges2Nodes_1[a1] == edge_start){H1 = a1;}
        else if (a1 != iedge && coordInfoVecs.edges2Nodes_2[a1] == edge_start){H1 = a1;}
        else if (a1 != iedge && coordInfoVecs.edges2Nodes_1[a1] == edge_end){H2 = a1;}
        else if (a1 != iedge && coordInfoVecs.edges2Nodes_2[a1] == edge_end){H2 = a1;}

        if (b1 != iedge && coordInfoVecs.edges2Nodes_1[b1] == edge_start){H1 = b1;}
        else if (b1 != iedge && coordInfoVecs.edges2Nodes_2[b1] == edge_start){H1 = b1;}
        else if (b1 != iedge && coordInfoVecs.edges2Nodes_1[b1] == edge_end){H2 = b1;}
        else if (b1 != iedge && coordInfoVecs.edges2Nodes_2[b1] == edge_end){H2 = b1;}

        if (c1 != iedge && coordInfoVecs.edges2Nodes_1[c1] == edge_start){H1 = c1;}
        else if (c1 != iedge && coordInfoVecs.edges2Nodes_2[c1] == edge_start){H1 = c1;}
        else if (c1 != iedge && coordInfoVecs.edges2Nodes_1[c1] == edge_end){H2 = c1;}
        else if (c1 != iedge && coordInfoVecs.edges2Nodes_2[c1] == edge_end){H2 = c1;}
        
        if (a2 != iedge && coordInfoVecs.edges2Nodes_1[a2] == edge_start){T1 = a2;}
        else if (a2 != iedge && coordInfoVecs.edges2Nodes_2[a2] == edge_start){T1 = a2;}
        else if (a2 != iedge && coordInfoVecs.edges2Nodes_1[a2] == edge_end){T2 = a2;}
        else if (a2 != iedge && coordInfoVecs.edges2Nodes_2[a2] == edge_end){T2 = a2;}

        if (b2 != iedge && coordInfoVecs.edges2Nodes_1[b2] == edge_start){T1 = b2;}
        else if (b2 != iedge && coordInfoVecs.edges2Nodes_2[b2] == edge_start){T1 = b2;}
        else if (b2 != iedge && coordInfoVecs.edges2Nodes_1[b2] == edge_end){T2 = b2;}
        else if (b2 != iedge && coordInfoVecs.edges2Nodes_2[b2] == edge_end){T2 = b2;}

        if (c2 != iedge && coordInfoVecs.edges2Nodes_1[c2] == edge_start){T1 = c2;}
        else if (c2 != iedge && coordInfoVecs.edges2Nodes_2[c2] == edge_start){T1 = c2;}
        else if (c2 != iedge && coordInfoVecs.edges2Nodes_1[c2] == edge_end){T2 = c2;}
        else if (c2 != iedge && coordInfoVecs.edges2Nodes_2[c2] == edge_end){T2 = c2;}

        //Now search for the associated 

        int CANDIDATE1_1 = coordInfoVecs.triangles2Nodes_1[H0];
        int CANDIDATE1_2 = coordInfoVecs.triangles2Nodes_2[H0];
        int CANDIDATE1_3 = coordInfoVecs.triangles2Nodes_3[H0];
        int CANDIDATE2_1 = coordInfoVecs.triangles2Nodes_1[T0];
        int CANDIDATE2_2 = coordInfoVecs.triangles2Nodes_2[T0];
        int CANDIDATE2_3 = coordInfoVecs.triangles2Nodes_3[T0];
        
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


        //int temp_edges2Nodes_2 = HEAD;
        //std::cout<<"head tail in loop: "<< HEAD << " "<< TAIL <<std::endl;
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
        if (generalParams.edges_in_upperhem[edges_iteration[0]] == 1){
                    linear_spring_constant = linearSpringInfoVecs.spring_constant_weak;
                }
                else{
                    linear_spring_constant = linearSpringInfoVecs.spring_constant;
        }
        edges_iteration[1] = H1;
        edges_iteration[2] = H2;
        edges_iteration[3] = T1;
        edges_iteration[4] = T2;
            for (int j = 0; j < 5; j++){
                if (generalParams.edges_in_upperhem[edges_iteration[j]] == 1){
                    bend_spring_constant = bendingTriangleInfoVecs.spring_constant_weak;
                }
                else{
                    bend_spring_constant = bendingTriangleInfoVecs.spring_constant;
                }
                int Tri1 = coordInfoVecs.edges2Triangles_1[edges_iteration[j]];//index of the 1st triangle
                int Tri2 = coordInfoVecs.edges2Triangles_2[edges_iteration[j]];
                //int id_k = coordInfoVecs.edges2Nodes_1[edges_iteration[j]];
                //int id_i = coordInfoVecs.edges2Nodes_2[edges_iteration[j]];

                double vec1x, vec1y, vec1z, vec2x, vec2y, vec2z;
                if (Tri1 != Tri2) {
                    int Tri1_n1 = coordInfoVecs.triangles2Nodes_1[Tri1];
                    int Tri1_n2 = coordInfoVecs.triangles2Nodes_2[Tri1];
                    int Tri1_n3 = coordInfoVecs.triangles2Nodes_3[Tri1];
                    vec1x = coordInfoVecs.nodeLocX[Tri1_n2] - coordInfoVecs.nodeLocX[Tri1_n1];
                    vec1y = coordInfoVecs.nodeLocY[Tri1_n2] - coordInfoVecs.nodeLocY[Tri1_n1];
                    vec1z = coordInfoVecs.nodeLocZ[Tri1_n2] - coordInfoVecs.nodeLocZ[Tri1_n1];
                    vec2x = coordInfoVecs.nodeLocX[Tri1_n3] - coordInfoVecs.nodeLocX[Tri1_n1];
                    vec2y = coordInfoVecs.nodeLocY[Tri1_n3] - coordInfoVecs.nodeLocY[Tri1_n1];
                    vec2z = coordInfoVecs.nodeLocZ[Tri1_n3] - coordInfoVecs.nodeLocZ[Tri1_n1];
                    std::vector<double> N1(3);
                    N1[0] = vec1y*vec2z - vec2y*vec1z;
                    N1[1] = -(vec1x*vec2z - vec2x*vec1z);
                    N1[2] = vec1x*vec2y - vec2x*vec1y;
                    double nN1 = sqrt(pow(N1[0],2)+pow(N1[1],2)+pow(N1[2],2));
					//std::cout<<"nN1 = "<<nN1<<std::endl;

                    int Tri2_n1 = coordInfoVecs.triangles2Nodes_1[Tri2];
                    int Tri2_n2 = coordInfoVecs.triangles2Nodes_2[Tri2];
                    int Tri2_n3 = coordInfoVecs.triangles2Nodes_3[Tri2];
                    vec1x = coordInfoVecs.nodeLocX[Tri2_n2] - coordInfoVecs.nodeLocX[Tri2_n1];
                    vec1y = coordInfoVecs.nodeLocY[Tri2_n2] - coordInfoVecs.nodeLocY[Tri2_n1];
                    vec1z = coordInfoVecs.nodeLocZ[Tri2_n2] - coordInfoVecs.nodeLocZ[Tri2_n1];
                    vec2x = coordInfoVecs.nodeLocX[Tri2_n3] - coordInfoVecs.nodeLocX[Tri2_n1];
                    vec2y = coordInfoVecs.nodeLocY[Tri2_n3] - coordInfoVecs.nodeLocY[Tri2_n1];
                    vec2z = coordInfoVecs.nodeLocZ[Tri2_n3] - coordInfoVecs.nodeLocZ[Tri2_n1];
                    std::vector<double> N2(3);
                    N2[0] = vec1y*vec2z - vec2y*vec1z;
                    N2[1] = -(vec1x*vec2z - vec2x*vec1z);
                    N2[2] = vec1x*vec2y - vec2x*vec1y; 
                    double nN2 = sqrt(pow(N2[0],2)+pow(N2[1],2)+pow(N2[2],2));;
					//std::cout<<"nN2 = "<<nN2<<std::endl;

                    double cosAngle = (N1[0]*N2[0] + N1[1]*N2[1] + N1[2]*N2[2])/(nN1*nN2);
					//std::cout<<"dotproduct = "<<cosAngle<<std::endl;
                   
    
                    if (cosAngle > 1.0) {
                        cosAngle = 1.0;
                    }
                    else if (cosAngle < -1.0){
                        cosAngle = -1.0;
                    }

                    double theta_current = acos( cosAngle );
					
                    
                    double local_energy = bend_spring_constant * (1 - cos(theta_current - bendingTriangleInfoVecs.initial_angle) );
                    
                    //bendingTriangleInfoVecs.spring_constant * (1 - cos(theta_current - bendingTriangleInfoVecs.initial_angle) );
                    temp_bend = temp_bend + local_energy;
					/*std::cout<<"bending energy "<<local_energy<<std::endl;
					for (int COUNT = 0; COUNT < 3; COUNT++){
					std::cout<<"unit normal 1 = "<<N1[COUNT]<<std::endl;
					std::cout<<"unit normal 2 = "<<N2[COUNT]<<std::endl;}
					std::cout<<"angle "<<theta_current<<std::endl;*/
                }
            }
               

        double bend_0 = temp_bend;
        //
        double linear_0;
        double DISTANCE = sqrt(
            pow(coordInfoVecs.nodeLocX[edge_end] - coordInfoVecs.nodeLocX[edge_start], 2.0) + 
            pow(coordInfoVecs.nodeLocY[edge_end] - coordInfoVecs.nodeLocY[edge_start], 2.0) + 
            pow(coordInfoVecs.nodeLocZ[edge_end] - coordInfoVecs.nodeLocZ[edge_start], 2.0));
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
            linear_0 = (linear_spring_constant/2)*(DISTANCE - generalParams.Rmin)*
			(DISTANCE - generalParams.Rmin);
        //}
        
        //else if (DISTANCE < generalParams.Rmin ){
        //    linear_0 = linearSpringInfoVecs.spring_constant_rep1*(1-exp(-linearSpringInfoVecs.spring_constant_rep2*(DISTANCE - generalParams.Rmin)))*(1-exp(-linearSpringInfoVecs.spring_constant_rep2*(DISTANCE - generalParams.Rmin)));
        //}
        
        /*double linear_0 = (linearSpringInfoVecs.spring_constant/2)*(sqrt(
            pow(coordInfoVecs.nodeLocX[edge_end] - coordInfoVecs.nodeLocX[edge_start], 2.0) + 
            pow(coordInfoVecs.nodeLocY[edge_end] - coordInfoVecs.nodeLocY[edge_start], 2.0) + 
            pow(coordInfoVecs.nodeLocZ[edge_end] - coordInfoVecs.nodeLocZ[edge_start], 2.0)) - linearSpringInfoVecs.edge_initial_length[0])*
			(sqrt(
            pow(coordInfoVecs.nodeLocX[edge_end] - coordInfoVecs.nodeLocX[edge_start], 2.0) + 
            pow(coordInfoVecs.nodeLocY[edge_end] - coordInfoVecs.nodeLocY[edge_start], 2.0) + 
            pow(coordInfoVecs.nodeLocZ[edge_end] - coordInfoVecs.nodeLocZ[edge_start], 2.0)) - linearSpringInfoVecs.edge_initial_length[0]);*/
			//std::cout<<"the energy of this edge is = "<<linear_0<<std::endl;
        
        int H0n1 = edge_start;//coordInfoVecs.triangles2Nodes_1[H0];
        int H0n2 = edge_end;//coordInfoVecs.triangles2Nodes_2[H0];
        int H0n3 = HEAD;//coordInfoVecs.triangles2Nodes_3[H0];
        int T0n1 = edge_start;//coordInfoVecs.triangles2Nodes_1[T0];
        int T0n2 = TAIL;//coordInfoVecs.triangles2Nodes_2[T0];
        int T0n3 = edge_end;//coordInfoVecs.triangles2Nodes_3[T0];
        double a = sqrt(pow((coordInfoVecs.nodeLocX[H0n2] - coordInfoVecs.nodeLocX[H0n1]),2.0) + 
                    pow((coordInfoVecs.nodeLocY[H0n2] - coordInfoVecs.nodeLocY[H0n1]),2.0) +
                    pow((coordInfoVecs.nodeLocZ[H0n2] - coordInfoVecs.nodeLocZ[H0n1]),2.0)
                    );
        double b = sqrt(pow((coordInfoVecs.nodeLocX[H0n3] - coordInfoVecs.nodeLocX[H0n1]),2.0) + 
                    pow((coordInfoVecs.nodeLocY[H0n3] - coordInfoVecs.nodeLocY[H0n1]),2.0) +
                    pow((coordInfoVecs.nodeLocZ[H0n3] - coordInfoVecs.nodeLocZ[H0n1]),2.0)
                    );
        double c = sqrt(pow((coordInfoVecs.nodeLocX[H0n3] - coordInfoVecs.nodeLocX[H0n2]),2.0) + 
                    pow((coordInfoVecs.nodeLocY[H0n3] - coordInfoVecs.nodeLocY[H0n2]),2.0) +
                    pow((coordInfoVecs.nodeLocZ[H0n3] - coordInfoVecs.nodeLocZ[H0n2]),2.0)
                    );
        double mean_abc = (a + b + c)/2;
        double d = sqrt(pow((coordInfoVecs.nodeLocX[T0n2] - coordInfoVecs.nodeLocX[T0n1]),2.0) + 
                    pow((coordInfoVecs.nodeLocY[T0n2] - coordInfoVecs.nodeLocY[T0n1]),2.0) +
                    pow((coordInfoVecs.nodeLocZ[T0n2] - coordInfoVecs.nodeLocZ[T0n1]),2.0)
                    );
        double e = sqrt(pow((coordInfoVecs.nodeLocX[T0n3] - coordInfoVecs.nodeLocX[T0n1]),2.0) + 
                    pow((coordInfoVecs.nodeLocY[T0n3] - coordInfoVecs.nodeLocY[T0n1]),2.0) +
                    pow((coordInfoVecs.nodeLocZ[T0n3] - coordInfoVecs.nodeLocZ[T0n1]),2.0)
                    );
        double f = sqrt(pow((coordInfoVecs.nodeLocX[T0n3] - coordInfoVecs.nodeLocX[T0n2]),2.0) + 
                    pow((coordInfoVecs.nodeLocY[T0n3] - coordInfoVecs.nodeLocY[T0n2]),2.0) +
                    pow((coordInfoVecs.nodeLocZ[T0n3] - coordInfoVecs.nodeLocZ[T0n2]),2.0)
                    );
        double mean_def = (d + e + f)/2;
        double area_H0 = sqrt(mean_abc*(mean_abc - a)*(mean_abc - b)*(mean_abc - c));
        double area_T0 = sqrt(mean_def*(mean_def - d)*(mean_def - e)*(mean_def - f));
        double area_0_energy = areaTriangleInfoVecs.spring_constant*pow((area_H0 - areaTriangleInfoVecs.initial_area),2.0)/(2*areaTriangleInfoVecs.initial_area) +
                            areaTriangleInfoVecs.spring_constant*pow((area_T0 - areaTriangleInfoVecs.initial_area),2.0)/(2*areaTriangleInfoVecs.initial_area);

        double E_0 = linear_0 + bend_0 + area_0_energy;
        /*std::cout<<"old linear energy: "<<linear_0<<std::endl;
        std::cout<<"old bend energy: "<<bend_0<<std::endl;
        std::cout<<"old area energy: "<<area_0_energy<<std::endl;
        std::cout<<"old total energy: "<<E_0<<std::endl;*/

        
        //Flip the edge, build the data structure for the smaller system.
        bool BAD_CHOICE = false;
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
        }


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
            int H1t1 = coordInfoVecs.edges2Triangles_1[H1];
            int H1t2 = coordInfoVecs.edges2Triangles_2[H1]; //These are the associated triangles to edge H1
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
            
            if (H1t1 == H0){H1t1_vec1x = coordInfoVecs.nodeLocX[edge_start] - coordInfoVecs.nodeLocX[HEAD];
                            H1t1_vec1y = coordInfoVecs.nodeLocY[edge_start] - coordInfoVecs.nodeLocY[HEAD];
                            H1t1_vec1z = coordInfoVecs.nodeLocZ[edge_start] - coordInfoVecs.nodeLocZ[HEAD];
                            H1t1_vec2x = coordInfoVecs.nodeLocX[TAIL] - coordInfoVecs.nodeLocX[HEAD];
                            H1t1_vec2y = coordInfoVecs.nodeLocY[TAIL] - coordInfoVecs.nodeLocY[HEAD];
                            H1t1_vec2z = coordInfoVecs.nodeLocZ[TAIL] - coordInfoVecs.nodeLocZ[HEAD];
                            H1t2_vec1x = coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_2[H1t2]] - coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_1[H1t2]];
                            H1t2_vec1y = coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_2[H1t2]] - coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_1[H1t2]];
                            H1t2_vec1z = coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_2[H1t2]] - coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_1[H1t2]];
                            H1t2_vec2x = coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_3[H1t2]] - coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_1[H1t2]];
                            H1t2_vec2y = coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_3[H1t2]] - coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_1[H1t2]];
                            H1t2_vec2z = coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_3[H1t2]] - coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_1[H1t2]];
							//std::cout<<"H1t1 = H0"<<std::endl;
                            }
            else if (H1t2 == H0){H1t2_vec1x = coordInfoVecs.nodeLocX[edge_start] - coordInfoVecs.nodeLocX[HEAD];
                            H1t2_vec1y = coordInfoVecs.nodeLocY[edge_start] - coordInfoVecs.nodeLocY[HEAD];
                            H1t2_vec1z = coordInfoVecs.nodeLocZ[edge_start] - coordInfoVecs.nodeLocZ[HEAD];
                            H1t2_vec2x = coordInfoVecs.nodeLocX[TAIL] - coordInfoVecs.nodeLocX[HEAD];
                            H1t2_vec2y = coordInfoVecs.nodeLocY[TAIL] - coordInfoVecs.nodeLocY[HEAD];
                            H1t2_vec2z = coordInfoVecs.nodeLocZ[TAIL] - coordInfoVecs.nodeLocZ[HEAD];
                            H1t1_vec1x = coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_2[H1t1]] - coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_1[H1t1]];
                            H1t1_vec1y = coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_2[H1t1]] - coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_1[H1t1]];
                            H1t1_vec1z = coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_2[H1t1]] - coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_1[H1t1]];
                            H1t1_vec2x = coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_3[H1t1]] - coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_1[H1t1]];
                            H1t1_vec2y = coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_3[H1t1]] - coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_1[H1t1]];
                            H1t1_vec2z = coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_3[H1t1]] - coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_1[H1t1]];
							//std::cout<<"H1t2 = H0"<<std::endl;
                            }
            int H2t1 = coordInfoVecs.edges2Triangles_1[H2];
            int H2t2 = coordInfoVecs.edges2Triangles_2[H2]; //These are the associated triangles to edge H2        
            if (H2t1 == H0){//In this case H2t1 turns into T0.
                            H2t1_vec1x = coordInfoVecs.nodeLocX[TAIL] - coordInfoVecs.nodeLocX[HEAD];
                            H2t1_vec1y = coordInfoVecs.nodeLocY[TAIL] - coordInfoVecs.nodeLocY[HEAD];
                            H2t1_vec1z = coordInfoVecs.nodeLocZ[TAIL] - coordInfoVecs.nodeLocZ[HEAD];
                            H2t1_vec2x = coordInfoVecs.nodeLocX[edge_end] - coordInfoVecs.nodeLocX[HEAD];
                            H2t1_vec2y = coordInfoVecs.nodeLocY[edge_end] - coordInfoVecs.nodeLocY[HEAD];
                            H2t1_vec2z = coordInfoVecs.nodeLocZ[edge_end] - coordInfoVecs.nodeLocZ[HEAD];
                            H2t2_vec1x = coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_2[H2t2]] - coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_1[H2t2]];
                            H2t2_vec1y = coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_2[H2t2]] - coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_1[H2t2]];
                            H2t2_vec1z = coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_2[H2t2]] - coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_1[H2t2]];
                            H2t2_vec2x = coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_3[H2t2]] - coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_1[H2t2]];
                            H2t2_vec2y = coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_3[H2t2]] - coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_1[H2t2]];
                            H2t2_vec2z = coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_3[H2t2]] - coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_1[H2t2]];
							//std::cout<<"H2t1 = H0"<<std::endl;
                            }
            else if (H2t2 == H0){//In this case H2t2 tunrs into T0
                            H2t2_vec1x = coordInfoVecs.nodeLocX[TAIL] - coordInfoVecs.nodeLocX[HEAD];
                            H2t2_vec1y = coordInfoVecs.nodeLocY[TAIL] - coordInfoVecs.nodeLocY[HEAD];
                            H2t2_vec1z = coordInfoVecs.nodeLocZ[TAIL] - coordInfoVecs.nodeLocZ[HEAD];
                            H2t2_vec2x = coordInfoVecs.nodeLocX[edge_end] - coordInfoVecs.nodeLocX[HEAD];
                            H2t2_vec2y = coordInfoVecs.nodeLocY[edge_end] - coordInfoVecs.nodeLocY[HEAD];
                            H2t2_vec2z = coordInfoVecs.nodeLocZ[edge_end] - coordInfoVecs.nodeLocZ[HEAD];
                            H2t1_vec1x = coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_2[H2t1]] - coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_1[H2t1]];
                            H2t1_vec1y = coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_2[H2t1]] - coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_1[H2t1]];
                            H2t1_vec1z = coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_2[H2t1]] - coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_1[H2t1]];
                            H2t1_vec2x = coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_3[H2t1]] - coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_1[H2t1]];
                            H2t1_vec2y = coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_3[H2t1]] - coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_1[H2t1]];
                            H2t1_vec2z = coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_3[H2t1]] - coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_1[H2t1]];
							//std::cout<<"H2t2 = H0"<<std::endl;
                            }
            int T1t1 = coordInfoVecs.edges2Triangles_1[T1];
            int T1t2 = coordInfoVecs.edges2Triangles_2[T1];
            if (T1t1 == T0){//In this case T1t1 turns into H0.
                            T1t1_vec1x = coordInfoVecs.nodeLocX[HEAD] - coordInfoVecs.nodeLocX[TAIL];
                            T1t1_vec1y = coordInfoVecs.nodeLocY[HEAD] - coordInfoVecs.nodeLocY[TAIL];
                            T1t1_vec1z = coordInfoVecs.nodeLocZ[HEAD] - coordInfoVecs.nodeLocZ[TAIL];
                            T1t1_vec2x = coordInfoVecs.nodeLocX[edge_start] - coordInfoVecs.nodeLocX[TAIL];
                            T1t1_vec2y = coordInfoVecs.nodeLocY[edge_start] - coordInfoVecs.nodeLocY[TAIL];
                            T1t1_vec2z = coordInfoVecs.nodeLocZ[edge_start] - coordInfoVecs.nodeLocZ[TAIL];
                            T1t2_vec1x = coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_2[T1t2]] - coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_1[T1t2]];
                            T1t2_vec1y = coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_2[T1t2]] - coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_1[T1t2]];
                            T1t2_vec1z = coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_2[T1t2]] - coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_1[T1t2]];
                            T1t2_vec2x = coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_3[T1t2]] - coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_1[T1t2]];
                            T1t2_vec2y = coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_3[T1t2]] - coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_1[T1t2]];
                            T1t2_vec2z = coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_3[T1t2]] - coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_1[T1t2]];
							//std::cout<<"T1t1 = T0"<<std::endl;
                            }
            else if (T1t2 == T0){//In this case T1t2 turns into H0.
                            T1t2_vec1x = coordInfoVecs.nodeLocX[HEAD] - coordInfoVecs.nodeLocX[TAIL];
                            T1t2_vec1y = coordInfoVecs.nodeLocY[HEAD] - coordInfoVecs.nodeLocY[TAIL];
                            T1t2_vec1z = coordInfoVecs.nodeLocZ[HEAD] - coordInfoVecs.nodeLocZ[TAIL];
                            T1t2_vec2x = coordInfoVecs.nodeLocX[edge_start] - coordInfoVecs.nodeLocX[TAIL];
                            T1t2_vec2y = coordInfoVecs.nodeLocY[edge_start] - coordInfoVecs.nodeLocY[TAIL];
                            T1t2_vec2z = coordInfoVecs.nodeLocZ[edge_start] - coordInfoVecs.nodeLocZ[TAIL];
                            T1t1_vec1x = coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_2[T1t1]] - coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_1[T1t1]];
                            T1t1_vec1y = coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_2[T1t1]] - coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_1[T1t1]];
                            T1t1_vec1z = coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_2[T1t1]] - coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_1[T1t1]];
                            T1t1_vec2x = coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_3[T1t1]] - coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_1[T1t1]];
                            T1t1_vec2y = coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_3[T1t1]] - coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_1[T1t1]];
                            T1t1_vec2z = coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_3[T1t1]] - coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_1[T1t1]];
							//std::cout<<"T1t2 = T0"<<std::endl;
                            }
            int T2t1 = coordInfoVecs.edges2Triangles_1[T2];
            int T2t2 = coordInfoVecs.edges2Triangles_2[T2];
            if (T2t1 == T0){T2t1_vec1x = coordInfoVecs.nodeLocX[edge_end] - coordInfoVecs.nodeLocX[TAIL];
                            T2t1_vec1y = coordInfoVecs.nodeLocY[edge_end] - coordInfoVecs.nodeLocY[TAIL];
                            T2t1_vec1z = coordInfoVecs.nodeLocZ[edge_end] - coordInfoVecs.nodeLocZ[TAIL];
                            T2t1_vec2x = coordInfoVecs.nodeLocX[HEAD] - coordInfoVecs.nodeLocX[TAIL];
                            T2t1_vec2y = coordInfoVecs.nodeLocY[HEAD] - coordInfoVecs.nodeLocY[TAIL];
                            T2t1_vec2z = coordInfoVecs.nodeLocZ[HEAD] - coordInfoVecs.nodeLocZ[TAIL];
                            T2t2_vec1x = coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_2[T2t2]] - coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_1[T2t2]];
                            T2t2_vec1y = coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_2[T2t2]] - coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_1[T2t2]];
                            T2t2_vec1z = coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_2[T2t2]] - coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_1[T2t2]];
                            T2t2_vec2x = coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_3[T2t2]] - coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_1[T2t2]];
                            T2t2_vec2y = coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_3[T2t2]] - coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_1[T2t2]];
                            T2t2_vec2z = coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_3[T2t2]] - coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_1[T2t2]];
							//std::cout<<"T2t1 = T0"<<std::endl;
                            }
            else if (T2t2 == T0){
                            T2t2_vec1x = coordInfoVecs.nodeLocX[edge_end] - coordInfoVecs.nodeLocX[TAIL];
                            T2t2_vec1y = coordInfoVecs.nodeLocY[edge_end] - coordInfoVecs.nodeLocY[TAIL];
                            T2t2_vec1z = coordInfoVecs.nodeLocZ[edge_end] - coordInfoVecs.nodeLocZ[TAIL];
                            T2t2_vec2x = coordInfoVecs.nodeLocX[HEAD] - coordInfoVecs.nodeLocX[TAIL];
                            T2t2_vec2y = coordInfoVecs.nodeLocY[HEAD] - coordInfoVecs.nodeLocY[TAIL];
                            T2t2_vec2z = coordInfoVecs.nodeLocZ[HEAD] - coordInfoVecs.nodeLocZ[TAIL];
                            T2t1_vec1x = coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_2[T2t1]] - coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_1[T2t1]];
                            T2t1_vec1y = coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_2[T2t1]] - coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_1[T2t1]];
                            T2t1_vec1z = coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_2[T2t1]] - coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_1[T2t1]];
                            T2t1_vec2x = coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_3[T2t1]] - coordInfoVecs.nodeLocX[coordInfoVecs.triangles2Nodes_1[T2t1]];
                            T2t1_vec2y = coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_3[T2t1]] - coordInfoVecs.nodeLocY[coordInfoVecs.triangles2Nodes_1[T2t1]];
                            T2t1_vec2z = coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_3[T2t1]] - coordInfoVecs.nodeLocZ[coordInfoVecs.triangles2Nodes_1[T2t1]];
							//std::cout<<"T2t2 = T0"<<std::endl;
                            }
            
            //First calculate the linear spring energy due to edge-swap.
        double linear_1;
        double DISTANCE = sqrt(
           pow(coordInfoVecs.nodeLocX[HEAD] - coordInfoVecs.nodeLocX[TAIL], 2.0) + 
           pow(coordInfoVecs.nodeLocY[HEAD] - coordInfoVecs.nodeLocY[TAIL], 2.0) + 
           pow(coordInfoVecs.nodeLocZ[HEAD] - coordInfoVecs.nodeLocZ[TAIL], 2.0));
        
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
            linear_1 = (linearSpringInfoVecs.spring_constant/2)*(DISTANCE - generalParams.Rmin)*
			(DISTANCE - generalParams.Rmin);
        //}
        
        //else if (DISTANCE < generalParams.Rmin){
         //   linear_1 =   linearSpringInfoVecs.spring_constant_rep1*(1-exp(-linearSpringInfoVecs.spring_constant_rep2*(DISTANCE - generalParams.Rmin)))*(1-exp(-linearSpringInfoVecs.spring_constant_rep2*(DISTANCE - generalParams.Rmin)));
        //}
            
        double prob;
        double random_number;
        double Edif;
        if (DISTANCE >= 0.0){
            //WARNING: RESET BENDING COUNTER
            temp_bend = 0.0;


    
            double N1_vec1x, N1_vec1y, N1_vec1z, N1_vec2x, N1_vec2y, N1_vec2z, N2_vec1x, N2_vec1y, N2_vec1z, N2_vec2x, N2_vec2y, N2_vec2z;
            bool THIS_SHOULD_NOT_HAPPEN = false;
            for (int j = 0; j < 5; j++){
                if (generalParams.edges_in_upperhem[edges_iteration[j]] == 1){
                    bend_spring_constant = bendingTriangleInfoVecs.spring_constant_weak;
                }
                else{
                    bend_spring_constant = bendingTriangleInfoVecs.spring_constant;
                }
                    if (j == 0){
                        N1_vec1x = coordInfoVecs.nodeLocX[HEAD] - coordInfoVecs.nodeLocX[TAIL];//x component of the 1st vector to calculate N1
                        N1_vec1y = coordInfoVecs.nodeLocY[HEAD] - coordInfoVecs.nodeLocY[TAIL];
                        N1_vec1z = coordInfoVecs.nodeLocZ[HEAD] - coordInfoVecs.nodeLocZ[TAIL];
                        N1_vec2x = coordInfoVecs.nodeLocX[edge_start] - coordInfoVecs.nodeLocX[TAIL];//x component of the 2nd vector to calculate N1
                        N1_vec2y = coordInfoVecs.nodeLocY[edge_start] - coordInfoVecs.nodeLocY[TAIL];
                        N1_vec2z = coordInfoVecs.nodeLocZ[edge_start] - coordInfoVecs.nodeLocZ[TAIL];
                        N2_vec1x = coordInfoVecs.nodeLocX[edge_end] - coordInfoVecs.nodeLocX[TAIL];
                        N2_vec1y = coordInfoVecs.nodeLocY[edge_end] - coordInfoVecs.nodeLocY[TAIL];
                        N2_vec1z = coordInfoVecs.nodeLocZ[edge_end] - coordInfoVecs.nodeLocZ[TAIL];
                        N2_vec2x = coordInfoVecs.nodeLocX[HEAD] - coordInfoVecs.nodeLocX[TAIL];
                        N2_vec2y = coordInfoVecs.nodeLocY[HEAD] - coordInfoVecs.nodeLocY[TAIL];
                        N2_vec2z = coordInfoVecs.nodeLocZ[HEAD] - coordInfoVecs.nodeLocZ[TAIL];
                    }
                    else if (j == 1){
                        N1_vec1x = H1t1_vec1x;
                        N1_vec1y = H1t1_vec1y;
                        N1_vec1z = H1t1_vec1z;
                        N1_vec2x = H1t1_vec2x;
                        N1_vec2y = H1t1_vec2y;
                        N1_vec2z = H1t1_vec2z;
                        N2_vec1x = H1t2_vec1x;
                        N2_vec1y = H1t2_vec1y;
                        N2_vec1z = H1t2_vec1z;
                        N2_vec2x = H1t2_vec2x;
                        N2_vec2y = H1t2_vec2y;
                        N2_vec2z = H1t2_vec2z;
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
                    std::vector<double> N1(3);
                    N1[0] = N1_vec1y*N1_vec2z - N1_vec2y*N1_vec1z;
                    N1[1] = -(N1_vec1x*N1_vec2z - N1_vec2x*N1_vec1z);
                    N1[2] = N1_vec1x*N1_vec2y - N1_vec2x*N1_vec1y;
                    double nN1 = sqrt(pow(N1[0],2)+pow(N1[1],2)+pow(N1[2],2));
					//std::cout<<"newnN1 = "<<nN1<<std::endl;

    
                    std::vector<double> N2(3);
                    N2[0] = N2_vec1y*N2_vec2z - N2_vec2y*N2_vec1z;
                    N2[1] = -(N2_vec1x*N2_vec2z - N2_vec2x*N2_vec1z);
                    N2[2] = N2_vec1x*N2_vec2y - N2_vec2x*N2_vec1y; 
                    double nN2 = sqrt(pow(N2[0],2)+pow(N2[1],2)+pow(N2[2],2));;
					//std::cout<<"newnN2 = "<<nN2<<std::endl;

                    double cosAngle = (N1[0]*N2[0] + N1[1]*N2[1] + N1[2]*N2[2])/(nN1*nN2);
					//std::cout<<"cosAngle = "<<cosAngle<<std::endl;
                    
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
                    
                    double local_energy = bend_spring_constant * (1 - cos(theta_current - bendingTriangleInfoVecs.initial_angle) );
                    temp_bend = temp_bend + local_energy;
					//std::cout<<"bending energy "<<local_energy<<std::endl;
					/*for (int COUNT = 0; COUNT < 3; COUNT++){
					std::cout<<"unit normal 1 = "<<N1[COUNT]<<std::endl;
					std::cout<<"unit normal 2 = "<<N2[COUNT]<<std::endl;}
					std::cout<<"angle "<<theta_current<<std::endl;*/
                }
            double bend_1 = temp_bend;

            int H0n1 = HEAD;
            int H0n2 = edge_start;
            int H0n3 = TAIL;
            int T0n1 = TAIL;
            int T0n2 = edge_end;
            int T0n3 = HEAD;
            double a = sqrt(pow((coordInfoVecs.nodeLocX[H0n2] - coordInfoVecs.nodeLocX[H0n1]),2.0) + 
                    pow((coordInfoVecs.nodeLocY[H0n2] - coordInfoVecs.nodeLocY[H0n1]),2.0) +
                    pow((coordInfoVecs.nodeLocZ[H0n2] - coordInfoVecs.nodeLocZ[H0n1]),2.0)
                    );
            double b = sqrt(pow((coordInfoVecs.nodeLocX[H0n3] - coordInfoVecs.nodeLocX[H0n1]),2.0) + 
                        pow((coordInfoVecs.nodeLocY[H0n3] - coordInfoVecs.nodeLocY[H0n1]),2.0) +
                        pow((coordInfoVecs.nodeLocZ[H0n3] - coordInfoVecs.nodeLocZ[H0n1]),2.0)
                        );
            double c = sqrt(pow((coordInfoVecs.nodeLocX[H0n3] - coordInfoVecs.nodeLocX[H0n2]),2.0) + 
                        pow((coordInfoVecs.nodeLocY[H0n3] - coordInfoVecs.nodeLocY[H0n2]),2.0) +
                        pow((coordInfoVecs.nodeLocZ[H0n3] - coordInfoVecs.nodeLocZ[H0n2]),2.0)
                        );
            double mean_abc = (a + b + c)/2;
            double d = sqrt(pow((coordInfoVecs.nodeLocX[T0n2] - coordInfoVecs.nodeLocX[T0n1]),2.0) + 
                        pow((coordInfoVecs.nodeLocY[T0n2] - coordInfoVecs.nodeLocY[T0n1]),2.0) +
                        pow((coordInfoVecs.nodeLocZ[T0n2] - coordInfoVecs.nodeLocZ[T0n1]),2.0)
                        );
            double e = sqrt(pow((coordInfoVecs.nodeLocX[T0n3] - coordInfoVecs.nodeLocX[T0n1]),2.0) + 
                        pow((coordInfoVecs.nodeLocY[T0n3] - coordInfoVecs.nodeLocY[T0n1]),2.0) +
                        pow((coordInfoVecs.nodeLocZ[T0n3] - coordInfoVecs.nodeLocZ[T0n1]),2.0)
                        );
            double f = sqrt(pow((coordInfoVecs.nodeLocX[T0n3] - coordInfoVecs.nodeLocX[T0n2]),2.0) + 
                        pow((coordInfoVecs.nodeLocY[T0n3] - coordInfoVecs.nodeLocY[T0n2]),2.0) +
                        pow((coordInfoVecs.nodeLocZ[T0n3] - coordInfoVecs.nodeLocZ[T0n2]),2.0)
                        );
            double mean_def = (d + e + f)/2;
            double area_H0 = sqrt(mean_abc*(mean_abc - a)*(mean_abc - b)*(mean_abc - c));
            double area_T0 = sqrt(mean_def*(mean_def - d)*(mean_def - e)*(mean_def - f));
            double area_1_energy = linear_spring_constant*pow((area_H0 - areaTriangleInfoVecs.initial_area),2.0)/(2*areaTriangleInfoVecs.initial_area) +
                                areaTriangleInfoVecs.spring_constant*pow((area_T0 - areaTriangleInfoVecs.initial_area),2.0)/(2*areaTriangleInfoVecs.initial_area);
            double E_1 = linear_1 + bend_1 + area_1_energy;
			/*std::cout<<"new linear energy = "<<linear_1<<std::endl;
			std::cout<<"new bend energy = "<<bend_1<<std::endl;
			std::cout<<"new area energy = "<<area_1_energy<<std::endl;
			std::cout<<"new total energy: "<<E_1<<std::endl;*/
        
        //Now compute the Boltzmann factor to determine if a swap occurs.
        std::random_device rd;  //Will be used to obtain a seed for the random number engine
        std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
        std::uniform_real_distribution<> dis(0.0, 1.0);
        random_number = dis(gen);
        //double random_number = 0;
        Edif = (E_1 - E_0);
        

         prob = generalParams.tau*exp(-Edif/generalParams.kT);
		//std::cout<<"P(swap): "<<prob<<std::endl;
        }
        else{
            prob = -1.0;
        }

        bool ACCEPT2;
        if (!isnan(prob)){
            if (prob >= 1){ACCEPT2 = true;}
            else if (prob < 1 && random_number <= prob){ACCEPT2 = true;}
            else if (prob < 1 && random_number > prob){ACCEPT2 = false;}
        }
        else{ACCEPT2 = false;}
        //std::cout<<"ACCEPT2 = "<<ACCEPT2<<std::endl;
        //Perform real update
        //if (ACCEPT2 == true){
        //if (Edif < 100.0){
        if (ACCEPT2 == true ){//&& THIS_SHOULD_NOT_HAPPEN == false){
            alpha = 1;
            coordInfoVecs.triangles2Nodes_1[H0] = HEAD;
            coordInfoVecs.triangles2Nodes_2[H0] = edge_start;
            coordInfoVecs.triangles2Nodes_3[H0] = TAIL;
            coordInfoVecs.triangles2Nodes_1[T0] = HEAD;
            coordInfoVecs.triangles2Nodes_2[T0] = TAIL;
            coordInfoVecs.triangles2Nodes_3[T0] = edge_end;
            coordInfoVecs.triangles2Edges_1[H0] = iedge;
            coordInfoVecs.triangles2Edges_2[H0] = H1;
            coordInfoVecs.triangles2Edges_3[H0] = T1;
            coordInfoVecs.triangles2Edges_1[T0] = iedge;
            coordInfoVecs.triangles2Edges_2[T0] = T2;
            coordInfoVecs.triangles2Edges_3[T0] = H2;
            coordInfoVecs.edges2Nodes_1[iedge] = TAIL;
            coordInfoVecs.edges2Nodes_2[iedge] = HEAD;
            H1t1 = coordInfoVecs.edges2Triangles_1[H1];
            H1t2 = coordInfoVecs.edges2Triangles_2[H1]; //These are the associated triangles to edge H1
            if (H1t1 == H0){coordInfoVecs.edges2Triangles_1[H1] = H0;}
            if (H1t2 == H0){coordInfoVecs.edges2Triangles_2[H1] = H0;}
            H2t1 = coordInfoVecs.edges2Triangles_1[H2];
            H2t2 = coordInfoVecs.edges2Triangles_2[H2]; //These are the associated triangles to edge H2        
            if (H2t1 == H0){coordInfoVecs.edges2Triangles_1[H2] = T0;}
            if (H2t2 == H0){coordInfoVecs.edges2Triangles_2[H2] = T0;}
            T1t1 = coordInfoVecs.edges2Triangles_1[T1];
            T1t2 = coordInfoVecs.edges2Triangles_2[T1];
            if (T1t1 == T0){coordInfoVecs.edges2Triangles_1[T1] = H0;}
            if (T1t2 == T0){coordInfoVecs.edges2Triangles_2[T1] = H0;}
            T2t1 = coordInfoVecs.edges2Triangles_1[T2];
            T2t2 = coordInfoVecs.edges2Triangles_2[T2];
            if (T2t1 == T0){coordInfoVecs.edges2Triangles_1[T2] = T0;}
            if (T2t2 == T0){coordInfoVecs.edges2Triangles_2[T2] = T0;}

            ////////////////////////////////////////////////////////////////////////////////////
            ///////////// UPDATING NEIGHBORING NODE INFO ///////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////

            ///////// DELETING CONNECTIVITY BETWEEN EDGE_START AND EDGE_END ////////////////////
            int data_id;
            if (coordInfoVecs.nndata1[edge_start] == edge_end){
                //data_id = 0;
                coordInfoVecs.nndata1[edge_start] = -2;
            }
            else if (coordInfoVecs.nndata2[edge_start] == edge_end){
                //data_id = 1;
                coordInfoVecs.nndata2[edge_start] = -2;
            }
            else if (coordInfoVecs.nndata3[edge_start] == edge_end){
            //    data_id = 2;
                coordInfoVecs.nndata3[edge_start] = -2;
            }
            else if (coordInfoVecs.nndata4[edge_start] == edge_end){
              //  data_id = 3;
                coordInfoVecs.nndata4[edge_start] = -2;
            }
            else if (coordInfoVecs.nndata5[edge_start] == edge_end){
                //data_id = 4;
                coordInfoVecs.nndata5[edge_start] = -2;
            }
            else if (coordInfoVecs.nndata6[edge_start] == edge_end){
            //    data_id = 5;
                coordInfoVecs.nndata6[edge_start] = -2;
            }
            else if (coordInfoVecs.nndata7[edge_start] == edge_end){
            //    data_id = 6;
                coordInfoVecs.nndata7[edge_start] = -2;
            }
            else if (coordInfoVecs.nndata8[edge_start] == edge_end){
             //   data_id = 7;
                coordInfoVecs.nndata8[edge_start] = -2;
            }
            else if (coordInfoVecs.nndata9[edge_start] == edge_end){
               // data_id = 8;
                coordInfoVecs.nndata9[edge_start] = -2;
            }
            else if (coordInfoVecs.nndata10[edge_start] == edge_end){
             //   data_id = 9;
                coordInfoVecs.nndata10[edge_start] = -2;
            }
            else if (coordInfoVecs.nndata11[edge_start] == edge_end){
               // data_id = 10;
                coordInfoVecs.nndata11[edge_start] = -2;
            }
            else if (coordInfoVecs.nndata12[edge_start] == edge_end){
             //   data_id = 11;
                coordInfoVecs.nndata12[edge_start] = -2;
            }
            else {}

            if (coordInfoVecs.nndata1[edge_end] == edge_start){
               // data_id = 0;
                coordInfoVecs.nndata1[edge_end] = -2;
            }
            else if (coordInfoVecs.nndata2[edge_end] == edge_start){
              //  data_id = 1;
                coordInfoVecs.nndata2[edge_end] = -2;
            }
            else if (coordInfoVecs.nndata3[edge_end] == edge_start){
              //  data_id = 2;
                coordInfoVecs.nndata3[edge_end] = -2;
            }
            else if (coordInfoVecs.nndata4[edge_end] == edge_start){
              //  data_id = 3;
                coordInfoVecs.nndata4[edge_end] = -2;
            }
            else if (coordInfoVecs.nndata5[edge_end] == edge_start){
              //  data_id = 4;
                coordInfoVecs.nndata5[edge_end] = -2;
            }
            else if (coordInfoVecs.nndata6[edge_end] == edge_start){
              //  data_id = 5;
                coordInfoVecs.nndata6[edge_end] = -2;
            }
            else if (coordInfoVecs.nndata7[edge_end] == edge_start){
              //  data_id = 6;
                coordInfoVecs.nndata7[edge_end] = -2;
            }
            else if (coordInfoVecs.nndata8[edge_end] == edge_start){
               // data_id = 7;
                coordInfoVecs.nndata8[edge_end] = -2;
            }
            else if (coordInfoVecs.nndata9[edge_end] == edge_start){
               // data_id = 8;
                coordInfoVecs.nndata9[edge_end] = -2;
            }
            else if (coordInfoVecs.nndata10[edge_end] == edge_start){
              //  data_id = 9;
                coordInfoVecs.nndata10[edge_end] = -2;
            }
            else if (coordInfoVecs.nndata11[edge_end] == edge_start){
             //   data_id = 10;
                coordInfoVecs.nndata11[edge_end] = -2;
            }
            else if (coordInfoVecs.nndata12[edge_end] == edge_start){
             //   data_id = 11;
                coordInfoVecs.nndata12[edge_end] = -2;
            }
            else {}


           
            ///////////// ESTABLISHING NEW CONNECTIVITY ////////////////////
           if (coordInfoVecs.nndata1[HEAD] < 0){
             //   data_id = 0;
                coordInfoVecs.nndata1[HEAD] = TAIL;
            }
            else if (coordInfoVecs.nndata2[HEAD] < 0){
             //   data_id = 1;
                coordInfoVecs.nndata2[HEAD] = TAIL;
            }
            else if (coordInfoVecs.nndata3[HEAD] < 0){
             //   data_id = 2;
                coordInfoVecs.nndata3[HEAD] = TAIL;
            }
            else if (coordInfoVecs.nndata4[HEAD] < 0){
              //  data_id = 3;
                coordInfoVecs.nndata4[HEAD] = TAIL;
            }
            else if (coordInfoVecs.nndata5[HEAD] < 0){
             //   data_id = 4;
                coordInfoVecs.nndata5[HEAD] = TAIL;
            }
            else if (coordInfoVecs.nndata6[HEAD] < 0){
             //   data_id = 5;
                coordInfoVecs.nndata6[HEAD] = TAIL;
            }
            else if (coordInfoVecs.nndata7[HEAD] < 0){
              //  data_id = 6;
                coordInfoVecs.nndata7[HEAD] = TAIL;
            }
            else if (coordInfoVecs.nndata8[HEAD] < 0){
              //  data_id = 7;
                coordInfoVecs.nndata8[HEAD] = TAIL;
            }
            else if (coordInfoVecs.nndata9[HEAD] < 0){
              //  data_id = 8;
                coordInfoVecs.nndata9[HEAD] = TAIL;
            }
            else if (coordInfoVecs.nndata10[HEAD] < 0){
              //  data_id = 9;
                coordInfoVecs.nndata10[HEAD] = TAIL;
            }
            else if (coordInfoVecs.nndata11[HEAD] < 0){
              //  data_id = 10;
                coordInfoVecs.nndata11[HEAD] = TAIL;
            }
            else if (coordInfoVecs.nndata12[HEAD] < 0){
              //  data_id = 11;
                coordInfoVecs.nndata12[HEAD] = TAIL;
            }
            else {}

            if (coordInfoVecs.nndata1[TAIL] < 0){
              //  data_id = 0;
                coordInfoVecs.nndata1[TAIL] = HEAD;
            }
            else if (coordInfoVecs.nndata2[TAIL] < 0){
              //  data_id = 1;
                coordInfoVecs.nndata2[TAIL] = HEAD;
            }
            else if (coordInfoVecs.nndata3[TAIL] < 0){
              //  data_id = 2;
                coordInfoVecs.nndata3[TAIL] = HEAD;
            }
            else if (coordInfoVecs.nndata4[TAIL] < 0){
              //  data_id = 3;
                coordInfoVecs.nndata4[TAIL] = HEAD;
            }
            else if (coordInfoVecs.nndata5[TAIL] < 0){
              //  data_id = 4;
                coordInfoVecs.nndata5[TAIL] = HEAD;
            }
            else if (coordInfoVecs.nndata6[TAIL] < 0){
              //  data_id = 5;
                coordInfoVecs.nndata6[TAIL] = HEAD;
            }
            else if (coordInfoVecs.nndata7[TAIL] < 0){
              //  data_id = 6;
                coordInfoVecs.nndata7[TAIL] = HEAD;
            }
            else if (coordInfoVecs.nndata8[TAIL] < 0){
              //  data_id = 7;
                coordInfoVecs.nndata8[TAIL] = HEAD;
            }
            else if (coordInfoVecs.nndata9[TAIL] < 0){
              //  data_id = 8;
                coordInfoVecs.nndata9[TAIL] = HEAD;
            }
            else if (coordInfoVecs.nndata10[TAIL] < 0){
              //  data_id = 9;
                coordInfoVecs.nndata10[TAIL] = HEAD;
            }
            else if (coordInfoVecs.nndata11[TAIL] < 0){
              //  data_id = 10;
                coordInfoVecs.nndata11[TAIL] = HEAD;
            }
            else if (coordInfoVecs.nndata12[TAIL] < 0){
              //  data_id = 11;
                coordInfoVecs.nndata12[TAIL] = HEAD;
            }
            else {}

            nndata[HEAD] += 1;
            nndata[TAIL] += 1;
            nndata[edge_start] -= 1;
            nndata[edge_end] -= 1;
            
    
        }
    } 
    
    };  
    return alpha;
//This completes the update (if necessary) of the following data structures: triangles2Nodes, edges2Nodes, edges2Triangles.
};




//copy configuration from device to host
void Edgeswap::transferDtoH(CoordInfoVecs& coordInfoVecs,
    HostSetInfoVecs& hostSetInfoVecs){
    thrust::copy(coordInfoVecs.nodeLocX.begin(),coordInfoVecs.nodeLocX.end(),hostSetInfoVecs.nodeLocX.begin());
    thrust::copy(coordInfoVecs.nodeLocY.begin(),coordInfoVecs.nodeLocY.end(),hostSetInfoVecs.nodeLocY.begin());
    thrust::copy(coordInfoVecs.nodeLocZ.begin(),coordInfoVecs.nodeLocZ.end(),hostSetInfoVecs.nodeLocZ.begin());

    
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

};

//copy configuration from host to device
void Edgeswap::transferHtoD(CoordInfoVecs& coordInfoVecs,
    HostSetInfoVecs& hostSetInfoVecs){
    thrust::copy(hostSetInfoVecs.nodeLocX.begin(),hostSetInfoVecs.nodeLocX.end(),coordInfoVecs.nodeLocX.begin());
    thrust::copy(hostSetInfoVecs.nodeLocY.begin(),hostSetInfoVecs.nodeLocY.end(),coordInfoVecs.nodeLocY.begin());
    thrust::copy(hostSetInfoVecs.nodeLocZ.begin(),hostSetInfoVecs.nodeLocZ.end(),coordInfoVecs.nodeLocZ.begin());

    
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

};
