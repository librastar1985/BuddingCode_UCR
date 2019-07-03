#ifndef LINEARSPRINGS_H_
#define LINEARSPRINGS_H_ 

#include "SystemStructures.h"

void ComputeLinearSprings(
    GeneralParams& generalParams,
    CoordInfoVecs& coordInfoVecs,
    LinearSpringInfoVecs& linearSpringInfoVecs,
    LJInfoVecs& ljInfoVecs);
    
struct LinearSpringFunctor {
    
    double spring_constant;
    double spring_constant_weak;
    int* edges_in_upperhem;
    double* locXAddr;
    double* locYAddr;
    double* locZAddr;

    int* idKey;
    double* forceXAddr;
    double* forceYAddr;
    double* forceZAddr;
    
	__host__ __device__ LinearSpringFunctor(
        double& _spring_constant,
        double& _spring_constant_weak,
        int* _edges_in_upperhem,
        double* _locXAddr,
        double* _locYAddr,
        double* _locZAddr,
        
        int* _idKey,
        double* _forceXAddr,
        double* _forceYAddr,
        double* _forceZAddr) :
        spring_constant(_spring_constant),
        spring_constant_weak(_spring_constant_weak),
        edges_in_upperhem(_edges_in_upperhem),
        locXAddr(_locXAddr),
        locYAddr(_locYAddr),
        locZAddr(_locZAddr),
        idKey(_idKey),
        forceXAddr(_forceXAddr),
        forceYAddr(_forceYAddr),
        forceZAddr(_forceZAddr) {}

	//hand in counting iterator and id of two edges and preferred length
	__device__ double operator()(const Tuuud &u3d) {
        		
        //counter ranges from 0 to num_edges. 
        int counter = thrust::get<0>(u3d);
		int place = 2 * counter;//represents location in write to vector.

        int edgeL = thrust::get<1>(u3d);
        int edgeR = thrust::get<2>(u3d);

        if (edgeL != INT_MAX && edgeR != INT_MAX){
            double what_spring_constant;
            if (edges_in_upperhem[counter] == 1){
                what_spring_constant = spring_constant_weak;
            }
            else if (edges_in_upperhem[counter] == 0){
                what_spring_constant = (spring_constant_weak + spring_constant)/2.0;
            }
            else{
                what_spring_constant = spring_constant;
            }

            
            double length_zero = thrust::get<3>(u3d);

            // compute forces.
            double xLoc_LR = locXAddr[edgeL] - locXAddr[edgeR];
            double yLoc_LR = locYAddr[edgeL] - locYAddr[edgeR];
            double zLoc_LR = locZAddr[edgeL] - locZAddr[edgeR];
    

            double length_current = sqrt( (xLoc_LR) * (xLoc_LR) + 
                                        (yLoc_LR) * (yLoc_LR)  + 
                                        (zLoc_LR) * (zLoc_LR) );

        double energy = 0.0;
        if (length_current != length_zero){
            double magnitude = -what_spring_constant * (length_current - length_zero);
            
            
                idKey[place] = edgeL;

                
            //issue here writing to vectors with force????
                forceXAddr[place] = magnitude * (xLoc_LR/length_current);
                forceYAddr[place] = magnitude * (yLoc_LR/length_current);
                forceZAddr[place] = magnitude * (zLoc_LR/length_current);

                idKey[place + 1] = edgeR;
                forceXAddr[place + 1] = -magnitude * (xLoc_LR/length_current);
                forceYAddr[place + 1] = -magnitude * (yLoc_LR/length_current);
                forceZAddr[place + 1] = -magnitude * (zLoc_LR/length_current);
            
            energy = (what_spring_constant/2.0) * (length_current - length_zero) * (length_current - length_zero);
        }
            return energy;
        }

        else{
            double energy = 0.0;
            return energy;
        }
    }
    
};

#endif