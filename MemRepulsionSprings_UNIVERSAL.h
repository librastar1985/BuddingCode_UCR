
#ifndef MEMREPULSIONSPRINGS_H_
#define MEMREPULSIONSPRINGS_H_ 

#include "SystemStructures.h"

void ComputeMemRepulsionSprings(
    CoordInfoVecs& coordInfoVecs,
	LinearSpringInfoVecs& linearSpringInfoVecs, 
	CapsidInfoVecs& capsidInfoVecs,
    GeneralParams& generalParams,
    AuxVecs& auxVecs);

struct MemRepulsionSpringFunctor : public thrust::unary_function<U2CVec3,CVec3> {
    double Rmin;
    double abs_Rmin;
    double epsilon_rep1;
    double epsilon_rep2;
    unsigned membraneMaxNode;

    double* membraneNodeXAddr;
    double* membraneNodeYAddr;
    double* membraneNodeZAddr;

    unsigned* id_value_expanded;
    unsigned* keyBegin;
    unsigned* keyEnd;
    
	__host__ __device__ 
    MemRepulsionSpringFunctor(
        double& _Rmin,
        double& _abs_Rmin,
        double& _epsilon_rep1,
        double& _epsilon_rep2,
        unsigned& _membraneMaxNode,
     
        double* _membraneNodeXAddr,
        double* _membraneNodeYAddr,
        double* _membraneNodeZAddr,

        unsigned* _id_value_expanded,
        unsigned* _keyBegin,
        unsigned* _keyEnd):

        Rmin(_Rmin),
        abs_Rmin(_abs_Rmin),
        epsilon_rep1(_epsilon_rep1),
        epsilon_rep2(_epsilon_rep2),
        membraneMaxNode(_membraneMaxNode),
    
        membraneNodeXAddr(_membraneNodeXAddr),
        membraneNodeYAddr(_membraneNodeYAddr),
        membraneNodeZAddr(_membraneNodeZAddr),

        id_value_expanded(_id_value_expanded),
        keyBegin(_keyBegin),
        keyEnd(_keyEnd) {}

	//hand in counting iterator id for membrane nodes
	__device__ 
    CVec3 operator()(const U2CVec3& u2d3) {

        unsigned node_id = thrust::get<0>(u2d3);
        
        unsigned bucketId = thrust::get<1>(u2d3);//bucket containing nodeId
	
		//beginning and end of attempted attachment id's in id_value_expanded
        //these indices are membrane
	
        double xLoc_LR;
        double yLoc_LR;
        double zLoc_LR;

        //extract current force
        double current_forceX = thrust::get<2>(u2d3);
        double current_forceY = thrust::get<3>(u2d3);
        double current_forceZ = thrust::get<4>(u2d3);

		double posX = membraneNodeXAddr[node_id];
        double posY = membraneNodeYAddr[node_id];
        double posZ = membraneNodeZAddr[node_id];

        //for now iterate through all membrane id's
       // unsigned begin = keyBegin[bucketId];
       // unsigned end = keyEnd[bucketId];

        for (unsigned memId_count = 0; memId_count < membraneMaxNode; memId_count++ ){
            
            unsigned memId = memId_count;//id_value_expanded[ memId_count ];

            if ((memId < membraneMaxNode) && (memId != node_id)) {
                //calculate distance
                xLoc_LR = -( posX - membraneNodeXAddr[memId] ); //This is like Rk - Rj, where Rj is the point we are trying to update
                yLoc_LR = -( posY - membraneNodeYAddr[memId] );
                zLoc_LR = -( posZ - membraneNodeZAddr[memId] );

                double R = sqrt( 
                            (xLoc_LR) * (xLoc_LR) + 
                            (yLoc_LR) * (yLoc_LR) + 
                            (zLoc_LR) * (zLoc_LR) );

                /*if (R < Rmin) {
                    double magnitude = 2*epsilon_rep1*
                                        (1-exp(-epsilon_rep2*(R-Rmin)))*
                                        (-exp(-epsilon_rep2*(R-Rmin)))*
                                        (epsilon_rep2/R);

                    current_forceX += -magnitude*xLoc_LR;
                    current_forceY += -magnitude*yLoc_LR;
                    current_forceZ += -magnitude*zLoc_LR;
                }*/
                if( R < abs_Rmin){
                   current_forceX += (epsilon_rep1)*(R - abs_Rmin)*(xLoc_LR)/R;
					current_forceY += (epsilon_rep1)*(R - abs_Rmin)*(yLoc_LR)/R;
					current_forceZ += (epsilon_rep1)*(R - abs_Rmin)*(zLoc_LR)/R;
                }
            }
        }
        return thrust::make_tuple(current_forceX, current_forceY, current_forceZ);
    }
};

#endif