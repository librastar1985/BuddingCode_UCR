#ifndef EDGESWAP_TEST_H_
#define EDGESWAP_TEST_H_

#include "SystemStructures.h"

class Edgeswap {
    std::vector<bool> boundary_node;
    std::vector<int> nndata;

    public:
    Edgeswap(CoordInfoVecs& coordInfoVecs, GeneralParams& generalParams);

	std::vector<bool> DomainBd (CoordInfoVecs& coordInfoVecs);
	std::vector<int> Number_of_Neighbor(CoordInfoVecs& coordInfoVecs);
    int edge_swap_device_vecs (int iedge, 
        GeneralParams& generalParams,
        CoordInfoVecs& coordInfoVecs,
        LinearSpringInfoVecs& linearSpringInfoVecs,
        BendingTriangleInfoVecs& bendingTriangleInfoVecs,
        AreaTriangleInfoVecs& areaTriangleInfoVecs);

    int edge_swap_host_vecs (int iedge, 
        GeneralParams& generalParams,
        HostSetInfoVecs& hostSetInfoVecs,
        LinearSpringInfoVecs& linearSpringInfoVecs,
        BendingTriangleInfoVecs& bendingTriangleInfoVecs,
        AreaTriangleInfoVecs& areaTriangleInfoVecs);
	
    void transferHtoD(CoordInfoVecs& coordInfoVecs,
        HostSetInfoVecs& hostSetInfoVecs);
    
    void transferDtoH(CoordInfoVecs& coordInfoVecs,
        HostSetInfoVecs& hostSetInfoVecs);

};

#endif