#include "System.h"
#include "VolumeSprings.h"



//WARNING: function must not reset coordInfoVecs.nodeForceX etc. 
void ComputeVolumeSprings(
    CoordInfoVecs& coordInfoVecs,
	LinearSpringInfoVecs& linearSpringInfoVecs, 
	CapsidInfoVecs& capsidInfoVecs,
    GeneralParams& generalParams,
    AuxVecs& auxVecs) {    
    
    /*thrust::fill(coordInfoVecs.tempNodeForceXReduced.begin(),coordInfoVecs.tempNodeForceXReduced.end(),0.0);
    thrust::fill(coordInfoVecs.tempNodeForceYReduced.begin(),coordInfoVecs.tempNodeForceYReduced.end(),0.0);
    thrust::fill(coordInfoVecs.tempNodeForceZReduced.begin(),coordInfoVecs.tempNodeForceZReduced.end(),0.0);
    thrust::fill(coordInfoVecs.tempNodeForceXUnreduced.begin(),coordInfoVecs.tempNodeForceXUnreduced.end(),0.0);
    thrust::fill(coordInfoVecs.tempNodeForceYUnreduced.begin(),coordInfoVecs.tempNodeForceYUnreduced.end(),0.0);
    thrust::fill(coordInfoVecs.tempNodeForceZUnreduced.begin(),coordInfoVecs.tempNodeForceZUnreduced.end(),0.0);*/

    //CVec4 init(0.0, 0.0, 0.0, 0.0); 
    thrust::counting_iterator<int> begin(0);

    thrust::transform(  
        thrust::make_zip_iterator(
            thrust::make_tuple(
                begin,
                auxVecs.id_bucket.begin(),
                coordInfoVecs.nodeForceX.begin(),
                coordInfoVecs.nodeForceY.begin(),
                coordInfoVecs.nodeForceZ.begin())),
        
        thrust::make_zip_iterator(
            thrust::make_tuple(
                begin,
                auxVecs.id_bucket.begin(),
                coordInfoVecs.nodeForceX.begin(),
                coordInfoVecs.nodeForceY.begin(),
                coordInfoVecs.nodeForceZ.begin())) + generalParams.maxNodeCount,

        thrust::make_zip_iterator(
            thrust::make_tuple(
                coordInfoVecs.nodeForceX.begin(),
                coordInfoVecs.nodeForceY.begin(),
                coordInfoVecs.nodeForceZ.begin())),

        VolumeSpringFunctor(
                generalParams.current_total_volume,
                generalParams.true_current_total_volume,
                generalParams.eq_total_volume,
                generalParams.volume_spring_constant,
                coordInfoVecs.num_triangles,
                generalParams.Rmin,

                thrust::raw_pointer_cast(coordInfoVecs.triangles2Nodes_1.data()),
                thrust::raw_pointer_cast(coordInfoVecs.triangles2Nodes_2.data()),
                thrust::raw_pointer_cast(coordInfoVecs.triangles2Nodes_3.data()),
                thrust::raw_pointer_cast(coordInfoVecs.nodeLocX.data()),
                thrust::raw_pointer_cast(coordInfoVecs.nodeLocY.data()),
                thrust::raw_pointer_cast(coordInfoVecs.nodeLocZ.data()), 
                
                          
                thrust::raw_pointer_cast(auxVecs.id_value_expanded.data()),
                thrust::raw_pointer_cast(auxVecs.keyBegin.data()),
                thrust::raw_pointer_cast(auxVecs.keyEnd.data()))
            );
            
            
    //for (int i = 0; i < generalParams.maxNodeCount; i++){
     //   std::cout<<"Force from volume on node 36 = "<<coordInfoVecs.nodeForceX[35]<<" "<<coordInfoVecs.nodeForceY[35]<<" "<<coordInfoVecs.nodeForceZ[35]<<std::endl;
    //}
    //std::cout<<"VOLUME"<<std::endl;
};

