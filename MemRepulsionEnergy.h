#ifndef MEMREPULSIONENERGY_H_
#define MEMREPULSIONENERGY_H_ 

#include "SystemStructures.h"

double ComputeMemRepulsionEnergy(
    CoordInfoVecs& coordInfoVecs,
	LinearSpringInfoVecs& linearSpringInfoVecs, 
	CapsidInfoVecs& capsidInfoVecs,
    GeneralParams& generalParams,
    AuxVecs& auxVecs);

struct MemRepulsionEnergyFunctor : public thrust::unary_function< UCVec3, double> {
    double Rmin;
    double abs_Rmin;
    double epsilon_rep1;
    double epsilon_rep2;
    int membraneMaxNode;

    int* nndata1;
    int* nndata2;
    int* nndata3;
    int* nndata4;
    int* nndata5;
    int* nndata6;
    int* nndata7;
    int* nndata8;
    int* nndata9;
    int* nndata10;
    int* nndata11;
    int* nndata12;

    double* membraneNodeXAddr;
    double* membraneNodeYAddr;
    double* membraneNodeZAddr;

    int* id_value_expanded;
    int* keyBegin;
    int* keyEnd;
    
	__host__ __device__ 
    MemRepulsionEnergyFunctor(
        double& _Rmin,
        double& _abs_Rmin,
        double& _epsilon_rep1,
        double& _epsilon_rep2,
        int& _membraneMaxNode,

        int* _nndata1,
        int* _nndata2,
        int* _nndata3,
        int* _nndata4,
        int* _nndata5,
        int* _nndata6,
        int* _nndata7,
        int* _nndata8,
        int* _nndata9,
        int* _nndata10,
        int* _nndata11,
        int* _nndata12,
     
        double* _membraneNodeXAddr,
        double* _membraneNodeYAddr,
        double* _membraneNodeZAddr,

        int* _id_value_expanded,
        int* _keyBegin,
        int* _keyEnd):

        Rmin(_Rmin),
        abs_Rmin(_abs_Rmin),
        epsilon_rep1(_epsilon_rep1),
        epsilon_rep2(_epsilon_rep2),
        membraneMaxNode(_membraneMaxNode),

        nndata1(_nndata1),
        nndata2(_nndata2),
        nndata3(_nndata3),
        nndata4(_nndata4),
        nndata5(_nndata5),
        nndata6(_nndata6),
        nndata7(_nndata7),
        nndata8(_nndata8),
        nndata9(_nndata9),
        nndata10(_nndata10),
        nndata11(_nndata11),
        nndata12(_nndata12),
    
        membraneNodeXAddr(_membraneNodeXAddr),
        membraneNodeYAddr(_membraneNodeYAddr),
        membraneNodeZAddr(_membraneNodeZAddr),

        id_value_expanded(_id_value_expanded),
        keyBegin(_keyBegin),
        keyEnd(_keyEnd) {}

	//hand in counting iterator id for membrane nodes
	__device__ 
    double operator()(const UCVec3& u1d3) {

        int node_id = thrust::get<0>(u1d3);
 
		double posX = thrust::get<1>(u1d3);
        double posY = thrust::get<2>(u1d3);
        double posZ = thrust::get<3>(u1d3);

        double xLoc_LR;
        double yLoc_LR;
        double zLoc_LR;

        
        int neighbor1 = nndata1[node_id];
        int neighbor2 = nndata2[node_id];
        int neighbor3 = nndata3[node_id];
        int neighbor4 = nndata4[node_id];
        int neighbor5 = nndata5[node_id];
        int neighbor6 = nndata6[node_id];
        int neighbor7 = nndata7[node_id];
        int neighbor8 = nndata8[node_id];
        int neighbor9 = nndata9[node_id];
        int neighbor10 = nndata10[node_id];
        int neighbor11 = nndata11[node_id];
        int neighbor12 = nndata12[node_id];

        //for now iterate through all membrane id's
        // int begin = keyBegin[bucketId];
        // int end = keyEnd[bucketId];
        double energy1 = 0.0;
        double energy2 = 0.0;
        double energy3 = 0.0;
        double energy4 = 0.0;
        double energy5 = 0.0;
        double energy6 = 0.0; 
        double energy7 = 0.0;
        double energy8 = 0.0;
        double energy9 = 0.0;
        double energy10 = 0.0;
        double energy11 = 0.0;
        double energy12 = 0.0;
        double energy = 0.0;
        int neighbor;
    for (int j = 0; j < 12; j ++){
        if (j == 0 && (nndata1[node_id] >= 0)){neighbor = nndata1[node_id];}
        //else{continue;}
        else if (j == 1 && (nndata2[node_id] >= 0)){neighbor = nndata2[node_id];}
        //else{continue;}
        else if (j == 2 && (nndata3[node_id] >= 0)){neighbor = nndata3[node_id];}
        //else{continue;}
        else if (j == 3 && (nndata4[node_id] >= 0)){neighbor = nndata4[node_id];}
        //else{continue;}
        else if (j == 4 && (nndata5[node_id] >= 0)){neighbor = nndata5[node_id];}
        //else{continue;}
        else if (j == 5 && (nndata6[node_id] >= 0)){neighbor = nndata6[node_id];}
        //else{continue;}
        else if (j == 6 && (nndata7[node_id] >= 0)){neighbor = nndata7[node_id];}
        //else{continue;}
        else if (j == 7 && (nndata8[node_id] >= 0)){neighbor = nndata8[node_id];}
        //else{continue;}
        else if (j == 8 && (nndata9[node_id] >= 0)){neighbor = nndata9[node_id];}
        //else{continue;}
        else if (j == 9 && (nndata10[node_id] >= 0)){neighbor = nndata10[node_id];}
        //else{continue;}
        else if (j == 10 && (nndata11[node_id] >= 0)){neighbor = nndata11[node_id];}
        //else{continue;}
        else if (j == 11 && (nndata12[node_id] >= 0)){neighbor = nndata12[node_id];}
        else{continue;}

        //for (int i = 0; i < 12; i++){
           if (nndata1[neighbor] >= 0 &&
                 nndata1[neighbor] != node_id &&
                 nndata1[neighbor] != neighbor1 &&
                 nndata1[neighbor] != neighbor2 &&
                 nndata1[neighbor] != neighbor3 &&
                 nndata1[neighbor] != neighbor4 &&
                 nndata1[neighbor] != neighbor5 &&
                 nndata1[neighbor] != neighbor6 &&
                 nndata1[neighbor] != neighbor7 &&
                 nndata1[neighbor] != neighbor8 &&
                 nndata1[neighbor] != neighbor9 &&
                 nndata1[neighbor] != neighbor10 &&
                 nndata1[neighbor] != neighbor11 &&
                 nndata1[neighbor] != neighbor12 
                 ){
                xLoc_LR = -( posX - membraneNodeXAddr[nndata1[neighbor]] ); //This is like Rk - Rj, where Rj is the point we are trying to update
                yLoc_LR = -( posY - membraneNodeYAddr[nndata1[neighbor]] );
                zLoc_LR = -( posZ - membraneNodeZAddr[nndata1[neighbor]] );
                double R = sqrt( 
                            (xLoc_LR) * (xLoc_LR) + 
                            (yLoc_LR) * (yLoc_LR) + 
                            (zLoc_LR) * (zLoc_LR) );
                if (R < abs_Rmin){
                    energy1 = epsilon_rep1*(1-exp(-epsilon_rep2*(R - abs_Rmin)))*(1-exp(-epsilon_rep2*(R - abs_Rmin)));
                }
                else{
                    energy1 = 0.0;
                }
           }
           else{
               energy1 = 0.0;
           }
           
           if (nndata2[neighbor] >= 0 &&
                 nndata2[neighbor] != node_id &&
                 nndata2[neighbor] != neighbor1 &&
                 nndata2[neighbor] != neighbor2 &&
                 nndata2[neighbor] != neighbor3 &&
                 nndata2[neighbor] != neighbor4 &&
                 nndata2[neighbor] != neighbor5 &&
                 nndata2[neighbor] != neighbor6 &&
                 nndata2[neighbor] != neighbor7 &&
                 nndata2[neighbor] != neighbor8 &&
                 nndata2[neighbor] != neighbor9 &&
                 nndata2[neighbor] != neighbor10 &&
                 nndata2[neighbor] != neighbor11 &&
                 nndata2[neighbor] != neighbor12 
                 ){
                xLoc_LR = -( posX - membraneNodeXAddr[nndata2[neighbor]] ); //This is like Rk - Rj, where Rj is the point we are trying to update
                yLoc_LR = -( posY - membraneNodeYAddr[nndata2[neighbor]] );
                zLoc_LR = -( posZ - membraneNodeZAddr[nndata2[neighbor]] );
                double R = sqrt( 
                            (xLoc_LR) * (xLoc_LR) + 
                            (yLoc_LR) * (yLoc_LR) + 
                            (zLoc_LR) * (zLoc_LR) );
                if (R < abs_Rmin){
                    energy2 = epsilon_rep1*(1-exp(-epsilon_rep2*(R - abs_Rmin)))*(1-exp(-epsilon_rep2*(R - abs_Rmin)));
                }
                else{
                    energy2 = 0.0;
                }
           }
           else{
               energy2 = 0.0;
           }
           
           if (nndata3[neighbor] >= 0 &&
                nndata3[neighbor] != node_id &&
                nndata3[neighbor] != neighbor1 &&
                nndata3[neighbor] != neighbor2 &&
                nndata3[neighbor] != neighbor3 &&
                nndata3[neighbor] != neighbor4 &&
                nndata3[neighbor] != neighbor5 &&
                nndata3[neighbor] != neighbor6 &&
                nndata3[neighbor] != neighbor7 &&
                nndata3[neighbor] != neighbor8 &&
                nndata3[neighbor] != neighbor9 &&
                nndata3[neighbor] != neighbor10 &&
                nndata3[neighbor] != neighbor11 &&
                nndata3[neighbor] != neighbor12 ){
                xLoc_LR = -( posX - membraneNodeXAddr[nndata3[neighbor]] ); //This is like Rk - Rj, where Rj is the point we are trying to update
                yLoc_LR = -( posY - membraneNodeYAddr[nndata3[neighbor]] );
                zLoc_LR = -( posZ - membraneNodeZAddr[nndata3[neighbor]] );
                double R = sqrt( 
                            (xLoc_LR) * (xLoc_LR) + 
                            (yLoc_LR) * (yLoc_LR) + 
                            (zLoc_LR) * (zLoc_LR) );
                if (R < abs_Rmin){
                    energy3 += epsilon_rep1*(1-exp(-epsilon_rep2*(R - abs_Rmin)))*(1-exp(-epsilon_rep2*(R - abs_Rmin)));
                    }
                    else{
                        energy3 = 0.0;
                    }
            }
           else{
                energy3 = 0.0;
           }
           
           if (nndata4[neighbor] >= 0 &&
                 nndata4[neighbor] != node_id &&
                 nndata4[neighbor] != neighbor1 &&
                 nndata4[neighbor] != neighbor2 &&
                 nndata4[neighbor] != neighbor3 &&
                 nndata4[neighbor] != neighbor4 &&
                 nndata4[neighbor] != neighbor5 &&
                 nndata4[neighbor] != neighbor6 &&
                 nndata4[neighbor] != neighbor7 &&
                 nndata4[neighbor] != neighbor8 &&
                 nndata4[neighbor] != neighbor9 &&
                 nndata4[neighbor] != neighbor10 &&
                 nndata4[neighbor] != neighbor11 &&
                 nndata4[neighbor] != neighbor12 ){
                xLoc_LR = -( posX - membraneNodeXAddr[nndata4[neighbor]] ); //This is like Rk - Rj, where Rj is the point we are trying to update
                yLoc_LR = -( posY - membraneNodeYAddr[nndata4[neighbor]] );
                zLoc_LR = -( posZ - membraneNodeZAddr[nndata4[neighbor]] );
                double R = sqrt( 
                            (xLoc_LR) * (xLoc_LR) + 
                            (yLoc_LR) * (yLoc_LR) + 
                            (zLoc_LR) * (zLoc_LR) );
            if (R < abs_Rmin){
                energy4 += epsilon_rep1*(1-exp(-epsilon_rep2*(R - abs_Rmin)))*(1-exp(-epsilon_rep2*(R - abs_Rmin)));
                    }
                    else{
                        energy4 = 0.0;
                    }
            }
           else{
               energy4 = 0.0;
           }
           
           if (nndata5[neighbor] >= 0 &&
                 nndata5[neighbor] != node_id &&
                 nndata5[neighbor] != neighbor1 &&
                 nndata5[neighbor] != neighbor2 &&
                 nndata5[neighbor] != neighbor3 &&
                 nndata5[neighbor] != neighbor4 &&
                 nndata5[neighbor] != neighbor5 &&
                 nndata5[neighbor] != neighbor6 &&
                 nndata5[neighbor] != neighbor7 &&
                 nndata5[neighbor] != neighbor8 &&
                 nndata5[neighbor] != neighbor9 &&
                 nndata5[neighbor] != neighbor10 &&
                 nndata5[neighbor] != neighbor11 &&
                 nndata5[neighbor] != neighbor12 ){
                xLoc_LR = -( posX - membraneNodeXAddr[nndata5[neighbor]] ); //This is like Rk - Rj, where Rj is the point we are trying to update
                yLoc_LR = -( posY - membraneNodeYAddr[nndata5[neighbor]] );
                zLoc_LR = -( posZ - membraneNodeZAddr[nndata5[neighbor]] );
                double R = sqrt( 
                            (xLoc_LR) * (xLoc_LR) + 
                            (yLoc_LR) * (yLoc_LR) + 
                            (zLoc_LR) * (zLoc_LR) );
            if (R < abs_Rmin){
                energy5 = epsilon_rep1*(1-exp(-epsilon_rep2*(R - abs_Rmin)))*(1-exp(-epsilon_rep2*(R - abs_Rmin)));
                    }
                    else{
                        energy5 = 0.0;
                    }
            }
           else{
               energy5 = 0.0;
           }
          
           if (nndata6[neighbor] >= 0 && 
                nndata6[neighbor] != node_id &&
                nndata6[neighbor] != neighbor1 &&
                 nndata6[neighbor] != neighbor2 &&
                 nndata6[neighbor] != neighbor3 &&
                 nndata6[neighbor] != neighbor4 &&
                 nndata6[neighbor] != neighbor5 &&
                 nndata6[neighbor] != neighbor6 &&
                 nndata6[neighbor] != neighbor7 &&
                 nndata6[neighbor] != neighbor8 &&
                 nndata6[neighbor] != neighbor9 &&
                 nndata6[neighbor] != neighbor10 &&
                 nndata6[neighbor] != neighbor11 &&
                 nndata6[neighbor] != neighbor12 ){
                xLoc_LR = -( posX - membraneNodeXAddr[nndata6[neighbor]] ); //This is like Rk - Rj, where Rj is the point we are trying to update
                yLoc_LR = -( posY - membraneNodeYAddr[nndata6[neighbor]] );
                zLoc_LR = -( posZ - membraneNodeZAddr[nndata6[neighbor]] );
                double R = sqrt( 
                            (xLoc_LR) * (xLoc_LR) + 
                            (yLoc_LR) * (yLoc_LR) + 
                            (zLoc_LR) * (zLoc_LR) );
            if (R < abs_Rmin){
                energy6 = epsilon_rep1*(1-exp(-epsilon_rep2*(R - abs_Rmin)))*(1-exp(-epsilon_rep2*(R - abs_Rmin)));
                    }
                    else{
                        energy6 = 0.0;
                    }
            }
           else{
               energy6 = 0.0;
           }
           
           if (nndata7[neighbor] >= 0 &&
                nndata7[neighbor] != node_id &&
                nndata7[neighbor] != neighbor1 &&
                 nndata7[neighbor] != neighbor2 &&
                 nndata7[neighbor] != neighbor3 &&
                 nndata7[neighbor] != neighbor4 &&
                 nndata7[neighbor] != neighbor5 &&
                 nndata7[neighbor] != neighbor6 &&
                 nndata7[neighbor] != neighbor7 &&
                 nndata7[neighbor] != neighbor8 &&
                 nndata7[neighbor] != neighbor9 &&
                 nndata7[neighbor] != neighbor10 &&
                 nndata7[neighbor] != neighbor11 &&
                 nndata7[neighbor] != neighbor12 ){
                xLoc_LR = -( posX - membraneNodeXAddr[nndata7[neighbor]] ); //This is like Rk - Rj, where Rj is the point we are trying to update
                yLoc_LR = -( posY - membraneNodeYAddr[nndata7[neighbor]] );
                zLoc_LR = -( posZ - membraneNodeZAddr[nndata7[neighbor]] );
                double R = sqrt( 
                            (xLoc_LR) * (xLoc_LR) + 
                            (yLoc_LR) * (yLoc_LR) + 
                            (zLoc_LR) * (zLoc_LR) );
            if (R < abs_Rmin){
                energy7 = epsilon_rep1*(1-exp(-epsilon_rep2*(R - abs_Rmin)))*(1-exp(-epsilon_rep2*(R - abs_Rmin)));
                    }
                    else{
                        energy7 = 0.0;
                    }
            }
           else{
               energy7 = 0.0;
           }
           
           if (nndata8[neighbor] >= 0 && 
                nndata8[neighbor] != node_id &&
                nndata8[neighbor] != neighbor1 &&
                nndata8[neighbor] != neighbor2 &&
                nndata8[neighbor] != neighbor3 &&
                nndata8[neighbor] != neighbor4 &&
                nndata8[neighbor] != neighbor5 &&
                nndata8[neighbor] != neighbor6 &&
                nndata8[neighbor] != neighbor7 &&
                nndata8[neighbor] != neighbor8 &&
                nndata8[neighbor] != neighbor9 &&
                nndata8[neighbor] != neighbor10 &&
                nndata8[neighbor] != neighbor11 &&
                nndata8[neighbor] != neighbor12 ){
                xLoc_LR = -( posX - membraneNodeXAddr[nndata8[neighbor]] ); //This is like Rk - Rj, where Rj is the point we are trying to update
                yLoc_LR = -( posY - membraneNodeYAddr[nndata8[neighbor]] );
                zLoc_LR = -( posZ - membraneNodeZAddr[nndata8[neighbor]] );
                double R = sqrt( 
                            (xLoc_LR) * (xLoc_LR) + 
                            (yLoc_LR) * (yLoc_LR) + 
                            (zLoc_LR) * (zLoc_LR) );
            if (R < abs_Rmin){
                energy8 = epsilon_rep1*(1-exp(-epsilon_rep2*(R - abs_Rmin)))*(1-exp(-epsilon_rep2*(R - abs_Rmin)));
                    }
                    else{
                        energy8 = 0.0;
                    }
            }
           else{
               energy8 = 0.0;
           }
           
           if (nndata9[neighbor] >= 0 && 
                nndata9[neighbor] != node_id &&
                nndata9[neighbor] != neighbor1 &&
                 nndata9[neighbor] != neighbor2 &&
                 nndata9[neighbor] != neighbor3 &&
                 nndata9[neighbor] != neighbor4 &&
                 nndata9[neighbor] != neighbor5 &&
                 nndata9[neighbor] != neighbor6 &&
                 nndata9[neighbor] != neighbor7 &&
                 nndata9[neighbor] != neighbor8 &&
                 nndata9[neighbor] != neighbor9 &&
                 nndata9[neighbor] != neighbor10 &&
                 nndata9[neighbor] != neighbor11 &&
                 nndata9[neighbor] != neighbor12 ){
                xLoc_LR = -( posX - membraneNodeXAddr[nndata9[neighbor]] ); //This is like Rk - Rj, where Rj is the point we are trying to update
                yLoc_LR = -( posY - membraneNodeYAddr[nndata9[neighbor]] );
                zLoc_LR = -( posZ - membraneNodeZAddr[nndata9[neighbor]] );
                double R = sqrt( 
                            (xLoc_LR) * (xLoc_LR) + 
                            (yLoc_LR) * (yLoc_LR) + 
                            (zLoc_LR) * (zLoc_LR) );
            if (R < abs_Rmin){
               energy9 = epsilon_rep1*(1-exp(-epsilon_rep2*(R - abs_Rmin)))*(1-exp(-epsilon_rep2*(R - abs_Rmin)));
                    }
                    else{
                        energy9 = 0.0;
                    }
            }
           else{
               energy9 = 0.0;
           }
           
           if (nndata10[neighbor] >= 0 && 
                nndata10[neighbor] != node_id &&
                nndata10[neighbor] != neighbor1 &&
                nndata10[neighbor] != neighbor2 &&
                nndata10[neighbor] != neighbor3 &&
                nndata10[neighbor] != neighbor4 &&
                nndata10[neighbor] != neighbor5 &&
                nndata10[neighbor] != neighbor6 &&
                nndata10[neighbor] != neighbor7 &&
                nndata10[neighbor] != neighbor8 &&
                nndata10[neighbor] != neighbor9 &&
                nndata10[neighbor] != neighbor10 &&
                nndata10[neighbor] != neighbor11 &&
                nndata10[neighbor] != neighbor12 ){
                xLoc_LR = -( posX - membraneNodeXAddr[nndata10[neighbor]] ); //This is like Rk - Rj, where Rj is the point we are trying to update
                yLoc_LR = -( posY - membraneNodeYAddr[nndata10[neighbor]] );
                zLoc_LR = -( posZ - membraneNodeZAddr[nndata10[neighbor]] );
                double R = sqrt( 
                            (xLoc_LR) * (xLoc_LR) + 
                            (yLoc_LR) * (yLoc_LR) + 
                            (zLoc_LR) * (zLoc_LR) );
            if (R < abs_Rmin){
                energy10 = epsilon_rep1*(1-exp(-epsilon_rep2*(R - abs_Rmin)))*(1-exp(-epsilon_rep2*(R - abs_Rmin)));
                    }
                    else{
                        energy10 = 0.0;
                    }
            }
           else{
               energy10 = 0.0;
           }
          
           if (nndata11[neighbor] >= 0 && 
                nndata11[neighbor] != node_id &&
                nndata11[neighbor] != neighbor1 &&
                nndata11[neighbor] != neighbor2 &&
                nndata11[neighbor] != neighbor3 &&
                nndata11[neighbor] != neighbor4 &&
                nndata11[neighbor] != neighbor5 &&
                nndata11[neighbor] != neighbor6 &&
                nndata11[neighbor] != neighbor7 &&
                nndata11[neighbor] != neighbor8 &&
                nndata11[neighbor] != neighbor9 &&
                nndata11[neighbor] != neighbor10 &&
                nndata11[neighbor] != neighbor11 &&
                nndata11[neighbor] != neighbor12 ){
                xLoc_LR = -( posX - membraneNodeXAddr[nndata11[neighbor]] ); //This is like Rk - Rj, where Rj is the point we are trying to update
                yLoc_LR = -( posY - membraneNodeYAddr[nndata11[neighbor]] );
                zLoc_LR = -( posZ - membraneNodeZAddr[nndata11[neighbor]] );
                double R = sqrt( 
                            (xLoc_LR) * (xLoc_LR) + 
                            (yLoc_LR) * (yLoc_LR) + 
                            (zLoc_LR) * (zLoc_LR) );
            if (R < abs_Rmin){
                energy11 = epsilon_rep1*(1-exp(-epsilon_rep2*(R - abs_Rmin)))*(1-exp(-epsilon_rep2*(R - abs_Rmin)));
                    }
                    else{
                        energy11 = 0.0;
                    }
            }
           else{
               energy11 = 0.0;
           }
           
           if (nndata12[neighbor] >= 0 && 
                nndata12[neighbor] != node_id &&
                nndata12[neighbor] != neighbor1 &&
                nndata12[neighbor] != neighbor2 &&
                nndata12[neighbor] != neighbor3 &&
                nndata12[neighbor] != neighbor4 &&
                nndata12[neighbor] != neighbor5 &&
                nndata12[neighbor] != neighbor6 &&
                nndata12[neighbor] != neighbor7 &&
                nndata12[neighbor] != neighbor8 &&
                nndata12[neighbor] != neighbor9 &&
                nndata12[neighbor] != neighbor10 &&
                nndata12[neighbor] != neighbor11 &&
                nndata12[neighbor] != neighbor12 ){
                xLoc_LR = -( posX - membraneNodeXAddr[nndata12[neighbor]] ); //This is like Rk - Rj, where Rj is the point we are trying to update
                yLoc_LR = -( posY - membraneNodeYAddr[nndata12[neighbor]] );
                zLoc_LR = -( posZ - membraneNodeZAddr[nndata12[neighbor]] );
                double R = sqrt( 
                            (xLoc_LR) * (xLoc_LR) + 
                            (yLoc_LR) * (yLoc_LR) + 
                            (zLoc_LR) * (zLoc_LR) );
            if (R < abs_Rmin){
                energy12 = epsilon_rep1*(1-exp(-epsilon_rep2*(R - abs_Rmin)))*(1-exp(-epsilon_rep2*(R - abs_Rmin)));
                    }
                    else{
                        energy12 = 0.0;
                    }
            }
           else{
               energy12 = 0.0;
           }

        energy += energy1 + energy2 + energy3 + energy4 + energy5 + energy6 + energy7 + energy8 + energy9 + energy10 + energy11 + energy12; 
       //}
    }
        return energy;
    }
};

#endif

