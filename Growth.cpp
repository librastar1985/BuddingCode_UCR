#include "System.h"
#include <random>
#include "Growth.h"
#include <math.h>
#include <vector>

void Growth::growth(
    int ielem, 
    GeneralParams& generalParams,
    CoordInfoVecs& coordInfoVecs,
    AreaTriangleInfoVecs& areaTriangleInfoVecs) {
        int current_num_of_triangles = coordInfoVecs.num_triangles;
        bool triggered = false;
        int triggered_counter = 0;
        for (int p = 0; p < current_num_of_triangles; p++){
        int k = p;
        k -= triggered_counter;
        int first = coordInfoVecs.triangles2Nodes_1[k];
        int second = coordInfoVecs.triangles2Nodes_2[k];
        int third = coordInfoVecs.triangles2Nodes_3[k];
        double v1x = coordInfoVecs.nodeLocX[second] - coordInfoVecs.nodeLocX[first];
        double v1y = coordInfoVecs.nodeLocY[second] - coordInfoVecs.nodeLocY[first];
        double v1z = coordInfoVecs.nodeLocZ[second] - coordInfoVecs.nodeLocZ[first];
        double v2x = coordInfoVecs.nodeLocX[third] - coordInfoVecs.nodeLocX[first];
        double v2y = coordInfoVecs.nodeLocY[third] - coordInfoVecs.nodeLocY[first];
        double v2z = coordInfoVecs.nodeLocZ[third] - coordInfoVecs.nodeLocZ[first];
        double This_area = sqrt((v1y*v2z - v2y*v1z)*(v1y*v2z - v2y*v1z) + 
                                ((-v1x*v2z) + (v2x*v1z))*((-v1x*v2z) + (v2x*v1z)) +
                                (v1x*v2y - v2x*v1y)*(v1x*v2y - v2x*v1y))/2.0;
        if (This_area/ areaTriangleInfoVecs.initial_area > 1.5){
            //triggered = true;
            triggered_counter += 1;
        }
        else{continue;}
        int ielem = k;
        int te1 = coordInfoVecs.triangles2Edges_1[ielem];
        int te2 = coordInfoVecs.triangles2Edges_2[ielem];
        int te3 = coordInfoVecs.triangles2Edges_3[ielem];

        int tn1 = coordInfoVecs.edges2Nodes_1[te1];
        int tn2 = coordInfoVecs.edges2Nodes_2[te1];
        int tn3;
        if (coordInfoVecs.edges2Nodes_1[te2] == tn2){
            tn3 = coordInfoVecs.edges2Nodes_2[te2];
        }
        else if (coordInfoVecs.edges2Nodes_1[te2] == tn1){
            tn3 = coordInfoVecs.edges2Nodes_2[te2];
        }
        else if (coordInfoVecs.edges2Nodes_2[te2] == tn2){
            tn3 = coordInfoVecs.edges2Nodes_1[te2];
        }
        else if (coordInfoVecs.edges2Nodes_2[te2] == tn1){
            tn3 = coordInfoVecs.edges2Nodes_1[te2];
        }
        //These extract the indices of vertices of the selected triangle "ielem"



        double newx = (coordInfoVecs.nodeLocX[tn1] + coordInfoVecs.nodeLocX[tn2] + coordInfoVecs.nodeLocX[tn3])/3.0;
        coordInfoVecs.nodeLocX.push_back(newx);
        double newy = (coordInfoVecs.nodeLocY[tn1] + coordInfoVecs.nodeLocY[tn2] + coordInfoVecs.nodeLocY[tn3])/3.0;
        coordInfoVecs.nodeLocY.push_back(newy);
        double newz = (coordInfoVecs.nodeLocZ[tn1] + coordInfoVecs.nodeLocZ[tn2] + coordInfoVecs.nodeLocZ[tn3])/3.0;
        coordInfoVecs.nodeLocZ.push_back(newz);
        //These are the coordinate of the new vertex. Its index is "coordInfoVecs.nodeLocX.size()+1"
        generalParams.maxNodeCount += 1;
        if ((generalParams.nodes_in_upperhem[coordInfoVecs.nodeLocX[tn1]] == 1) && (generalParams.nodes_in_upperhem[coordInfoVecs.nodeLocX[tn2]] == 1) &&
        (generalParams.nodes_in_upperhem[coordInfoVecs.nodeLocX[tn3]] == 1)){
            generalParams.nodes_in_upperhem.push_back(1);
        }
        else{
            generalParams.nodes_in_upperhem.push_back(-1);
        }



        coordInfoVecs.triangles2Nodes_1.erase(coordInfoVecs.triangles2Nodes_1.begin() + ielem);
        coordInfoVecs.triangles2Nodes_2.erase(coordInfoVecs.triangles2Nodes_2.begin() + ielem);
        coordInfoVecs.triangles2Nodes_3.erase(coordInfoVecs.triangles2Nodes_3.begin() + ielem);
        //Delete the associated vertices information of the selected triangle.
        //Since we delete the chosen triangles, any triangle indexed lower than the deleted one will have its index reduced (or moved up) by 1.
        //Hence, we need to sweep through all data structures using the triangle index to change the index accordingly.
        for (int i = 0; i < coordInfoVecs.num_edges; i++){
            if (coordInfoVecs.edges2Triangles_1[i] > ielem){
                coordInfoVecs.edges2Triangles_1[i] = coordInfoVecs.edges2Triangles_1[i] - 1;
            }
            if (coordInfoVecs.edges2Triangles_2[i] > ielem){
                coordInfoVecs.edges2Triangles_2[i] = coordInfoVecs.edges2Triangles_2[i] - 1;
            }
        }
        //This completes the sweep. After this, the indices of triangle used in edges2Triangles data structure should be the correct one.
        coordInfoVecs.num_triangles += -1;


        int NODESIZE= generalParams.maxNodeCount;//coordInfoVecs.nodeLocX.size();
        coordInfoVecs.triangles2Nodes_1.push_back(tn1);
        coordInfoVecs.triangles2Nodes_2.push_back(tn2);
        coordInfoVecs.triangles2Nodes_3.push_back(NODESIZE-1);
        //This is a new triangle associated with (tn1, tn2, newnode). Its index is "coordinfoVecs.triangles2Nodes_1.size()-2".
        coordInfoVecs.triangles2Nodes_1.push_back(tn2);
        coordInfoVecs.triangles2Nodes_2.push_back(tn3);
        coordInfoVecs.triangles2Nodes_3.push_back(NODESIZE-1);
        //This is a new triangle associated with (tn2, tn3, newnode). Its index is "coordinfoVecs.triangles2Nodes_1.size()-1".
        coordInfoVecs.triangles2Nodes_1.push_back(tn3);
        coordInfoVecs.triangles2Nodes_2.push_back(tn1);
        coordInfoVecs.triangles2Nodes_3.push_back(NODESIZE-1);
        //This is a new triangle associated with (tn3, tn1, newnode). Its index is "coordinfoVecs.triangles2Nodes_1.size()".
        coordInfoVecs.num_triangles += 3;



        //Now we add new edges formed by the addition of the new node.
        coordInfoVecs.edges2Nodes_1.push_back(NODESIZE-1);
        coordInfoVecs.edges2Nodes_2.push_back(tn1);
        //This is a new edge associated with (newnode, tn1). Its index is "edges2Nodes_1.size()-2".
        coordInfoVecs.edges2Nodes_1.push_back(NODESIZE-1);
        coordInfoVecs.edges2Nodes_2.push_back(tn2);
        //This is a new edge associated with (newnode, tn2). Its index is "edges2Nodes_1.size()-1".
        coordInfoVecs.edges2Nodes_1.push_back(NODESIZE-1);
        coordInfoVecs.edges2Nodes_2.push_back(tn3);
        //This is a new edge associated with (newnode, tn3). Its index is "edges2Nodes_1.size()".
        coordInfoVecs.num_edges += 3;



        //Now we update the edges2Triangles data structure with new edges.
        int edgeindex, a, a1, a2, a3, temp1, temp2;
        int TRIANGLESIZE = coordInfoVecs.num_triangles;//coordInfoVecs.triangles2Nodes_1.size();
        if (coordInfoVecs.edges2Triangles_1[te1] == ielem){
            coordInfoVecs.edges2Triangles_1[te1] = TRIANGLESIZE-3;
        }
        else if (coordInfoVecs.edges2Triangles_2[te1] == ielem){
            coordInfoVecs.edges2Triangles_2[te1] = TRIANGLESIZE-3;
        }
        else{}
        if (coordInfoVecs.edges2Triangles_1[te2] == ielem){
            coordInfoVecs.edges2Triangles_1[te2] = TRIANGLESIZE-2;
        }
        else if (coordInfoVecs.edges2Triangles_2[te2] == ielem){
            coordInfoVecs.edges2Triangles_2[te2] = TRIANGLESIZE-2;
        }
        else{}
        if (coordInfoVecs.edges2Triangles_1[te3] == ielem){
            coordInfoVecs.edges2Triangles_1[te3] = TRIANGLESIZE-1;
        }
        else if (coordInfoVecs.edges2Triangles_2[te3] == ielem){
            coordInfoVecs.edges2Triangles_2[te3] = TRIANGLESIZE-1;
        }
        else{}
        //The above change the existing edges2Triangles data structure to accomodate new triangles added.

        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        for (int j = 0; j < 3; j++){
            //Now we will add the new data associated with both new edges and triangles.
            coordInfoVecs.edges2Triangles_1.push_back(TRIANGLESIZE-(3-j));
            if (j == 0){
                coordInfoVecs.edges2Triangles_2.push_back(TRIANGLESIZE-1);
            }
            else if (j == 1){
                coordInfoVecs.edges2Triangles_2.push_back(TRIANGLESIZE-3);
            }
            else if (j == 2){
                coordInfoVecs.edges2Triangles_2.push_back(TRIANGLESIZE-2);
            }
            //Now we check to see if the order of "push_back" done is correct, i.e. are edges2Triangles data in correct orientation.
            //This is crucial in the bendingspring computation.
            edgeindex = (coordInfoVecs.num_edges - (3-j));
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
                //std::cout<<"SWITCHED?"<<std::endl;
            }
            else{}
            //This checks if the orientation is correct or not, if not, flip the ordering.
        }

        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        //Now we will take care of the last unedited data structure "triangles2Edges".
        //int aa, bb;
        int EDGESIZE = coordInfoVecs.num_edges;//coordInfoVecs.edges2Nodes_1.size();
        for (int j = 0; j < 3; j++){
            if (j == 0){
                coordInfoVecs.triangles2Edges_1.push_back(EDGESIZE-3);
                coordInfoVecs.triangles2Edges_2.push_back(te1);
                coordInfoVecs.triangles2Edges_3.push_back(EDGESIZE-2);   
            }
            else if (j == 1){
                coordInfoVecs.triangles2Edges_1.push_back(EDGESIZE-2);
                coordInfoVecs.triangles2Edges_2.push_back(te2);
                coordInfoVecs.triangles2Edges_3.push_back(EDGESIZE-1);   
            }
            else if (j ==2){
                coordInfoVecs.triangles2Edges_1.push_back(EDGESIZE-1);
                coordInfoVecs.triangles2Edges_2.push_back(te3);
                coordInfoVecs.triangles2Edges_3.push_back(EDGESIZE-3);   
            }
            
        }
        coordInfoVecs.triangles2Edges_1.erase(coordInfoVecs.triangles2Edges_1.begin() + ielem);
        coordInfoVecs.triangles2Edges_2.erase(coordInfoVecs.triangles2Edges_2.begin() + ielem);
        coordInfoVecs.triangles2Edges_3.erase(coordInfoVecs.triangles2Edges_3.begin() + ielem);
//Erase the edge infomation related to the deleted triangle.

        }
//This should completes the dreadful data structure update associated with cell (membrane) growth. Have fun modifying it if you need more function!
};





