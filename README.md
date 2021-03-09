# BuddingCode_UCR
Budding model. The initial condition can be a closed system or an open system.

Case 1. yeast budding

      Starting from an initial spherical (or non-spherical) triangulated model cell, a selected region can undergo changes in mechanical properties leading to deformation. 
      Additional triangles can be introduced into the selected region to simulate local growth.

Case 2. Viral budding

      Starting from an initial flat triangulated sheet representing a piece of membrane, a selected region can undergo changes in mechanical properties leading to deformation.
      Interaction with a single viral particle or a cluster of viral particles can also be accommodated.


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

The general flow of simulation steps can be found at "void System::solveSystem()" in System.cu.

For particular functions such as linear spring, please refer to LinearSpring.cu and LinearSpring.h files. The same applies to bending spring and area springs.

For edge-swap algorithm and general data structure manipulation functions, please refer to Edgeswap_test.cpp and Edgeswap_test.h.

To change the name of saved animation and data output, please refer to Storage.cpp and Storage.h.

To change the simulation job title and to some extend simulation time step size, please refer to SBATCH.sh.

Initial data structure (built via MATLAB functions) is located in Data_Structure.xml.

Overall flow of the simulation steps:

I. Initialization of global parameters and data structures.

II. Run a predetermined number of relaxation steps of the model system to attain quasi-steady state.

III. Start the actual simulation:

      1. Update parameters (if necessary).
      2. Run a predetermined number of relaxation steps or dynamical number of relaxation steps depending on simulation types (molecular dynamics vs energy minimization).
      3. Run edge-swap algorithm.
      4. Repeat (a-c) for a number of times, and then test for growth (if applicable).
      5. Repeat (a-d) until simulation terminates.

To run the simulation on UCR HPCC:
1. Make sure you have every file & folder in this repository on your HPCC account.
2. type: cd [name of the folder holding the file SBATCH.sh] (System.cu, etc should also be there as well).
3. type: module load extra; module load GCC; module load cuda/9.1
4. type: make (Or make -j N, N can be 2,3,4,....,or 12. But this should only be done if you are using an interactive gpu session. See UCR HPCC website for detail)
5. type: sbatch -p gpu --gres=gpu:1 --time=144:00:00 SBATCH.sh 
