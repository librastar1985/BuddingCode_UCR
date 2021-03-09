# BuddingCode_UCR
Budding model. The initial condition can be a closed system or an open system.

The general flow of simulation steps can be found at "void System::solveSystem()" in System.cu.
For particular functions such as linear spring, please refer to LinearSpring.cu and LinearSpring.h files. The same applies to bending spring and area springs.
For edge-swap algorithm and general data structure manipulation functions, please refer to Edgeswap_test.cpp and Edgeswap_test.h.
To change the name of saved animation and data output, please refer to Storage.cpp and Storage.h.
To change the simulation job title and to some extend simulation time step size, please refer to SBATCH.sh.
Initial data structure (built via MATLAB functions) is located in Data_Structure.xml.

Overall flow of the simulation steps:
1. Initialization of global parameters and data structures.
2. Run a predetermined number of relaxation steps of the model system to attain quasi-steady state.
3. Start the actual simulation:

 a. Update parameters (if necessary).
 b. Run a predetermined number of relaxation steps or dynamical number of relaxation steps depending on simulation types (molecular dynamics vs energy minimization).
 c. Run edge-swap algorithm.
 d. Repeat (a-c) for a number of times, and then test for growth (if applicable).
 e. Repeat (a-d) until simulation terminates.

To run the simulation on UCR HPCC:
1. Make sure you have every file & folder in this repository on your HPCC account.
2. type: cd [name of the folder holding the file SBATCH.sh] (System.cu, etc should also be there as well).
3. type: module load extra; module load GCC; module load cuda/9.1
4. type: make (Or make -j N, N can be 2,3,4,....,or 12. But this should only be done if you are using an interactive gpu session. See UCR HPCC website for detail)
5. type: sbatch -p gpu --gres=gpu:1 --time=144:00:00 SBATCH.sh 
