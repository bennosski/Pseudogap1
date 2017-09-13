This code calculates pseudogap selfenergies.
Uses MPI
This version uses Norman Bandstructure


Compile with
mpic++ -fopenmp -I ../Eigen/eigen... main.cpp selfenergy.cpp globals.cpp

On Sherlock:
load module openmpi/1.8.4-nodlopen/gcc


main.cpp
calculate pseudogap selfenergy at one k point

cut.cpp 
calculate pseudogap selfenergy at cut of k points

FermiSurfaceMap.cpp
calculate pseudogap selfenergy at w=0 and all k points


####MPI problem
mpirun -np 24 a.out will run the code on only one node
which can only have up to 8 cores on corn and 24 cores on barley! 

Also remember to submit a job to barley specify 48 in script.sh and run:
qsub -pe orte 48 script.sh
