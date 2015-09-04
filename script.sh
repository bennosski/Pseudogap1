#!/bin/bash

#tell grid engine to use current directory
#$ -cwd


## the "meat" of the script

#just print the name of this machine

mpirun -np 24 a.out>log
