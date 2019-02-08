#!/bin/sh
### Job name
#PBS -N Vortex_dxdt_cpp16
### Number of nodes and processors per node.
#PBS -l nodes=1:ppn=1,mem=2G,walltime=300:00:00
#PBS -j oe
#PBS -m n


cd $PBS_O_WORKDIR
echo Host: $HOSTNAME
####echo Date: $(date)
echo Dir: $PWD
echo This jobs runs on the following processors:
cat $PBS_O_WORKDIR
# openmpi test

#matlab -nojvm -nodisplay -r computeOscMass >& log
#matlab -nojvm -nodisplay -r reinitscript >& log
#matlab -nojvm -nodisplay -r PGALS2DTest >& log
g++ FinalImplementation_lakshmanboundarycondition.cpp
./a.out >& log
