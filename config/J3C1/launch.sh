#!/bin/sh

TEMPLATE=run
xCPU=$1
yCPU=$2
zCPU=$3
OMP_NTHREADS=$4
TASKxNODE=$5
MPI=`expr $xCPU \* $yCPU \* $zCPU`

echo "xCPUs: $xCPU"
echo "yCPUs: $yCPU"
echo "zCPUs: $zCPU"
echo "OMP NUM THREADS: $OMP_NTHREADS"
echo "MPI: $MPI"

TASKxNODE=1 

DIMs=3

if test $OMP_NTHREADS -eq 0
then
TEMPLATE_CMD=$TEMPLATE"""_mpi.cmd"
else
TEMPLATE_CMD=$TEMPLATE"""_omp.cmd"
fi
echo "Template used: $TEMPLATE_CMD" 

run_tmp="run_tmp_"$MPI"_"$OMP_NTHREADS

sed 's/$1/'"$xCPU"'/
s/$2/'"$yCPU"'/
s/$3/'"$zCPU"'/
s/$4/'"$MPI"'/ 
s/$5/'"$OMP_NTHREADS"'/
s/$6/'"$TASKxNODE"'/ 
s/$7/'"$DIMs"'/' <$TEMPLATE_CMD> $run_tmp 


sbatch $run_tmp
