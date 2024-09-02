#! /bin/bash
#SBATCH --job-name="J3C1"                               
#SBATCH --workdir=.                                                                
#SBATCH --output=mpi_%j.out                                                                          
#SBATCH --error=mpi_%j.err                                                                               
#SBATCH --ntasks=$4                                                                                              
#SBATCH --cpus-per-task=$5                                                                                  
#SBATCH --tasks-per-node=$6                                                                                       
#SBATCH --time=01:00:00  
#SBATCH --qos=thin_astro

outRoot="/storage/scratch/lv72/lv72805/D10/J3C1/"
path_wkd="/home/lv72/lv72805/RHD_GC/J3C1/DATA"

export MP_IMPL=anl2
export LD_LIBRARY_PATH=/storage/apps/local/lib/:/storage/apps/SZIP/2.1.1/lib/:$LD_LIBRARY_PATH
export PATH=/storage/apps/HDF5/gcc/1.8.20/bin:$PATH

#module load gcc/4.6.1
#module load hdf5/1.8.20 intel/2018.3.222  impi/2018.3.222 mkl/2018.3.222
#export PATH=/storage/apps/HDF5/gcc/1.8.20/bin
module load hdf5/1.8.22_intel intel/2018.3.222  impi/2018.3.222 mkl/2018.3.222
export OMP_NUM_THREADS=$5

#date 
time /usr/bin/srun ./RATPENAT $outRoot $path_wkd $1 $2 $3 >out_$4_$5.dat
#date



