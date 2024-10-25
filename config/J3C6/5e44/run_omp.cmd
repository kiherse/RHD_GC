#! /bin/bash
#SBATCH --job-name="J3C6_544"
#SBATCH --workdir=.
#SBATCH --output=/home/lv72/lv72805/RHD_GC/output/J3C6/5e44/mpi_%j.out
#SBATCH --error=/home/lv72/lv72805/RHD_GC/output/J3C6/5e44/mpi_%j.err
#SBATCH --ntasks=$4
#SBATCH --cpus-per-task=$5
#SBATCH --tasks-per-node=$6
#SBATCH --time=96:00:00
#SBATCH --mem=100G
#SBATCH --partition=genoa_s 
#SBATCH --qos=hera
# #SBATCH --qos=thin_astro
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kiara.hervella@uv.es

outRoot="/storage/scratch/lv72/lv72805/D06/J3C6/5e44/"
path_output="/home/lv72/lv72805/RHD_GC/output/J3C6/5e44/"
path_wkd="/home/lv72/lv72805/RHD_GC/config/J3C6/5e44"

# THIN QUEUE
# module load hdf5/1.14.1-2_ompi_gcc13.2 hwloc/2.7.1
# HERA QUEUE
module load hdf5/1.14.1-2_gcc13.2_ompi_rhel8 

make clean
make

export OMP_NUM_THREADS=$5

#date 
time /usr/bin/srun ./RATPENAT $outRoot $path_wkd $1 $2 $3 > ${path_output}/out_${SLURM_JOB_ID}_$4_$5.dat
#date
