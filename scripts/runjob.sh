#!/bin/bash

# this script submits a job on Discover


# job naming convention: ATLAS_[crot]_[cmlon]_[current grid]_[total grids]

#SBATCH --job-name=ATLAS_2099_180_EX3
#SBATCH -o output.%j
#SBATCH -e error.%j
#SBATCH --account=s1430

# what type of nodes to use
#SBATCH -C west
# how many nodes
#SBATCH -N 64
# how many tasks per node, used in scheduling
#SBATCH --ntasks-per-node=2
# how long
#SBATCH --time=1:00:00

# so i can run any commands
. /usr/share/modules/init/bash

# load appropriate modules
module load comp/intel-14.0.0.080
module load other/mpi/openmpi/1.7.3-intel-14.0.0.080

# move to correct directory
cd /discover/nobackup/bjgreer/2099-180/

# set number of threads (cores per node / mpi tasks per node)
export OMP_NUM_THREADS=6

echo "Starting job.."
mpirun -np 128 -npernode 2 ./atlas -v -r 1 -ml doplist -ts16 -grid grid_file -maxtiles 2 -loaddops 32 -outdir OUT/ -bk hmi.V_avg120.2099.180.mean.fits
echo "Ending job."
exit 0


# timing notes:
# 32 nodes, 2 tasks per node, 6 threads per task, 2 tiles per task = ~3200 tiles in 6 hours (530/hr, 8.9/min, 0.15/sec)
# 32 nodes, 2 tasks per node, 6 threads per task, 3 tiles per task = ~2100 tiles in 6 hours (350/hr)
# 64 nodes, 2 tasks per node, 6 threads per task, 2 tiles per task = ~5400 tiles in 6 hours (900/hr) <= best so far
# 128nodes, 2 tasks per node, 6 threads per task, 2 tiles per task = ~2600 tiles in 3 hours (887/hr)
