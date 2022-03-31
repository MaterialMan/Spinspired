#!/bin/bash
#SBATCH --job-name=qgr250        # Job name
#SBATCH --mail-type=END,FAIL             # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=matt.dale@york.ac.uk     # Where to send mail  
#SBATCH --ntasks=1                          # Run a single task  
#SBATCH --cpus-per-task=10                  # Number of CPU cores per task 
#SBATCH --mem=20gb                        # Job memory request
#SBATCH --time=48:00:00                  # Time limit hrs:min:sec
#SBATCH --output=ch_%A_%a.log        # Standard output and error log
#SBATCH --account=cs-charcsi-2019        # Project account
#SBATCH --array=1

echo My working directory is `pwd`
echo Running job on host:
echo -e '\t'`hostname` at `date`
echo $SLURM_CPUS_ON_NODE CPU cores available

echo Running array job index $SLURM_ARRAY_TASK_ID, on host:
echo

cd '/mnt/lustre/users/md596/Magnets/MassMM/Support files/Simulators/vampire/Batches/batch0/Cores/core0'
chmod u+x vampire-serial 
chmod -R ugo+rwx /mnt/lustre/users/md596/Magnets/MassMM/

cd '/mnt/lustre/users/md596/Magnets/MassMM/'

module load toolchain/foss/2019a
module load math/MATLAB/2018a

matlab -r 'testRectWithQuad($SLURM_ARRAY_TASK_ID, "non_chan_eq",250)'

echo
echo Job completed at `date`
