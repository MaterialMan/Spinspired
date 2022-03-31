#!/bin/bash
#SBATCH --job-name=STOhousing       # Job name
#SBATCH --mail-type=END,FAIL             # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=matt.dale@york.ac.uk     # Where to send mail  
#SBATCH --ntasks=1                          # Run a single task  
#SBATCH --cpus-per-task=1                  # Number of CPU cores per task 
#SBATCH --mem=20gb                        # Job memory request
#SBATCH --time=48:00:00                  # Time limit hrs:min:sec
#SBATCH --output=ch_%A_%a.log        # Standard output and error log
#SBATCH --account=cs-charcsi-2019        # Project account
#SBATCH --array=1-10

echo My working directory is `pwd`
echo Running job on host:
echo -e '\t'`hostname` at `date`
echo $SLURM_CPUS_ON_NODE CPU cores available

echo Running array job index $SLURM_ARRAY_TASK_ID, on host:
echo

cd '/mnt/lustre/users/md596/Magnets/multilayer/Support files/Simulators/vampire/Batches/batch0/Cores/core0'
chmod u+x vampire-serial 

cd '/mnt/lustre/users/md596/Magnets/multilayer/Evolve to Task/Single objective'

module load toolchain/foss/2020b
module load math/MATLAB/2020b

matlab -r 'MicrobialGAFcn($SLURM_ARRAY_TASK_ID,"STO",{36},"california_housing",10,1000,0)'
matlab -r 'MicrobialGAFcn($SLURM_ARRAY_TASK_ID,"STO",{64},"california_housing",10,1000,0)'

#matlab -r 'MicrobialGAFcn($SLURM_ARRAY_TASK_ID,"STO",{36},"california_housing",10,1000,1)'
#matlab -r 'MicrobialGAFcn($SLURM_ARRAY_TASK_ID,"STO",{64},"california_housing",10,1000,1)'

echo
echo Job completed at `date`
