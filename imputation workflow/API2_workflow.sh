#!/bin/bash
#SBATCH -p medium
#SBATCH -c 4
#SBATCH -C scratch2
#SBATCH --mem=8G
#SBATCH -t 0-02:00:00
#SBATCH -J API2_workflow
#SBATCH -o API2_workflow.out
#SBATCH -e API2_workflow.err
#SBATCH --array=1-7

### Initialization
# Get Array ID
i=${SLURM_ARRAY_TASK_ID}
NAME_PARAMETER_FILE="API2_example_file.xlsx"

# Load modules and environment
module load r/4.0.3 
module load python/3.8.6
source /scratch2/users/tniehof/imputation_env/bin/activate

# Run the script
Rscript /usr/users/tniehof/imputation/API2_workflow.R $i $NAME_PARAMETER_FILE
echo '******************** FINISHED ***********************'
echo
