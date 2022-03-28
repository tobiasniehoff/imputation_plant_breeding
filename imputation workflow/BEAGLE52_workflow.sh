#!/bin/bash
#SBATCH -p medium
#SBATCH -c 4
#SBATCH -C scratch2
#SBATCH --mem=8G
#SBATCH -t 0-02:00:00
#SBATCH -J BEAGLE52_workflow
#SBATCH -o BEAGLE52_workflow.out
#SBATCH -e BEAGLE52_workflow.err
#SBATCH --array=1-7

### Initialization
# Get Array ID
i=${SLURM_ARRAY_TASK_ID}
NAME_PARAMETER_FILE="BEAGLE52_example_file.xlsx"

# Run the script
module load r/4.0.3
module load openjdk/1.8.0_222-b10
source /scratch2/users/tniehof/imputation_env/bin/activate

Rscript /usr/users/tniehof/imputation/BEAGLE52_workflow.R $i $NAME_PARAMETER_FILE
echo $SLURM_ARRAY_TASK_ID
echo '******************** FINISHED ***********************'
echo
