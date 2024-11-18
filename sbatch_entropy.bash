#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=64

# Bash script for running MSE on the supercomputer for the Entropy project
#
# 20241113SM - created; run 6 and fix  (as entropy_onSubject_PSC.sh)
# 20241116WF - remove double run; add EPOCH and LENGTH (as sbatch_entropy.bash)

module load matlab
[ ! -r $INPUTFILE ] && echo "ERROR: no '$INPUTFILE'" && exit 1
[[ "${EPOCH}" =~ delay|fix ]] && echo "ERROR: epoch '$EPOCH' not 'delay' or 'fix'" && exit 1
[ -z "$LENGTH" ] && echo "ERROR: missing \$LENGHT not set" && exit 1

# matlab will test if file doesn't exit
matlab -nodesktop -r "run_Entropy_PSC('$INPUTFILE','MGS', '$EPOCH', $LENGTH);"

