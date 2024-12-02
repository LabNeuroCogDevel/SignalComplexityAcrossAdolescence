#!/bin/bash
#SBATCH --time=03:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24

# Bash script for running MSE on the supercomputer for the Entropy project
#
# 20241113SM - created; run 6 and fix  (as entropy_onSubject_PSC.sh)
# 20241116WF - remove double run; add EPOCH and LENGTH (as sbatch_entropy.bash)

[ ! -r "$INPUTFILE" ] && echo "ERROR: $0: no '$INPUTFILE'" && exit 1
! [[ "${EPOCH}" =~ delay|fix ]] && echo "ERROR: $0: epoch '$EPOCH' not 'delay' or 'fix'" && exit 1
[ -z "$LENGTH" ] && echo "ERROR: $0: \$LENGHT not set" && exit 1

# matlab will test if file doesn't exit
cmd="try, run_Entropy_PSC('$INPUTFILE','MGS', '$EPOCH', $LENGTH); catch e, disp(e), end; quit"
echo "# [$(date)] $0 running $cmd in matlab"
module load matlab
matlab -nosplash -nodesktop -r "$cmd"
echo "# [$(date)] finished matlab"

