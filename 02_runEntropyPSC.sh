#!/usr/bin/env bash
set -euo pipefail

# 20241116WF - rewrite of allSubs_PSC.sh to loop over epohc and length

# to see what this would do but not actually do it, run like
#  DRYRUN=1 ./02_runEntropyPSC.sh
[ -n "${DRYRUN:-}" ] && DRYRUN="echo" || DRYRUN=

outdir=$PWD/Results/MGS_Entropy/individual_subject_files
logname=$PWD/log/%j-%x.log # job number and job name

echo "# [$(date)] collecting set file list"
ALLFILES=(data/ICAwholeClean_homogenize/*.set)
echo "# [$(date)] running for ${#ALLFILES[@]} set files"

for inputfile in "${ALLFILES[@]}"; do 
  ld8=$(grep -oP "\d{5}_\d{8}" <<< "$inputfile")

  for EPOCH in "delay" "fix"; do

     # epoch determins what lengths there are
     case $EPOCH in 
       fix)   lengths=(0);;
       delay) lengths=(6);; # TODO: maybe add 8 and 10
     esac

     for LENGTH in "${lengths[@]}"; do
         job_name=${ld8}_${EPOCH}_$LENGTH
         output=$outdir/${ld8}_MultiScaleEntropy_$EPOCH$LENGTH.csv

         #TODO: does fix0 need to be just 'fix'. if fix0 is okay, don't need below
         [ "$EPOCH$LENGTH" == fix0 ] && output=$outdir/${ld8}_MultiScaleEntropy_$EPOCH.csv

         test -r "$output" && echo "# skipping $ld8 $EPOCH $LENGTH: have '$output'" && continue
         squeue -A soc230004p -o %j | grep "$job_name" &&
             echo "# skipping $ld8 $EPOCH $LENGTH already in queue as '$job_name'" && continue

         echo "# [$(date)] submitting '$job_name' to make '$output'"
         $DRYRUN sbatch -J $job_name -e $logname o $logname \
           --export="INPUTFILE=$inputfile,EPOCH=$EPOCH,LENGTH=$LENGTH" \
           sbatch_entropy.bash

         #TODO: remove this break when ready to run for all
         break 3
     done
  done
done

