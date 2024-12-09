#!/usr/bin/env bash
set -euo pipefail

# 20241116WF - rewrite of allSubs_PSC.sh to loop over epohc and length
# 20241206WF - add EPOCHS; add $sbatch_args but not used
EPOCHS=("fix")

# to see what this would do but not actually do it, run like
#  DRYRUN=1 ./02_runEntropyPSC.sh
[ -n "${DRYRUN:-}" ] && DRYRUN="echo" || DRYRUN=
# run from script directory
cd -P "$(dirname "$0")"

case $(hostname) in
   *rhea*) 
      ALLFILES=(/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/preprocessed_data/MGS/AfterWhole/ICAwholeClean_homogenize/*icapru.set);;
   *)
      ALLFILES=(data/ICAwholeClean_homogenize/*.set) ;;
esac
outdir=$PWD/Results/MGS_Entropy/individual_subject_files
log_file=/ocean/projects/soc230004p/shared/SignalComplexityAcrossAdolescene/log/%j-%x.log # job number and job name

echo "# [$(date)] collecting set file list"
echo "# [$(date)] running for ${#ALLFILES[@]} set files"

for inputfile in "${ALLFILES[@]}"; do 
   [ ! -r "$inputfile" ] &&
      echo "cannot read $inputfile" >&2 && continue

  ld8=$(grep -oP "\d{5}_\d{8}" <<< "$inputfile")

  for EPOCH in "${EPOCHS[@]}"; do

     # epoch determins what lengths there are
     case $EPOCH in 
       fix)   lengths=(0);
	       # TODO: this is unused and should maybe stay that way (20241206WF)
	       sbatch_args=(--ntasks-per-node=1 -p RM-shared --time=06:00:00);;
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
         $DRYRUN sbatch -J "$job_name" \
		 -e "$log_file" -o "$log_file" \
                 --export="INPUTFILE=$inputfile,EPOCH=$EPOCH,LENGTH=$LENGTH" \
                 $PWD/sbatch_entropy.bash
	         # TODO add "${sbatch_args[@]}" to force fix to one core on RM-shared (12Gb ram)

         [ -n "${ONLYONE:-}" ] && break 3
     done
  done
done

echo "# [$(date)] finished submitting jobs. $0 done"
