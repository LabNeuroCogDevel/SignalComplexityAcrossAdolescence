for inputfile in data/ICAwholeClean_homogenize/*.set; do 
	ld8=$(grep -P "\d{5}_\d{8}" <<< "$inputfile")
	finalFile=Results/Rest_Entropy/individual_subject_files/${ld8}_MultiScaleEntropy_eyesClosed.csv 
	test -r $finalFile && echo found $_ &&  continue
	echo do not have $finalFile
	# 20241116WF - skip an actual output file (delay)?
	#test -r Results/MGS_Entropy/individual_subject_files/${ld8}_MultiScaleEntropy_delay8.csv && continue
	dryrun sbatch -J $ld8  --export=INPUTFILE=$inputfile entropy_onSubject_PSC.sh


done

