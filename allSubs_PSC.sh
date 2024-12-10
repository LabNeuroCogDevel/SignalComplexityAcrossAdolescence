for inputfile in data/ICAwholeClean_homogenize/*.set; do 
	ld8=$(grep -P "\d{5}_\d{8}" <<< "$inputfile")
	test -r Results/MGS_Entropy/individual_subject_files/${ld8}_MultiScaleEntropy_eyesClosed.csv && continue
	# 20241116WF - skip an actual output file (delay)?
	#test -r Results/MGS_Entropy/individual_subject_files/${ld8}_MultiScaleEntropy_delay8.csv && continue
	sbatch -J $ld8  --export=INPUTFILE=$inputfile entropy_onSubject_PSC.sh


done

