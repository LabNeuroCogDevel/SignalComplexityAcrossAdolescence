#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=2
module load matlab
matlab -nodisplay -r 'try, matlab_01_entry(4), catch(e), e, end, quit'
