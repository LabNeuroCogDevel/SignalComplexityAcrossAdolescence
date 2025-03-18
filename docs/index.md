<br>

# Intrinsic Neural Timescales via the Autocorrelation Window Throughout Adolescence 

### Project Lead
Shane D. McKeon

### Faculty Lead
Beatriz Luna 

### Project Start Date
July 2024

### Current Project Status
In progress 

### Datasets
LNCD 7T

### Github Repository
[https://github.com/LabNeuroCogDevel/7T_EEG/tree/main/Entropy](https://github.com/LabNeuroCogDevel/SignalComplexityAcrossAdolescence/tree/main)

## Code Documentation
**Preprocessing:**

Preprocessing can be run using [01_Entropy_Preprocessing.sh](/LabNeuroCogDevel/7T_EEG/blob/main/Entropy/01_Entropy_Preprocessing.sh)

Note this is the same preprocessing as the resting state data in the [Aperiodic EEG Project](https://labneurocogdevel.github.io/7T_EEG/fooofMRS.html)

  ```matlab -nodesktop -r "addpath(genpath('../Preprocessing_Functions/')); run_preprocessing_pipeline('Resting_State')" ```
  
* Initial preprocessing was done using matlab code run_preprocessing_pipeline.m (../Preprocessing_Functions) which first pulls in raw data from 'Raw/EEG/7TBrainMech'.
* Set the task as 'Resting_State' to select the resting state data.
* Bandpass filter between 0.5 Hz and 70 Hz
* Downsamples the data from 1024Hz to 150 Hz
* Removes bad channels (the following criterion is used)
  - arg_flatline: 8
    - Maximum tolerated flatline duration. In seconds. If a channel has a longer flatline than this, it will be considered abnormal.
  - arg_highpass: [0.25 0.75]
    - Transition band for the initial high-pass filter in Hz. This is formatted as [transition-start, transition-end]
  - arg_channel: 0.7
    - Minimum channel correlation. If a channel is correlated at less than this value to a reconstruction of it based on other channels, it is considered abnormal in the given time window. This method requires that channel locations are available and roughly correct; otherwise a fallback criterion will be used.
  - arg_noisy: 5
    - If a channel has more line noise relative to its signal than this value, in standard deviations based on the total channel population, it is considered abnormal.
  - arg_burst: 15
    - Standard deviation cutoff for removal of bursts (via ASR). Data portions whose variance is larger than this threshold relative to the calibration data are considered missing data and will be removed. 
  - arg_window: 0.3
    - Criterion for removing time windows that were not repaired completely. This may happen if the artifact in a window was composed of too many simultaneous uncorrelated sources (for example, extreme movements such as jumps). This is the maximum fraction of contaminated channels that are tolerated in the final output data for each considered window.
* Interpolates missing channels
  - Dataset includes a few subjects that used a 128 cap as opposed to a 64 channel cap. The code removes the 4 channels that are found in 128 but not 64 and reinterpolates the missing channels that were removed from above
* Run ICA to identify eye movements and blinks
* Homogenize Channel locations
  - Read in the channel locations and make sure all files have the correct locations, especially the few subjects who were ran using a 128 channel cap
* Filter out 60 Hz artifact from line noise

<br> 

### Calculate Autocorrelation Window (ACW)
Preprocessed MGS and resting state EEG data was used to calculate the autocorrelation window during the delay and fixation epochs, as well as, the eyes open and eyes closed segments of resting state using [signalComplexityCalculations.m](https://github.com/LabNeuroCogDevel/SignalComplexityAcrossAdolescence/blob/main/signalComplexityCalculations.m). 
This function read in what task (resting state or MGS), epoch (delay, fix, restEyesOpen, or restEyesClosed), lengthValue (how long the epoch was, really only for delay and length 6,8, or 10), segmentValue (run on the entire epoch (0) or in 2 second segments (2)), and calculation (ACW or MSE). For the purpose of this project, a segment length of 2 was set. The function then sets the corresponding datapath and path of where to the save the results, loads in the EEG, filters for line noise, and then correctly epochs the data. It will then call [calculate_ACW.m](https://github.com/LabNeuroCogDevel/SignalComplexityAcrossAdolescence/blob/main/code/calculate_ACW.m) which runs acw on every trial, for every channel, broken into 2 second chunks (if desired) and saves out the results for every subject. The individual data frames are then read and combined using [CombineSubjectDataframes.R](https://github.com/LabNeuroCogDevel/SignalComplexityAcrossAdolescence/blob/main/CombineSubjectDataframes.R) 

### Results: Increases in ACW Emerge in WM Delay Compared to Resting State in Parietal and Occipital Regions 
**Associations Between MGS Delay Period ACWs and Eyes Open Resting State** <br>
Intrinsic Neural Timescales (INTs) via the ACW using EEG data across different cognitive states (rest vs. task-related delays) and brain regions were assessed in [INT_by_task.Rmd](https://github.com/LabNeuroCogDevel/SignalComplexityAcrossAdolescence/blob/main/INT_by_task.Rmd). Briefly, this markdown processes and visualizes autocorrelation window (ACW) data from whole-brain and regional EEG recordings, comparing INT between a resting state and a delay task. Statistical tests (paired t-tests) assess differences in INT across conditions, both globally and within specific regions (frontal, parietal, occipital). The script generates multiple figures, including boxplots and EEG electrode maps, to illustrate findings across age groups (adolescents vs. adults). Results are saved as PDFs for further analysis and presentation.

### Results: Maturation of ACW in Frontal, Parietal, and Occipital Regions During WM Delay and Rest 
We examine the developmental trajectory of intrinsic neural timescales (INT) using EEG autocorrelation window (ACW) data across different brain regions and task conditions in [developmentalEffects.Rmd](https://github.com/LabNeuroCogDevel/SignalComplexityAcrossAdolescence/blob/main/developmentalEffects.Rmd) <br>

**Developmental decreases in the autocorrelation window** <br>
It visualizes ACW changes with age in frontal, parietal, and occipital regions during rest and delay periods, employing generalized additive models (GAMs) to analyze maturation trends. Statistical analyses assess age-related INT differences, controlling for sex and subject variability.

**Rate and Age of Maturation** <br>
We next calculate and visualize the rate and age of maturation for the ACW across adolescence in the frontal, parietal, and occipital regions. By iterating over each epoch and region, the script applies the growthrate() (found in [code/rateOfMaturation.R](https://github.com/LabNeuroCogDevel/SignalComplexityAcrossAdolescence/blob/main/code/rateOfMaturation.R) function to estimate GAM-based statistics, including the first derivative of ACW with respect to age, which represents the rate of maturation. This approach allows for identifying periods of significant ACW change across adolescence, providing insights into the developmental trajectories of intrinsic neural timescales. The results are stored in lists for statistical values, derivative estimates, and fitted GAM models, facilitating further analysis of brain maturation patterns. In this specific analysis, the age of maturation is derived by calculating the smooth decrease offset from the GAM model, aka the youngest age with a significant negative derivative.

### Results: ACW Across Different WM Epochs vs Resting State
In this analysis, the focus was on examining the differences in ACW across various conditions and regions, particularly comparing task states (Delays and Fix) with rest (Eyes Open), as seen in [taskVsRest.Rmd](https://github.com/LabNeuroCogDevel/SignalComplexityAcrossAdolescence/blob/main/taskVsRest.Rmd). The data were filtered to include only relevant epochs (Delays, restEyesOpen, and Fix) and analyzed using a series of statistical tests and visualizations. Boxplots were used to compare ACW values between conditions, with significant differences identified between task and rest in specific regions like Parietal and Occipital. The results showed that ACW values decrease with age across the task states, and no significant differences were found between the Delays and Fix epochs. Additionally, comparisons between task and rest in different age groups revealed significant differences in the Parietal and Occipital regions, but these differences were consistent across both adolescents and adults. The GAM models further supported these findings, indicating significant age-related changes and differences between task and rest conditions.


### Results: Decreases in the ACW Support Improved WM Performance 
**ACW Associations with WM performance measures in the frontal, parietal, and occipital regions** <br>
In order to assess behavioral relationships with the ACW, the following functions found here [code/behavioralStats.R](https://github.com/LabNeuroCogDevel/SignalComplexityAcrossAdolescence/blob/main/code/behavioralStats.R) were used. Each function applies the GAM model using the mgcv package, followed by extracting and formatting the summary results (fixed and smooth terms). 

  - behavioralInteraction: This function investigates the interaction between a brain measure (e.g., ACW) and a behavioral measure (e.g., task performance) across different epochs, considering age, sex, and participant-specific random effects. It fits a generalized additive model (GAM) to the data and returns a list containing the model, the fixed effects results, and the smooth terms.
  - behavioralAgeInteraction: Similar to the first function, this one also explores the interaction between brain and behavioral measures but with a focus on how the relationship varies with age. The model includes smooth terms for age and age-by-brain measure interaction, allowing for more flexibility in age-related effects.
  - behavioralRegionInteraction: This function extends the analysis to examine how brain-behavior interactions differ across regions. It incorporates region as a factor in the model to understand if the relationship between brain measures and behavioral performance varies by region, in addition to accounting for age, sex, and participant-specific random effects.
  - behavioralMainEffect: This function focuses on the main effect of a brain measure on a behavioral outcome, considering age, sex, and random effects, without including interactions. It serves to assess the direct relationship between the brain measure and the behavioral outcome.

For the main analysis, behavioralInteraction and behavioralMainEffect were used to assess the main effect of ACW on WM performance, as well as any epoch (task vs rest) by ACW interactions. We then further assessed the significant results with a Time-Varying Effect Model (TVEM) which used the behavioralAgeInteraction function. Predictions for WM performance are generated and visualized using ggpredict, with confidence intervals included. The derivatives of the model are then assessed to identify where the derivative is no longer significant, which highlights the points where the time-varying effect model (TVEM) transitions from being significant to non-significant. This analysis helps pinpoint the age ranges or brain measure values at which the relationship between age and behavioral performance becomes statistically meaningful. When the confidence intervals around the derivatives no longer exclude zero, it indicates that the relationship between age and behavior is no longer significantly changing, marking the shift from a significant to non-significant effect. This provides insight into the dynamics of the brain-behavior interaction across age.

