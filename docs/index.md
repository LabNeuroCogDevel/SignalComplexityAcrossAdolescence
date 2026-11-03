<img width="468" height="19" alt="image" src="https://github.com/user-attachments/assets/a863d67b-fb61-4a64-a717-8b009f647bed" /><br>

# Intrinsic Neural Timescales via the Autocorrelation Window Throughout Adolescence 

### Project Lead
Shane D. McKeon

### Faculty Lead
Beatriz Luna 

### Project Start Date
July 2024

### Current Project Status
Under Review

### Datasets
LNCD 7T

### Github Repository
[https://github.com/LabNeuroCogDevel/SignalComplexityAcrossAdolescence](https://github.com/LabNeuroCogDevel/SignalComplexityAcrossAdolescence)

## Code Documentation
**Preprocessing:**

Preprocessing can be run using [01_Entropy_Preprocessing.sh](https://github.com/LabNeuroCogDevel/SignalComplexityAcrossAdolescence/blob/main/01_Entropy_Preprocessing.sh)

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

### EEG Source Localization 
Preprocessed MGS and resting state EEG data was source localized using Brainstorm (version 3; Tadel et al., 2011) using unique scripts for each task state found in [the source localization folder](https://github.com/LabNeuroCogDevel/SignalComplexityAcrossAdolescence/tree/main/sourceLocalization). Each function reads in the corresponding epoch (delay, fixation, resting state), reads in the continuous EEG recordings and coregisters the EEG sensor locations to the participant-specific anatomy, and head-point registration was refined using a 2mm tolerance. . A three-layer boundary element model (BEM) consisting of the scalp, skull, and brain compartments was generated using Brainstorm. The EEG forward model was calculated using the OpenMEEG BEM method, with 1082 vertices the scalp, skull and brain, with tissue conductivities of 1 for the scalp, 0.0125 for the skull, and 1 for the brain.

Source reconstruction was performed on the cortical surface using a linearly constrained minimum variance (LCMV) beamformer implemented in Brainstorm. Data and noise covariance matrices were calculated from the epoched EEG recordings. The noise covariance was estimated over the 0–300ms time window, with the data covariance estimated over the 0–400ms time window. Covariance estimation was performed block-wise to reduce the influence of slow shifts in the data. The LCMV beamformer used a neural activity index (NAI) as the source measure, with free source orientation and a loose orientation constraint of 0.2. Depth weighting was applied to reduce bias toward superficial cortical sources, with a weighting exponent of 0.5 and a maximum weighting limit of 10. Noise regularization was set to 0.1, and the source kernel was computed for the EEG sensor configuration.

The resulting source estimates were parcellated according to the Human Connectome Project multimodal parcellation (HCP-MMP; Glasser et al., 2016), which divides the cerebral cortex into 360 regions (180 per hemisphere). Regional source time series were extracted by averaging source activity across all vertices within each Glasser parcel. Source orientations were not sign-flipped, and the extracted regional time series were not additionally normalized. These regional cortical time series were subsequently used for calculation of intrinsic neural timescales and spectral measures.

### Calculate Autocorrelation Window (ACW)
Source-localized resting-state EEG data were parcellated according to the 360-region HCP multimodal cortical parcellation, yielding one time series per cortical region for each epoch (resting state, fixation, and all three delay period lengths). The time series were extracted using [extractTimeSeries_beamforming.m](https://github.com/LabNeuroCogDevel/SignalComplexityAcrossAdolescence/blob/main/sourceLocalization/extractTimeSeries_beamforming.m). 

The autocorrelation window (ACW) was then derived from the autocorrelation function (ACF) within [ACW_glasser.m](https://github.com/LabNeuroCogDevel/SignalComplexityAcrossAdolescence/blob/main/ACW_glasser.m), which quantifies the similarity of a signal with itself at different time lags. The time lag when the ACF delays to 50% of its maximum was calculated and deemed ACW-50 (as this is the only ACW measure we report in this study; we will further refer to it as simply the ACW). This represents the temporal "half-life" of the signal's self-similarity and has been previous reported in the literature. A threshold of 2 standard deviations from the mean was used to exclude statistical outliers. The ACW for every trial is calculated and then averaged together, and are saved out per subject. The individual data frames are then read and combined using [CombineSubjectDataframes_Glasser.R](https://github.com/LabNeuroCogDevel/SignalComplexityAcrossAdolescence/blob/main/CombineSubjectDataframes_Glasser.R).

### FOOOF/Specparam
Spectral parameterization was performed on source-localized resting-state EEG data using the FOOOF/specparam algorithm implemented in Brainstorm here [brainstormSourceLocal_beamforming_fooof.m](https://github.com/LabNeuroCogDevel/SignalComplexityAcrossAdolescence/blob/main/sourceLocalization/brainstormSourceLocal_beamforming_fooof.m). Source time series were obtained for each cortical region of the Human Connectome Project (HCP) multimodal parcellation and analyzed independently. For each region, power spectral density (PSD) was estimated using Welch's method with 2-s windows and 50% overlap. Power was expressed in physical units (U^2/Hz).
The resulting power spectra were parameterized over the 1–60 Hz frequency range using the MATLAB implementation of specparam. The model decomposed each power spectrum into an aperiodic component and superimposed periodic oscillatory peaks. The aperiodic component was modeled using a fixed aperiodic mode, and the periodic component was characterized using up to three oscillatory peaks with peak widths constrained to 0.5–12 Hz. A minimum peak height of 3 dB was used for peak detection, with a proximity threshold of 2. Sixty-Hz power-line activity was excluded from the parameterization. For each cortical region, the fitted aperiodic component was characterized by its aperiodic exponent and offset. 

Spectral parameterization was performed independently for each of the 360 cortical regions of the HCP multimodal parcellation (180 per hemisphere). Regional aperiodic parameters were subsequently used to characterize the spectral properties of cortical activity and to examine their relationship with intrinsic neural timescales.





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

