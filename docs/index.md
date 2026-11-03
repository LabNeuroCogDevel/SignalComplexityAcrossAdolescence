
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

### Segmenting EEG data
MGS EEG data was segmented into the different delay periods (6, 8 and 10 seconds) using the following script [selectDelayPeriods.m](https://github.com/LabNeuroCogDevel/SignalComplexityAcrossAdolescence/blob/main/selectDelayPeriods.m)

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

### Associations between working memory vs resting state
Given prior work suggesting that intrinsic neural timescales (INTs) are modulated by task state, we obtained ACWs for 4 minutes rest and during a memory guided saccade working memory task (MGS). We first assessed whether ACW differed across the three delay lengths (6, 8, and 10 seconds). As no significant differences in ACW were observed across delay lengths, ACW values were averaged across the three delays within each Glasser region, yielding a single delay-epoch ACW measure per region for each participant. We next tested whether ACW differed between the delay and fixation epochs of the MGS task. As no significant differences were observed between these epochs, delay and fixation ACW values were averaged within each region to derive a unified task-state ACW measure for comparison with resting-state ACW. This work can be found in [INT_by_task_glasser.Rmd](https://github.com/LabNeuroCogDevel/SignalComplexityAcrossAdolescence/blob/main/INT_by_task_glasser.Rmd)

### Results: Intrinsic neural timescales decrease through adolescence to adulthood in rest and task states. 
We examine the developmental trajectory of intrinsic neural timescales (INT) using EEG autocorrelation window (ACW) data across different brain regions and task conditions in [INT_developmental_effects_glasser.Rmd](https://github.com/LabNeuroCogDevel/SignalComplexityAcrossAdolescence/blob/main/INT_developmental_effects_glasser.Rmd) <br>

***Developmental Analyses***

To assess developmental changes in the ACW across adolescence, we employed semiparametric GAMMs (generalized additive mixed models) with penalized splines using the mgcv package in R. Independent GAMMs were run for each Glasser region parcel. The developmental models were fit with ACW as the dependent variable modelling age as a smooth term. Fixed effects for sex were also included in the model, as well as a random intercept per participant. To prevent overfitting, a maximum basis complexity (k) of 3 was used. Significance values were FDR corrected (n=360). In addition to the main development models, independent GAMMs per region with an age-by-task interaction term were fit to assess whether age-related patterns differed across task states. These GAMMs were fit with ACW as the dependent variable, smooth terms for age, a factor-smooth interaction term between age and task, as well as a random intercept per participant. Fixed effects of sex were also included in the model.
We then sought to assess the rate of ACW change with age in each of region. A GAMM (generalized additive mixed model) model was fit with the ACW as the dependent variable, age as a smooth term, controlling for sex, and with a random intercept per participant. The first derivative of the age spline was computed to estimate the rate of developmental change using the gratia package in R, using 200 equally spaced age intervals. To index the average rate of ACW change across the age spline, the average first derivative was computed, where a more negative average derivative is indicative of greater maturational decline.

***Associations with the S-A Axis***

Parcel-wise ACW values were averaged by age group (adolescent, adult) and cognitive state (rest, task) and correlated with ranked sensorimotor–association (S–A) axis values from the Glasser 360-parcel atlas. Pearson correlations assessed spatial gradients in ACW across the cortical hierarchy. To examine developmental organization, mean derivatives of parcel-wise age-related ACW trajectories from GAMs were similarly correlated with S–A axis values separately for rest and task. Parcels failing established signal-to-noise criteria were excluded. Linear regression lines were used for visualization.

***Task- versus Resting-State ACW***

To examine whether intrinsic neural timescales differed between cognitive states, we compared ACW between task and resting-state epochs within each Glasser parcel. Analyses were conducted separately for adolescents and adults to determine whether state-dependent differences in ACW varied across developmental groups. For each Glasser parcel and age group, task and resting-state ACW values were paired within participant and compared using paired-samples t-tests. The resulting t-statistics were used to characterize the magnitude and direction of task–rest differences across the cortex, with positive t-values indicating greater ACW during task relative to rest and negative t-values indicating greater ACW during rest relative to task. To account for multiple comparisons across cortical parcels, p-values were corrected using the Bonferroni method, with corrected p < 0.05 considered statistically significant.

### Results: Shorter intrinsic neural timescales are associated with EEG-derived measures of E/I balance and aggregate spiking activity. 
To examine the relationship between intrinsic neural timescales and aperiodic spectral properties, we assessed associations between ACW and the aperiodic exponent and offset derived from resting-state EEG using [ACW_fooof_glasser.R](https://github.com/LabNeuroCogDevel/SignalComplexityAcrossAdolescence/blob/main/code/ACW_fooof_glasser.R). Aperiodic measures were first averaged across available estimates within each participant and Glasser parcel and subsequently aligned with the corresponding parcel-wise ACW values. Analyses were restricted to the resting-state epoch. Associations between ACW and each aperiodic measure were assessed separately within each Glasser parcel using linear mixed-effects models. Aperiodic exponent or offset was entered as the dependent variable, ACW as the predictor of interest, and age as a covariate, with participant included as a random intercept to account for repeated observations. Separate models were performed for the aperiodic exponent and offset across cortical parcels. To visualize these relationships in representative cortical regions, bilateral parcel values were averaged within the superior temporal sulcus (STSvp), visual cortex (V4), primary motor cortex (area 4), and dorsolateral prefrontal cortex (area 46) for each participant and session. Linear mixed-effects models of the same form were then used to assess ACW–aperiodic relationships within each representative region, controlling for age and including participant as a random intercept. Corresponding analyses were performed for both the aperiodic exponent and offset.

### Results: Shorter intrinsic neural timescales are associated with improved working memory performance in adolescence.
During EEG data collection, participants completed a memory guided saccade task. The primary behavioral outputs from the MGS task were the average accuracy, accuracy variability, latency, and latency variability. To assess the effects of the ACW within each region on working memory performance, we used GAM models, independently run for each performance metric, here: [behavioralEffects_glasser.R](https://github.com/LabNeuroCogDevel/SignalComplexityAcrossAdolescence/blob/main/behavioralEffects_glasser.R). We first tested for significant main effects, with the behavioral measure as the dependent variable, age as a smooth term, fixed effects of ACW, task state (resting or task), and sex. Participant specific random intercepts were modeled using a random effect term to account for repeated measures. In addition to the main effect models, independent GAMs per region with an ACW-by-task interaction term were fit to assess whether ACW-related patterns on performance differed across task states. These GAMs were fit independently with each behavioral measure as the dependent variable, smooth terms for age, an interaction smooth term between ACW and task, as well as a random intercept per participant. Fixed effects of sex and task state were also included in the model. For each behavioral measure, p-values were Bonferroni corrected based on the 3 cortical regions of interest.

Additional analyses were conducted on the associations that yielded significant results to further interrogate how the relationship between ACW and behavior changes across development. To do so we modeled the age-varying effect of ACW on each behavioral measure using a semiparametric GAM framework. This is conceptually equivalent to a time-varying effect model (TVEM), allowing the strength of the ACW-behavior relationship to change across age. Specifically, we modeled behavior as a function of smooth age terms, with an interaction smooth term between age and ACW. The model also included fixed effects for task state and sex, as well as a participant-specific random intercept. To characterize the age-varying effect of the ACW-behavior association, predicted values were visualized using ggpredict in R. The first derivative of the age smooth was computed to find the significant points of developmental change, identified as the age ranges where the derivative confidence interval did not cross zero. 


### Results: Associations between intrinsic neural timescales and cortical myelination. 
MRI data acquisition and preprocessing have been previously described in Sydnor et al., 2025. The corresponding github can be found here [cortical myelination github](https://github.com/LabNeuroCogDevel/corticalmyelin_maturation)

To study the associations between cortical R1 and EEG-derived ACW we leveraged developmental GAMs with R1 (superficial or deep) as the dependent variable, age as a smooth term, ACW as a smooth term, and a random intercept per participant ([R1_ACW50_glasser.R](https://github.com/LabNeuroCogDevel/SignalComplexityAcrossAdolescence/blob/main/R1_ACW50_glasser.R)). GAMMs were fit separately for the 3 cortical regions (frontal, parietal, and occipital) for task and resting state ACW values. To assess the strength and significance of the associations between R1 and the ACW in superficial and deep layers, extracted the F-statistic and corresponding p-value for the ACW smooth term from each model. For task and resting state conditions, p-values were Bonferroni corrected based on the 3 cortical regions of interest. 

We then sought to interrogate for differences in the strength of the R1-ACW relationship by depth. Thus, we fit GAMMs with a depth interaction term, separately for the task and resting state conditions, as well as, for frontal, parietal, and occipital regions. The models included R1 (superficial or deep) as the dependent variable, age as a smooth term, and an ordered factor smooth interaction between depth and ACW. Depth was also included as a linear covariate along with a random intercept per participant. 


### Results: Developmental decreases in intrinsic neural timescales replicate in human iEEG 

***Participants***

iEEG recordings were performed in 14 neurosurgical patients with drug-resistant epilepsy as part of clinical evaluation for epilepsy surgery. 

***Acquisition and preprocessing***

Each patient underwent implantation of depth electrodes as previously described Abel et. al 2018. Recording sites were selected based on clinical necessity for seizure localization. Briefly, Dixi Medical Microdeep electrodes were employed, featuring a diameter of 0.8 mm, a contact length of 2 mm, and a center-to-center spacing of 3.5 mm. Each electrode comprised between 8 and 18 contacts, which are referred to as channels. Neural recordings were collected during 5-minute eyes open resting state, recorded at 1kHz using the Ripple Grapevine Nomad system (model R02000) with real-time notch filtering applied at 60, 120, and 180 Hz. Data was resampled to 512 Hz and a common average reference (CAR) filter was used to remove noise shared across channels. 

***Autocorrelation Window Calculation***

The autocorrelation window (ACW) was calculated using the same approach described previously ([sEEG_ACW.m](https://github.com/LabNeuroCogDevel/SignalComplexityAcrossAdolescence/blob/main/sEEG_ACW.m)). In brief, resting-state recordings from each channel were segmented into 2-second epochs. For each epoch, the autocorrelation function (ACF)—which quantifies the similarity of a signal with itself across varying time lags—was computed. The ACW was defined as the time lag at which the ACF decayed to 50% of its maximum value, a measure referred to as ACW-50. As this is the only ACW metric reported in the current study, it will be referred to hereafter simply as the ACW. This metric represents the temporal "half-life" of the signal’s self-similarity and has been previously described in the literature. Only channels in non-seizure sites, determined by clinical evaluation, were used in the main analysis. A threshold of 2 standard deviations from the mean was used to exclude statistical outliers.

***Developmental analyses***

As in our longitudinal developmental cohort, we assessed age-related changes in ACW across adolescence using semiparametric GAMMs (generalized additive mixed models) with penalized splines using the mgcv package in R ([sEEG_ACW.Rmd](https://github.com/LabNeuroCogDevel/SignalComplexityAcrossAdolescence/blob/main/sEEG_ACW.Rmd)). Independent GAMMs were run for each regional analysis (frontal, parietal, and occipital). The developmental models were fit with ACW as the dependent variable, age as a smooth term. Fixed effects for a random intercept per participant. To prevent overfitting, a maximum basis complexity (k) of 3 was used. Significance values were Bonferroni corrected (n=3). In addition to the main development models, independent GAMMs per region with an age-by-task interaction term were fit to assess whether age-related patterns differed across task states. These GAMMs were fit with ACW as the dependent variable, smooth terms for age, an interaction smooth term between age and task, as well as a random intercept per participant. 


### Results: Developmental decreases in intrinsic neural timescales replicate in macaque LFPs

***Subjects***

Neurophysiological data were collected from 8 rhesus monkeys (Macaca mulatta), comprised of 6 males and 2 females. All surgical procedures and animal use protocols were reviewed and approved by the Institutional Animal Care and Use Committees at Wake Forest University and Vanderbilt University, in compliance with the U.S. Public Health Service Policy on Humane Care and Use of Laboratory Animals and the National Research Council’s Guide for the Care and Use of Laboratory Animals. 

***Autocorrelation Window Calculation***

The autocorrelation window (ACW) was calculated using the same approach described previously ([ACW_MacaqueLFPs.m](https://github.com/LabNeuroCogDevel/SignalComplexityAcrossAdolescence/blob/main/ACW_MacaqueLFPs.m)). In brief, recordings from each channel, for each behavioral session, for every trial, for each monkey, were separated into the delay and fixation epochs of the ODR task. For each epoch, the autocorrelation function (ACF)—which quantifies the similarity of a signal with itself across varying time lags—was computed. The ACW was defined as above, as the time lag at which the ACF decayed to 50% of its maximum value, a measure referred to as ACW-50. As this is the only ACW metric reported in the current study, it will be referred to hereafter simply as the ACW. This metric represents the temporal "half-life" of the signal’s self-similarity and has been previously described in the literature20,22. A threshold of 2 standard deviations from the mean was used to exclude statistical outliers.

***Developmental modeling***
 
To examine developmental changes in the autocorrelation window during task performance, we fit a generalized additive mixed model (GAMM) using the bam() function (gam for large datasets) from the mgcv package in R ([monkeyDevelopmentalEffect.R](https://github.com/LabNeuroCogDevel/SignalComplexityAcrossAdolescence/blob/main/monkeyDevelopmentalEffect.R)). The model included smooth terms for age nested within individual monkeys to capture subject-specific developmental trajectories, as well as separate smooths for behavioral session, trial number, and channel-by-monkey interactions to account for within-session and within-subject variability. Additional random effects were included for neuron recording area, monkey identity, and recording epoch (delay or fix). A smooth interaction between age and epoch was also included to model age-related effects that varied across experimental epochs. 

To visualize developmental trajectories of ACW across adolescence, we generated smooth fits using generalized additive mixed models (GAMMs) with ggplot. Individual trajectories were plotted for each monkey using subject-specific smooths with fixed basis complexity (k = 3), providing within-subject developmental trends. A group-level trend line was overlaid using a GAM with an unfixed smoothing parameter, capturing the overall developmental pattern across the cohort. Line color indicated monkey identity, while line type distinguished between experimental epochs. Shaded confidence intervals were included for both individual and group fits.



