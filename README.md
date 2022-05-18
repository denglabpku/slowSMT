# slowSMT
Calculating residence time etc using slowSMT trajectory data
------------------------------------------------------------
This repository contains Matlab code for performing the survival analysis and calculating residence time using slowSMT trajectory data.

## Overview
Using long exposure time (500 ms in our case), freely-diffusing molecules will be motion-blurred into the background that can not be detected by MTT localization algorithm. Then, the leftover motion-less molecules that detected by MTT algorithm will be tracked, and the recorded trajectory length of each 'bound' molecules can be used to generate a survival curve (1-CDF). The generated curve can then be fitted by kinetic models such as exponential model to extract the unbinding rate constant or k<sub>off</sub> <br><br>

## SurvivalCurve
Under this folder, run 'script_SurvivalAnalysis.m' to generate the bootstrapped survival curve from the slowSMT trajectory data. The input of this script is a MAT-file that contains 'trackedPar' matlab variable that generated from MTT localization and tracking software (https://gitlab.com/tjian-darzacq-lab/SPT_LocAndTrack/ParallelProcess_slowSPT_JF646.m). In the script, we implement the matlab built-in _bootstrp_ function to perform bootstrap on the whole trajectory lengths calculated from 'trackedPar' when generating survival curve. For details of how bootstrp works, please check the matlab help documents (https://www.mathworks.com/help/stats/bootstrp.html?s_tid=srchtitle_bootstrp_1).<br>

The main output of the script is 'trackedEnsemble', which contains the major info and value of the survival curve, you can find its definition in the description section of the script.

Dependent functions:
1) shadedErrorBar.m <br>
Creates an attractive shaded error region rather than discrete bars.<br>
You can download it from matlab file exchange (https://www.mathworks.com/matlabcentral/fileexchange/26311-raacampbell-shadederrorbar?s_tid=srchtitle_shadedErrorbar_1) or from Github https://github.com/raacampbell/shadedErrorBar.
2) wrapped_ksdensity.m <br>
A wrapped kernal smooth function to estimate the y value of the survival curve.<br><br>


## BootstrapResidenceTime
Under this folder, run 'script_ResidenceTimesMultipleReps_MergedData_Normalization.m' to perform bootstrapped 2-exponential fitting on survival curve of the slowSMT trajectory data. Different from the script in <b>SurvivalCurve</b>, this script will not plot survival curve directly but will perform fitting on the survival curve generated from trajectory lengths calculated from 'trackedPar'. The script implements the built-in matlab _bootstrp_ function to bootstrap on trajectory lengths and perform survival curve fitting. The user is also allowed to perform photobleach correction by providing the slowSMT data of H2B. <br><br>

Another feature of this script is that it allows the user to start curve fitting only on trajectory lengths longer than the defined minimal lengths (see the definition of _MinFrames_ in the script). This feature allows user to calculate residence time using the longer trajectory lengths if short trajectory lengths in the slowSMT dataset is a concern. <br><br>

The main output of the script - _trackedEnsemble_ - contains below fields:
* bs_ResidenceTime : specific residence time calculated from bootstrapped data
* bs_nsResidenceTime : nonspecific residence time calculated from bootstrapped data
* bs_SlowFraction : fraction of specific residence time portion in 2-exponential fitting
* bs_SlowRate : slow unbinding rate, reciprocal of the specific residence time
* bs_FastRate : fast unbinding rate, reciprocal of the nonspecific residence time

Another script 'script_ResidenceTimesMultipleReps_PlotOnly.m' provides some useful codes for plotting figures, which is not a part of the analysis.

For detailed explanation on calculating residence time and survival analysis using slowSMT data, we recommend readers to go through the "Residence time measurements from SMT" parts in Material and Methods sections published by Hansen et al, 2017 (https://elifesciences.org/articles/25776). <br><br>

Citation:<br>
Hansen, Anders S., Iryna Pustova, Claudia Cattoglio, Robert Tjian, and Xavier Darzacq. “CTCF and Cohesin Regulate Chromatin Loop Stability with Distinct Dynamics.” ELife 6 (May 3, 2017): e25776.


