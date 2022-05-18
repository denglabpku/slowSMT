function [bootstat] = ResidenceTimesMultipleReps_Bootstrap(ResidenceFrames, Exposure, doPhotoBleach, H2B_bs_SPT_merged_SlowRate_mean)
%ResidenceTimesMultipleReps_Bootstrap 
%The goal of this function is to analyze residence times for multiple replicates of a single experiment whithin one folder. 

    %-----HISTORY -----% 
    % Originally written by Anders Sejr Hansen, Feb 2015
    % Modified by Zuhui Wang 2020/04/07, convert this funcion to be compatible with Script_ResidenceTimesMultipleReps_MergedData_Normalization.
    
    if nargin <3
        doPhotoBleach = 0;
        H2B_bs_SPT_merged_SlowRate_mean = {};
    end

    StartTimeForFit = [Exposure*2,Exposure*3, Exposure*4, Exposure*5, Exposure*6, Exposure*7, Exposure*8, Exposure*9, Exposure*10];
    
    %Store results for merged data
    SPT_merged_SlowFraction = zeros(length(StartTimeForFit),1);  % zeros: Create array of all zeros
    SPT_merged_FastRate = zeros(length(StartTimeForFit),1);
    SPT_merged_SlowRate = zeros(length(StartTimeForFit),1);
%     SPT_merged_OneRate = zeros(length(StartTimeForFit),1);
%     SPT_merged_OneExpA = zeros(length(StartTimeForFit),1);
    SPT_merged_TwoExpA = zeros(length(StartTimeForFit),1);

    SPT_merged_FractionBound = struct([]); 

    %% Repeat Res-Time fitting for MERGED trackedPar
    for k=1:length(StartTimeForFit)
        StartFrameForFit = round(StartTimeForFit(1,k) / Exposure);
        [SlowFraction, FastRate, SlowRate, TwoExpPrefixA, TwoExpFractionBound, FittingTrajNum] = ResTime2ExpFitter_v3_bootstrap(ResidenceFrames, StartFrameForFit, Exposure);
%         [OneRate, OneExpPrefixA, OneExpFractionBound] = ResTime1ExpFitter(SPT_merged_trackedPar, StartFrameForFit, Exposure, MinFrames);
        %store the data
        SPT_merged_SlowFraction(k,1) = SlowFraction;
        SPT_merged_FastRate(k,1) = FastRate;
        SPT_merged_SlowRate(k,1) = SlowRate;
%         SPT_merged_OneRate(k,1) = OneRate;
%         SPT_merged_OneExpA(k,1) = OneExpPrefixA;
        SPT_merged_TwoExpA(k,1) = TwoExpPrefixA;
%         SPT_merged_FractionBound(n).OneExpFractionBound = OneExpFractionBound;
        SPT_merged_FractionBound(1).TwoExpFractionBound = TwoExpFractionBound;   
    end

    %% Calcualte Residence Time for MERGED or MEAN (this is the output; the default is to use MERGED value)
    % WD Aug 2018    
    % Calculate MERGED residence time
    if doPhotoBleach == 0
        merged_ResTime2Exp_uncorrected= 1./SPT_merged_SlowRate;
    else
        merged_ResTime2Exp_uncorrected= 1./(SPT_merged_SlowRate-cell2mat(H2B_bs_SPT_merged_SlowRate_mean));
    end%if
    % merged_ResTime1Exp_uncorrected = 1/SPT_merged_OneRate(StartTimeForFit);

    % add an output structure containing multiple arguments; add Sep. 2021 WD
      bootstat = struct;
      bootstat.merged_ResTime2Exp_uncorrected = merged_ResTime2Exp_uncorrected;
      bootstat.SPT_merged_SlowFraction = SPT_merged_SlowFraction;
      bootstat.SPT_merged_FastRate = SPT_merged_FastRate;
      bootstat.SPT_merged_SlowRate = SPT_merged_SlowRate;       
end
