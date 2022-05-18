%% This is a script to do survival analysis of slow SPT data 
% Zuhui Wang 2020/04/09
% Wulan Version April 14th. 2021 Sep. 

% Purpose: Do survival analysis with the option of using bootstrapped samples and
% can draw shadow-like errorbar overlaid on the survival function.
% Original file name: Script_SurvivalAnalysis_multipleSampleOverlay_yValue_fitcurver.m

% Dependency:
%   -shadedErrorBar.m
%   -wrapped_ksdensity.m

% Important variable explanation:
%   trackedEnsemble.Filenames: name of input MTT tracked mat file
%   trackedEnsemble.ResidenceFrames: residence frames after filter out
%   minimal frame length
%   trackedEnsemble.ResidenceFramesGoneFreq: the proper spaced x value (residence frames) of 1-CDF
%   trackedEnsemble.cFrames: logical value to check if there is any censor
%   data in ResidenceFrames (>= total frame length)
%   trackedEnsemble.bs_ksdensity: 
%       trackedEnsemble.bs_ksdensity.f: 1-CDF value of bootstrapped data
%       trackedEnsemble.bs_ksdensity.xi: the proper spaced x value
%       (residence frames) of 1-CDF, same value as trackedEnsemble.ResidenceFramesGoneFreq
%       trackedEnsemble.bs_ksdensity.ResidenceTime: ResidenceFrame*Exposure
%       trackedEnsemble.bs_ksdensity.f_abovecutOff: the max possible yValue matrix (a vertically concatenated f
%       value returned from ksdensity from bootstapped sample)
%   trackedEnsemble.bs_ResidenceFramesIdx: bootstrapped ResidenceFrame
%   index of the input residence frame sample
%   trackedEnsemble.OneExpFit: 
%       trackedEnsemble.OneExpFit.fitObject: returned fit object of 1-exp
%       trackedEnsemble.ObeExpFit.yValue: evaluated y value of the 1-exp fit function at trackedEnsemble.bs_ksdensity.xi
%   trackedEnsemble.TwoExpFit: similar to trackedEnsemble.OneExpFit, except using 2-exp fit model

clc; clear all; close all;

%% shared parameters
Exposure = 0.5; %seconds
nboot = 20; %number of bootstrap
MinFrames = 2; %MinFrames: only consider trajectories with at least this many frames
yValue_cutOff = 10^(-5); % cutoff of 1-CDF, of which the 1-CDF below this value won't be plot to save errorbar plot from error caused by possible negative std. 
output_path = '/home/dell/Documents/METHOD/SPT/SPT_Analysis/denglabpku-repository/slowSMT-denglabpku/ExampleFigures/';
trackedEnsemble = struct;

%%  -- a section to add fit curves of H2B (iter == 1), or FOXA2 DE (iter == 3)---
plot_fit_curve = 0;

%% Input files to be analyzed
% % %---- 500ms ----%
% input_path = '/Volumes/AnyShare/MyDocuments.localized/邓伍兰_1906184199/2_课题/1_FOXA2_SPT/DATA_Analysis_sep2019_current/slowSPT_JF646_MaxD_007_loc_error_6pt25/hESC_trajectories_633laser4pert250mW_MaxD_007_loc_error_6pt25_500ms/trajectory_matlab_files_Merged/';
% mat_files = {'500ms_H2B_DE_Tracked_SPT_merged_trackedPar.mat', '500ms_KI14_APS_Tracked_SPT_merged_trackedPar.mat',  '500ms_KI14_DE_Tracked_SPT_merged_trackedPar.mat', '500ms_dC11_DE_Tracked_SPT_merged_trackedPar.mat', '500ms_dCdN10_DE_Tracked_SPT_merged_trackedPar.mat'};
% legend_label = {'H2B-Halo DE', 'FOXA2-Halo APS', 'FOXA2-Halo DE', '\DeltaCTD-Halo DE', '\DeltaC\DeltaN-Halo DE'};
% figure_name = '1-CDF_all_nboot20_minframe2_1';
% % color_plate = {[0.4118 0.4118 0.4118], [0.2588 0.5451 0.7922],[0.8510 0.3255 0.3098], [0.3608 0.7216 0.3608], [0.3922 0.1804 0.4863]}; %  black blue red gree purple 
% color_plate = {[0.1, 0.1,0.1], [0.8510 0.3255 0.3098], [0.3608 0.7216 0.3608], [0.5, 0, 0.5], [0.3569 0.7529 0.8706]}; 
% % match color with RTvsThreshold plot, black, red, green, purple, cyan

% %---- 500ms slowSMT pooled Track MAT----% %
input_path = '/mnt/f58069a5-1cf3-43b8-bb9b-ea74327327c9/WZH-DataCenter/Submission/2022_FOXA2Data/FOXA2_SMT_data/slowSMT_trajectories_pooled';
mat_files = {'DE_H2B-Halo_2Hz_trackedPar.mat', 'APS_FOXA2-Halo_C14_2Hz_trackedPar.mat',  'DE_FOXA2-Halo_C14_2Hz_trackedPar.mat', 'DE_dCTD-Halo_C11_2Hz_trackedPar.mat', 'DE_dCdN-Halo_C10_2Hz_trackedPar.mat'};
legend_label = {'H2B-Halo DE', 'FOXA2-Halo APS', 'FOXA2-Halo DE', '\DeltaCTD-Halo DE', '\DeltaC\DeltaN-Halo DE'};
figure_name = 'survival_1-CDF_all_nboot20_minframe2_1';
% color_plate = {[0.4118 0.4118 0.4118], [0.2588 0.5451 0.7922],[0.8510 0.3255 0.3098], [0.3608 0.7216 0.3608], [0.3922 0.1804 0.4863]}; %  black blue red gree purple 
color_plate = {[0.1, 0.1,0.1], [0.8510 0.3255 0.3098], [0.3608 0.7216 0.3608], [0.5, 0, 0.5], [0.3569 0.7529 0.8706]}; 
% match color with RTvsThreshold plot, black, red, green, purple, cyan


%% Data read and bootstrap

for iter = 1:length(mat_files)
    
    disp('bootstrapping sample trackedPar mat files....');
    disp(mat_files(iter));
    
    % Read the input track mat file generated from MTT
    file_name = char(mat_files(iter));
    trackedEnsemble(iter).Filenames = file_name(1:end-4);
    load([input_path filesep trackedEnsemble(iter).Filenames '.mat'],'trackedPar');
    iiter = 1;
    
    % Filter out the track length that are less than MinFrames
    for i=1:length(trackedPar)  
        %Account for missed frames/gap frames
        TempFrames = trackedPar(i).Frame;
        if length(TempFrames) >= MinFrames
            trackedEnsemble(iter).ResidenceFrames(iiter) = max(TempFrames)-min(TempFrames)+1;% Convert to column vector to be valid bootstrp input
            iiter = iiter + 1;
        end
    end
    
    %To make ResidenceFrames a valid input of bootstrp
    trackedEnsemble(iter).ResidenceFrames = trackedEnsemble(iter).ResidenceFrames';
    
    %Just to generate a good xbin frequency to evaluate cdf
    [~,temp_xi] = ecdf(trackedEnsemble(iter).ResidenceFrames, 'Function', 'survivor');
    trackedEnsemble(iter).ResidenceFramesGoneFreq = {temp_xi'};
    
    
    % Bootstrap the residence time of a sample
    % Each column in bs_ResidenceFramesIdx contains indices of the
    % values that were drawn from the original data sets to constitute
    % the corresponding bootstrap sample.
    [trackedEnsemble(iter).bs_ksdensity,trackedEnsemble(iter).bs_ResidenceFramesIdx] = bootstrp(nboot, @wrapped_ksdensity,...
        trackedEnsemble(iter).ResidenceFrames, trackedEnsemble(iter).ResidenceFramesGoneFreq);
    trackedEnsemble(iter).bs_ksdensity(1).ResidenceTime = trackedEnsemble(iter).bs_ksdensity(1).xi.*Exposure;
    trackedEnsemble(iter).bs_ksdensity(1).mean_f = mean(vertcat(trackedEnsemble(iter).bs_ksdensity(:).f),1);
    % Return the f matrix above the cutoff
    [trackedEnsemble(iter).bs_ksdensity(1).f_abovecutOff, maxColIdx_abovecutOff] = ...
        getMax_abovecutOff(vertcat(trackedEnsemble(iter).bs_ksdensity(:).f), yValue_cutOff, nboot);
    
    % Plot survival curve with shaded error bar
    s = shadedErrorBar(trackedEnsemble(iter).bs_ksdensity(1).ResidenceTime(1:maxColIdx_abovecutOff), ...
        trackedEnsemble(iter).bs_ksdensity(1).f_abovecutOff, {@mean, @std},...
        'lineProps',{'-', 'color', color_plate{iter}, 'linewidth', 3},'transparent',1);   %{lineColors{iter}, 'linewidth', 2},'transparent',false);
    set(s.edge, 'color', 'none');
    hold on;
        
    
end

%%  -- a section to add fit curves of H2B (iter == 1), or FOXA2 DE (iter == 3)---

if plot_fit_curve == 1
    
    iter = 3;
    StartTimeforFit =1; 
    legend_label = {'H2B-Halo DE', 'FOXA2-Halo APS', 'FOXA2-Halo DE',' 2-Exp','1-Exp'};  %'FOXA2-\DeltaDBD-Halo', 

    disp('bootstrapping sample trackedPar mat files....');
    disp(mat_files(iter));
    
    % Read the input track mat file generated from MTT
    file_name = char(mat_files(iter));
    trackedEnsemble(iter).Filenames = file_name(1:end-4);
    load([input_path filesep trackedEnsemble(iter).Filenames '.mat'],'trackedPar');
    iiter = 1;
    
    % Filter out the track length that are less than MinFrames
    for i=1:length(trackedPar)  
        %Account for missed frames/gap frames
        TempFrames = trackedPar(i).Frame;
        if length(TempFrames) >= MinFrames
            trackedEnsemble(iter).ResidenceFrames(iiter) = max(TempFrames)-min(TempFrames)+1;% Convert to column vector to be valid bootstrp input
            iiter = iiter + 1;
        end
    end
    
    %To make ResidenceFrames a valid input of bootstrp
    trackedEnsemble(iter).ResidenceFrames = trackedEnsemble(iter).ResidenceFrames';
    
    %Just to generate a good xbin frequency to evaluate cdf
    [~,temp_xi] = ecdf(trackedEnsemble(iter).ResidenceFrames, 'Function', 'survivor');
    trackedEnsemble(iter).ResidenceFramesGoneFreq = {temp_xi'};
    

    % Bootstrap the residence time of a sample
    % Each column in bs_ResidenceFramesIdx contains indices of the
    % values that were drawn from the original data sets to constitute
    % the corresponding bootstrap sample.
    [trackedEnsemble(iter).bs_ksdensity,trackedEnsemble(iter).bs_ResidenceFramesIdx] = bootstrp(nboot, @wrapped_ksdensity,...
        trackedEnsemble(iter).ResidenceFrames, trackedEnsemble(iter).ResidenceFramesGoneFreq);
    trackedEnsemble(iter).bs_ksdensity(1).ResidenceTime = trackedEnsemble(iter).bs_ksdensity(1).xi.*Exposure;
    trackedEnsemble(iter).bs_ksdensity(1).mean_f = mean(vertcat(trackedEnsemble(iter).bs_ksdensity(:).f),1);

        
    %StartTimeforFit ResidenceTime index
    StartTimeforFitIndx = ...
    find(trackedEnsemble(iter).bs_ksdensity(1).ResidenceTime == StartTimeforFit);
    %1-exp fitting, fitting return as a fit object
    trackedEnsemble(iter).OneExpFit.fitObject = fit(trackedEnsemble(iter).bs_ksdensity(1).ResidenceTime(StartTimeforFitIndx:end)', ...
        trackedEnsemble(iter).bs_ksdensity(1).mean_f(StartTimeforFitIndx:end)', 'exp1'); %input must be column vectors
    trackedEnsemble(iter).OneExpFit.yValue = feval(trackedEnsemble(iter).OneExpFit.fitObject, ...
        trackedEnsemble(iter).bs_ksdensity(1).ResidenceTime); %Evaluate the fit at given x (ResidenceTime)

    %2-exp fitting, fitting return as a fit object 
    trackedEnsemble(iter).TwoExpFit.fitObject = fit(trackedEnsemble(iter).bs_ksdensity(1).ResidenceTime(StartTimeforFitIndx:end)', ...
        trackedEnsemble(iter).bs_ksdensity(1).mean_f(StartTimeforFitIndx:end)', 'exp2'); %input must be column vectors
        %small coefficient is the slow rate
    trackedEnsemble(iter).TwoExpFit.slowRate = min(...
        -trackedEnsemble(iter).TwoExpFit.fitObject.b, -trackedEnsemble(iter).TwoExpFit.fitObject.d);
        %Evaluate the fit at given x (ResidenceTime)
    trackedEnsemble(iter).TwoExpFit.yValue = feval(trackedEnsemble(iter).TwoExpFit.fitObject, ...
        trackedEnsemble(iter).bs_ksdensity(1).ResidenceTime);        

    %Plot the 2-exp fitting curve, swap ResidenceTime to column vector
    plot(trackedEnsemble(iter).bs_ksdensity(1).ResidenceTime', ...
        trackedEnsemble(iter).TwoExpFit.yValue, 'Color', 'k', 'LineWidth', 1); %[0.8510 0.3255 0.3098
    %Plot the 1-exp fitting curve, swap ResidenceTime to column vector
    plot(trackedEnsemble(iter).bs_ksdensity(1).ResidenceTime', ...
        trackedEnsemble(iter).OneExpFit.yValue, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--');  %[0.8510 0.3255 0.3098
end

hold off;

%% linear: - Adjust plot appearance 
ax = gca; ax.LineWidth = 2; ax.Box = 'on';
ax.FontSize = 20; ax.XScale = 'linear'; 
ax.YScale = 'log'; %YScale set to log can see error more easily. (optional)
ax.YLim = [0.00005 0.3]; ax.XLim = [0 800]; 

title('SPT Survival Curve: 1-CDF', 'FontSize', 20);
xlabel('Trajectory length (sec)', 'FontSize', 20); ylabel ('1-CDF', 'FontSize', 20);
legend(legend_label, 'Location' ,'southwest', 'FontSize', 20);
legend boxoff;

hold off;

% Save figure as PDF
h=gcf;
set(h,'Units','Inches');
pos = get(h,'Position');

set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])  %pos(3) and pos(4) are width and height of this windows
print(h,[output_path, figure_name, '_linear.pdf'], '-dpdf', '-r0');


%% log: - Adjust plot appearance
ax = gca; ax.LineWidth = 2; ax.Box = 'on';
ax.FontSize = 20; ax.XScale = 'log'; 
ax.YScale = 'log'; %YScale set to log can see error more easily. (optional)
ax.YLim = [0.00005 1]; ax.XLim = [1 1000]; 


title('SPT Survival Curve: 1-CDF', 'FontSize', 20);
xlabel('Trajectory length (sec)', 'FontSize', 20); ylabel ('1-CDF', 'FontSize', 20);
legend(legend_label, 'Location' ,'southwest', 'FontSize', 20);
legend boxoff;

hold off;

% Save figure as PDF;

% figure 
h=gcf; h.OuterPosition = [ 1 1 4.8 8]; % unit is inch, so the value is small % [left bottom width height]
set(h,'Units','Inches');
pos = get(h,'Position');

set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])  %pos(3) and pos(4) are width and height of this windows
print(h,[output_path, figure_name, '_log_2.pdf'], '-dpdf', '-r0');

%%
disp('---- script completed ----');


%% Dependent local function of this script
function [f_abovecutOff, maxColIdx_abovecutOff] = getMax_abovecutOff (x, yValue_cutOff, nboot)

%getMax_abovecutOff Get the max matrix and column index of yValue matrix above cutoff value 
%   Return the max possible yValue matrix (a vertically concatenated f
%   value returned from ksdensity from bootstapped sample) and its max
%   column index above the defined cutoff value.
% Input:
%   x: the vertically concatenated f value returned from ksdensity from
%   bootstrapped samples.
%   yValue_cutOff: cutoff value that any 1-CDF below this number will be
%   discarded
%   nboot: bootnumber when using bootstrp()
% Output:
%   f_abovecutOff: the max possible yValue matrix (a vertically concatenated f
%   value returned from ksdensity from bootstapped sample).
%   maxColIdx_abovecutOff: max column index of the f_abovecutOff.

yValue_cutOff_Idx = x >= yValue_cutOff;
maxColIdx_abovecutOff = 1;
while maxColIdx_abovecutOff <= size(yValue_cutOff_Idx,2)
    if isequal(ones([nboot,1]),yValue_cutOff_Idx(:,maxColIdx_abovecutOff))
        maxColIdx_abovecutOff = maxColIdx_abovecutOff+1;
    else
        maxColIdx_abovecutOff = maxColIdx_abovecutOff-1;
        break;
    end
end%while
f_abovecutOff = x(:,1:maxColIdx_abovecutOff);
end%function



