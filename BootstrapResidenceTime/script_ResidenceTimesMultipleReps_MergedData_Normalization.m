% Script_ResidenceTimesMultipleReps_MergedData_Normalization
% Zuhui Wang 2020/04/07
% Wulan Version 2020 April 9th. 

% Purpose: (1) bootstraping data with or without photobleach correction; 
%          (2) to generate threshold fitting of residence time 
% Input: slowSPT merged trajectories files
% Output: (1) figure of residence time as a function of threshold fitting
%         (2) trackedEnsemble mat files with all data after 2-exp fitting,  (this is the input for  Script_ResidenceTimesMultipleReps_PlotOnly.m) 
% Note: This scrip will call function:  ResidenceTimesMultipleReps_Bootstrap

clc; clear all; close all;

%% shared parameters


Exposure = 0.5; %Second
nboot = 20; %number of bootstrap; default = 20
MinFrames = 2; %MinFrames: only consider trajectories with at least this many frames
doPhotoBleach = 0; % 1 for do photobleaching based on H2B slow SPT. 0 for not doing photobleacing correction
single_uncorrect = 0; % default is 0; change to 1 if want unique formatting of the figure when i want to plot single sample without correction. 
save_variables = 1; % whether you want to save the bootstrapped matfile, if you want to save it, change it from 0 to 1. 

H2B.bs_SPT_merged_SlowRate_mean = {NaN}; %must be 1x1 cell array to be valid bootstrp input;

% Make sure StartTimeForFit is the same as the ResidenceTimesMultipleReps_Bootstrap defined
StartTimeForFit = [Exposure*2,Exposure*3, Exposure*4, Exposure*5, Exposure*6, Exposure*7, Exposure*8, Exposure*9, Exposure*10];
StartTimeForPlot = Exposure*2;  % fixed start time for fitting. unit is Second.  StartFrameForFit = round(StartTimeForPlot / Exposure);
PlotIndex = find(StartTimeForFit==StartTimeForPlot); 

output_path = '/home/dell/Documents/METHOD/SPT/SPT_Analysis/denglabpku-repository/slowSMT-denglabpku/ExampleFigures/';

color_plate = {[0.1, 0.1,0.1], [0.8510 0.3255 0.3098], [0.3608 0.7216 0.3608], [0.3490 0.3490 0.4784], [0.3569 0.7529 0.8706]}; % color_plate{i}
% black [0, 0,0]
% grey [0.4118 0.4118 0.4118]
% light red  [0.8510 0.3255 0.3098]
% light blue [0.2588 0.5451 0.7922]
% light green [0.3608 0.7216 0.3608]
% purple [0.3922 0.1804 0.4863]

%% [photobleach corrected] Files to be analyzed and photo-bleaching corrected. 
% change legend label and figure name accordingly
if doPhotoBleach == 1
%---- H7 500ms ----%
input_path = '/mnt/f58069a5-1cf3-43b8-bb9b-ea74327327c9/WZH-DataCenter/Submission/2022_FOXA2Data/FOXA2_SMT_data/slowSMT_trajectories_pooled/';
H2B_mat = 'DE_H2B-Halo_2Hz_trackedPar.mat'; 

mat_files = {'APS_FOXA2-Halo_C14_2Hz_trackedPar.mat',  'DE_FOXA2-Halo_C14_2Hz_trackedPar.mat', 'DE_dCTD-Halo_C11_2Hz_trackedPar.mat', 'DE_dCdN-Halo_C10_2Hz_trackedPar.mat'};
legend_label = {'FOXA2-Halo APS', 'FOXA2-Halo DE', '\DeltaCTD-Halo DE', '\DeltaC\DeltaN-Halo DE'};
figure_name = 'ResTime_de_nboot20_minframe2_Corrected';
end

%% [photobleach NOT corrected]  500ms ----%
if doPhotoBleach == 0
%---- H7 500ms ----%
input_path = '/mnt/f58069a5-1cf3-43b8-bb9b-ea74327327c9/WZH-DataCenter/Submission/2022_FOXA2Data/FOXA2_SMT_data/slowSMT_trajectories_pooled/';
mat_files = {'DE_H2B-Halo_2Hz_trackedPar.mat', 'APS_FOXA2-Halo_C14_2Hz_trackedPar.mat',  'DE_FOXA2-Halo_C14_2Hz_trackedPar.mat', 'DE_dCTD-Halo_C11_2Hz_trackedPar.mat', 'DE_dCdN-Halo_C10_2Hz_trackedPar.mat'};
legend_label = {'H2B-Halo DE', 'FOXA2-Halo APS', 'FOXA2-Halo DE', '\DeltaCTD-Halo DE', '\DeltaC\DeltaN-Halo DE'};
figure_name = 'ResTime_de_nboot20_minframe2_Uncorrected';
end

%% Specify the matfile used for correction, usually use H2B
if doPhotoBleach == 1 
    disp('bootstrapping the H2B trackedPar mat files....');
    H2B_mat_file  = [input_path, H2B_mat];
    load(H2B_mat_file,'trackedPar');
    H2B = struct;

    iiter = 1;
    for i=1:length(trackedPar)  
        %Account for missed frames/gap frames
        TempFrames = trackedPar(i).Frame;
        if length(TempFrames) >= MinFrames
            H2B.ResidenceFrames(iiter) = max(TempFrames)-min(TempFrames)+1;% Convert to column vector to be valid bootstrp input
            iiter = iiter + 1;
        end
    end
    
    disp([num2str(iiter), ' out of ', num2str(length(trackedPar)), ' trajectories are used for bootstrapping and fitting.']);             
    bootstat = bootstrp(nboot,...
        @ResidenceTimesMultipleReps_Bootstrap, H2B.ResidenceFrames',Exposure); 
    for i = 1:length(bootstat)
        H2B.bs_ResidenceTime(:,i)=bootstat(i).merged_ResTime2Exp_uncorrected;
    end
    H2B.bs_ResidenceTime_mean = mean(H2B.bs_ResidenceTime, 2); 
    H2B.bs_SPT_merged_SlowRate = 1./H2B.bs_ResidenceTime;
    H2B.bs_SPT_merged_SlowRate_mean = {mean(H2B.bs_SPT_merged_SlowRate,2)}; %must be 1x1 cell array to be valid bootstrp input;
%             
end

%% Analyze samples with bootstraping and photo-bleach correction with H2B. (if doPhotoBleach = 0, correction will be omitted automatically)
trackedEnsemble = struct; %saving track ensemble data for downstream analysis
for iter = 1:length(mat_files)
    
    disp('bootstrapping sample trackedPar mat files....');
    disp(mat_files(iter));

    file_name = char(mat_files(iter));
    trackedEnsemble(iter).Filenames = file_name(1:end-4);
    load([input_path filesep trackedEnsemble(iter).Filenames '.mat'],'trackedPar');
    iiter = 1;
    for i=1:length(trackedPar)  
        %Account for missed frames/gap frames
        TempFrames = trackedPar(i).Frame;
        if length(TempFrames) >= MinFrames
            trackedEnsemble(iter).ResidenceFrames(iiter) = max(TempFrames)-min(TempFrames)+1;% Convert to column vector to be valid bootstrp input
            iiter = iiter + 1;
        end
    end
    trackedEnsemble(iter).note = strcat([num2str(iiter), ' out of ', num2str(length(trackedPar)), ' trajectories are used for bootstrapping and fitting.']);  % save this info for future reference WD
    disp(trackedEnsemble(iter).note);
 
    % export all bootstrap arguments instead of just the residence time. Sep. 2021, WD
    bootstat = bootstrp(nboot,...
        @ResidenceTimesMultipleReps_Bootstrap, trackedEnsemble(iter).ResidenceFrames',... % ResidenceFrames': Convert to column vector to be valid bootstrp input
        Exposure, doPhotoBleach,  H2B.bs_SPT_merged_SlowRate_mean);                                                                                                  % [merged_ResTime2Exp_uncorrected] = ResidenceTimesMultipleReps_Bootstrap (ResidenceFrames, Exposure)
    for i = 1:length(bootstat)

        %Swap the matrix to make one column corresponding to one bootstrapped dataset.
        trackedEnsemble(iter).bs_ResidenceTime(:,i)=bootstat(i).merged_ResTime2Exp_uncorrected; % if photobleach == 1, it is actually corrected value
        trackedEnsemble(iter).bs_SlowFraction(:,i)=bootstat(i).SPT_merged_SlowFraction;
        trackedEnsemble(iter).bs_FastRate(:,i)=bootstat(i).SPT_merged_FastRate;
        trackedEnsemble(iter).bs_SlowRate(:,i)=bootstat(i).SPT_merged_SlowRate;
      
        trackedEnsemble(iter).bs_nsResidenceTime(:,i)=1./bootstat(i).SPT_merged_FastRate;  % ns = non specific, it is the fast component 
        
    end
    
    % 2-Exp fit results
    trackedEnsemble(iter).bs_SlowFraction_mean = mean(trackedEnsemble(iter).bs_SlowFraction, 2);
        trackedEnsemble(iter).bs_SlowFraction_se = std(trackedEnsemble(iter).bs_SlowFraction,0, 2);
    trackedEnsemble(iter).bs_FastRate_mean = mean(trackedEnsemble(iter).bs_FastRate, 2);
        trackedEnsemble(iter).bs_FastRate_se = std(trackedEnsemble(iter).bs_FastRate,0,2);
    trackedEnsemble(iter).bs_SlowRate_mean = mean(trackedEnsemble(iter).bs_SlowRate, 2);
        trackedEnsemble(iter).bs_SlowRate_se = std(trackedEnsemble(iter).bs_SlowRate,0,2);
        
    % calculate specific and nonspecific Residence time
    trackedEnsemble(iter).bs_ResidenceTime_mean = mean(trackedEnsemble(iter).bs_ResidenceTime, 2); 
        trackedEnsemble(iter).bs_ResidenceTime_se = std(trackedEnsemble(iter).bs_ResidenceTime,0,2);
    trackedEnsemble(iter).bs_nsResidenceTime_mean = mean(trackedEnsemble(iter).bs_nsResidenceTime, 2); 
        trackedEnsemble(iter).bs_nsResidenceTime_se = std(trackedEnsemble(iter).bs_nsResidenceTime,0, 2); 
end

%% Plot specific residence time vs fitting threshold,  dot plot with errorbar

%color_plate = {'#d9534f', '#428bca', '#5cb85c', '#59597a', '#5bc0de'}; % color to plot each sample, use: color_plate{i} or char(color_plate(i))  matlab2019 can recognize hexadecimal color code, but earlier version can not. 
   
figure; %
for i = 1:length(trackedEnsemble)
    errorbar(StartTimeForFit, trackedEnsemble(i).bs_ResidenceTime_mean, trackedEnsemble(i).bs_ResidenceTime_se, '.', 'color',color_plate{i} , 'markersize', 50, 'linewidth', 2);
    hold on;
end
hold off
%plot = errorbar(plot_StartTimeForFit(PlotIndex:end,:), plot_bs_ResidenceTime_mean, plot_bs_ResidenceTime_se, '.', 'markersize', 50, 'linewidth', 2);

%% Adjust plot style

ax = gca; ax.FontSize = 20; ax.LineWidth = 2; 

xlabel('Fitting threshold (sec)', 'FontSize',20); 
[~, objh] = legend(legend_label,'Location' ,'southeast', 'FontSize',20);  % i don't really know what is the objh magic... 
legend boxoff

xlim([0.5 5.5]); %xlim([0.6 3.3]);
YLabeltemp = get(gca, 'YTickLabel');
ylim([0 min((str2double(YLabeltemp{end})+5), 500)]);
%ylim([0 60]);  % change the y axis limit according to your needs

% label the title and ylabel differently 
if doPhotoBleach == 0
    title(sprintf(['Uncorrected residence time']), 'FontSize',20);  % \n (inferred residence time as a function of fitting threshold)
    ylabel('Uncorrected residence time (sec)', 'FontSize',20);
end
if doPhotoBleach == 1 
    title(sprintf(['Photo-bleaching corrected residence time']), 'FontSize',20);
    ylabel('Corrected residence time (sec)', 'FontSize',20);
end

if single_uncorrect == 1
    title(sprintf(['Uncorrected residence time \n (number of trajectories: ', num2str(length(trackedEnsemble(1).ResidenceFrames)), ')']), 'FontSize',20);
    ax = gca; ax.FontSize = 20;
    [~, objh] = legend(legend_label,'Location' ,'southeast', 'FontSize',20);  % i don't really know what is the objh magic... 
    legend boxoff

end

%% Save figure as PDF; save variables

% save figure 
h=gcf;
set(h,'Units','Inches');
pos = get(h,'Position');

set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])  %pos(3) and pos(4) are width and height of this windows
print(h,[output_path, figure_name, '.pdf'], '-dpdf', '-r0');

% save variables  (default is not saving, 0) 
if save_variables == 1
    save([output_path, figure_name, '_info.mat'], 'StartTimeForFit', 'trackedEnsemble', 'figure_name', 'legend_label', 'nboot', 'MinFrames', 'doPhotoBleach');
end

disp('---- script completed ----');

