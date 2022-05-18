% Script_ResidenceTimesMultipleReps_PlotOnly.m
% Wulan Deng, Sep. 2021

% Purpose: freely plot different samples on the same plot. 
% Input: file is generated from Script_ResidenceTimesMultipleReps_MergedData_Normalization; It
% contains bootstrapped RT with variable StartTimeForFit; it could be
% either photobleach corrected or uncorrected (different files)
% Output: (1) RT vs fitting threshold plot
%         (2) bar graph of specific RT at a specific fitting threshold, with errorbar 
% Note: this script does not call any function

%% if load matfile, run this, otherwise, omit these two lines. 
clc; close all; clear all
load('/Volumes/T7/WORKON-MANUSCRIPT/1-DATA_Analysis_Figure/slowSPT_Bootstrap_RT_curve_2021Sep/RT_uncorrected/ResTime_de_nboot20_minframe2_UNcorrected_1_info.mat');

%% Inputs (modify these variables)
ori_figure_name = figure_name;
figure_name = [figure_name '_de_mutant'];  % '_aps_de'    '_de_mutant'  '_aps_de_dc'
title_tag = 'Non-corrected ';  % Corrected or Uncorrected
plot_sample = [1,2,3,4, 5]; % this is to select the index of samples to be plotted

%% shared parameters that are not in the matfile 

single_uncorrect = 0;
output_path = '/Users/wdeng/Downloads/';
color_plate = {[0.1, 0.1,0.1], [0.8510 0.3255 0.3098], [0.3608 0.7216 0.3608], [0.5, 0, 0.5], [0.3569 0.7529 0.8706]}; % color_plate{i}
% black [0, 0,0]
% grey [0.4118 0.4118 0.4118]
% light red  [0.8510 0.3255 0.3098]
% light blue [0.2588 0.5451 0.7922]
% light green [0.3608 0.7216 0.3608]
% purple [0.3922 0.1804 0.4863]
% purple [0.5, 0, 0.5]
%color_plate = {'#d9534f', '#428bca', '#5cb85c', '#59597a', '#5bc0de'}; % color to plot each sample, use: color_plate{i} or char(color_plate(i))  matlab2019 can recognize hexadecimal color code, but earlier version can not. 
   
%% Plot specific residence time vs fitting threshold,  dot plot with errorbar
disp('---- plot specific RT vs Fitting threshold ----');

figure; 

for i = plot_sample
    disp(legend_label(i));
    errorbar(StartTimeForFit, trackedEnsemble(i).bs_ResidenceTime_mean, trackedEnsemble(i).bs_ResidenceTime_se, '.', 'color',color_plate{i} , 'markersize', 50, 'linewidth', 2);
    hold on;
end
hold off
%plot = errorbar(plot_StartTimeForFit(PlotIndex:end,:), plot_bs_ResidenceTime_mean, plot_bs_ResidenceTime_se, '.', 'markersize', 50, 'linewidth', 2);

%######## Adjust plot style

ax = gca; ax.FontSize = 20; ax.LineWidth = 2; 

xlabel('Fitting threshold (sec)', 'FontSize',20); 
[~, objh] = legend(legend_label(plot_sample),'Location' ,'southeast', 'FontSize',20);  % i don't really know what is the objh magic... 
legend boxoff

xlim([0.5 5.5]); %xlim([0.6 3.3]);
YLabeltemp = get(gca, 'YTickLabel');
ylim([0 min((str2double(YLabeltemp{end})+5), 500)]);
%ylim([0 60]);  % change the y axis limit according to your needs

% ###### label the title and ylabel differently 
if doPhotoBleach == 0
    title(sprintf(['Residence time without correction']), 'FontSize',20, 'FontWeight', 'Normal');  % \n (inferred residence time as a function of fitting threshold)
    ylabel('Residence time (sec)', 'FontSize',20);
end
if doPhotoBleach == 1 
    title(sprintf(['Photo-bleaching corrected residence time']), 'FontSize',20, 'FontWeight', 'Normal');
    ylabel('Corrected residence time (sec)', 'FontSize',20);
end

if single_uncorrect == 1
    title(sprintf(['Uncorrected residence time \n (number of trajectories: ', num2str(length(trackedEnsemble(1).ResidenceFrames)), ')']), 'FontSize',20, 'FontWeight', 'Normal');
    ax = gca; ax.FontSize = 20;
    [~, objh] = legend(legend_label,'Location' ,'southeast', 'FontSize',20);  % i don't really know what is the objh magic... 
    legend boxoff

end

%######## Save figure as PDF;
h=gcf;
set(h,'Units','Inches');
pos = get(h,'Position');

set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])  %pos(3) and pos(4) are width and height of this windows
print(h,[output_path, figure_name, '.pdf'], '-dpdf', '-r0');


%% Plot bar graph of specific RT at a specific fitting threshold, with errorbar

%% threshold = 1 sec
disp('---- plot bar graph of threshold = 1 ----');

plot_name = [figure_name, '_bar_t1'];

width = 60*length(label)+60;  % width of plot expand with more samples

% get relevant data of selected samples
j = 0;
for i = plot_sample
    j = j +1; 
    y1 (j) = trackedEnsemble(i).bs_ResidenceTime_mean(1);
    y1_se (j) = trackedEnsemble(i).bs_ResidenceTime_se(1);
    y2 (j) = trackedEnsemble(i).bs_nsResidenceTime_mean(1);
    y2_se (j) = trackedEnsemble(i).bs_nsResidenceTime_se(1);
    label(j) = legend_label(i);
    %b.CData(j,:) = color_plate{i}; % color bars 
end
x = categorical(label); % default is alphabetical
x = reordercats(x,cellstr(x)');  % change it to the original order

% plot the specific RT
figure; hold all;
% bar graph 
b = bar(x,y1,'FaceColor','flat');  % categorical(label) is the x label
ylabel('Residence time (sec)', 'FontSize', 20); title(sprintf(['Residence time without correction \n Fitting threshold = 1 sec']), 'FontSize',20, 'FontWeight', 'Normal');
% overlay with error bar
errorbar(x,y1, y1_se, '.', 'color', 'k', 'linewidth', 2); hold off % errorbar(x,y,err)
% adjust format
ax = gca; ax.FontSize = 20; ax.LineWidth = 2; ax.XAxis.TickLabelRotation = 90;  %ax.PlotBoxAspectRatio = [1 1.5 1];   %Relative length of each axis, specified as a three-element vector of the form [px py pz] 
h=gcf; h.OuterPosition = [ 100 100 width 600]; % [left bottom width height]  width = 266
% save figure 
h=gcf;set(h,'Units','Inches'); pos = get(h,'Position'); set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])  %pos(3) and pos(4) are width and height of this windows
print(h,[output_path, plot_name, '_RTs', '.pdf'], '-dpdf', '-r0');

% plot the nonspecific RT
figure; hold all;
% bar graph 
b = bar(x,y2,'FaceColor','flat');  % categorical(label) is the x label
ylabel('Residence time (sec)', 'FontSize', 20); title(sprintf(['Nonspecific residence time \n Fitting threshold = 1 sec']), 'FontSize',20, 'FontWeight', 'Normal');
% overlay with error bar
errorbar(x,y2, y2_se, '.', 'color', 'k', 'linewidth', 2); hold off % errorbar(x,y,err)
% adjust format
ax = gca; ax.FontSize = 20; ax.LineWidth = 2; ax.XAxis.TickLabelRotation = 90;  %ax.PlotBoxAspectRatio = [1 1.5 1];   %Relative length of each axis, specified as a three-element vector of the form [px py pz] 
h=gcf; h.OuterPosition = [ 100 100 width 600]; % [left bottom width height]
% save figure 
h=gcf;set(h,'Units','Inches'); pos = get(h,'Position'); set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])  %pos(3) and pos(4) are width and height of this windows
print(h,[output_path, plot_name, '_RTns', '.pdf'], '-dpdf', '-r0');

%% threshold = 2 sec
disp('---- plot bar graph of threshold = 2 ----');

plot_name = [figure_name, '_bar_t2'];
% get relevant data of selected samples
j = 0;
for i = plot_sample
    j = j +1; 
    y3 (j) = trackedEnsemble(i).bs_ResidenceTime_mean(3);
    y3_se (j) = trackedEnsemble(i).bs_ResidenceTime_se(3);
    y4 (j) = trackedEnsemble(i).bs_nsResidenceTime_mean(3);
    y4_se (j) = trackedEnsemble(i).bs_nsResidenceTime_se(3);
    label(j) = legend_label(i);
   % b.CData(j,:) = color_plate{i}; % color bars 
end
x = categorical(label); % default is alphabetical
x = reordercats(x,cellstr(x)');  % change it to the original order

% plot the specific RT
figure; hold all;
% bar graph 
b = bar(x,y3,'FaceColor','flat');  % categorical(label) is the x label
ylabel('Residence time (sec)', 'FontSize', 20); title(sprintf(['Residence time without correction \n Fitting threshold = 2 sec']), 'FontSize',20, 'FontWeight', 'Normal');
% overlay with error bar
errorbar(x,y3, y3_se, '.', 'color', 'k', 'linewidth', 2); hold off % errorbar(x,y,err)
% adjust format 
ax = gca; ax.FontSize = 20; ax.LineWidth = 2; ax.XAxis.TickLabelRotation = 90;  %ax.PlotBoxAspectRatio = [1 1.5 1];   %Relative length of each axis, specified as a three-element vector of the form [px py pz] 
h=gcf; h.OuterPosition = [ 100 100 width 600]; % [left bottom width height][ 100 100 366 600];
% save figure 
h=gcf;set(h,'Units','Inches'); pos = get(h,'Position'); set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])  %pos(3) and pos(4) are width and height of this windows
print(h,[output_path, plot_name, '_RTs', '.pdf'], '-dpdf', '-r0');

% plot the nonspecific RT
figure; hold all;
% bar graph 
b = bar(x,y4,'FaceColor','flat');  % categorical(label) is the x label
ylabel('Residence time (sec)', 'FontSize', 20); title(sprintf([ 'Nonspecific residence time \n Fitting threshold = 2 sec']), 'FontSize',20, 'FontWeight', 'Normal');
% overlay with error bar
errorbar(x,y4, y4_se, '.', 'color', 'k', 'linewidth', 2); hold off % errorbar(x,y,err)
% adjust format
ax = gca; ax.FontSize = 20; ax.LineWidth = 2; ax.XAxis.TickLabelRotation = 90;  %ax.PlotBoxAspectRatio = [1 1.5 1];   %Relative length of each axis, specified as a three-element vector of the form [px py pz] 
h=gcf; h.OuterPosition = [ 100 100 width 600]; % [left bottom width height]
% save figure 
h=gcf;set(h,'Units','Inches'); pos = get(h,'Position'); set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])  %pos(3) and pos(4) are width and height of this windows
print(h,[output_path, plot_name, '_RTns', '.pdf'], '-dpdf', '-r0');

%% RTs with correction (calculated from slow rate, not bootstrap)

% get data 
N = 1; % select fitting time, default is t = 1sec
for i = 1:5  % this is num of cell samples 
 
%     RTs (i) = trackedEnsemble(i).bs_ResidenceTime_mean(N);
%     RTs_se (i) = trackedEnsemble(i).bs_ResidenceTime_se(N);
%     RTns (i) = trackedEnsemble(i).bs_nsResidenceTime_mean(N);
%     RTns_se (i) = trackedEnsemble(i).bs_nsResidenceTime_se(N);

    slow(i)=  trackedEnsemble(i).bs_SlowRate_mean(N);
    slow_se(i) = trackedEnsemble(i).bs_SlowRate_se(N);
    fast(i) = trackedEnsemble(i).bs_FastRate_mean(N);
    fast_se(i) = trackedEnsemble(i).bs_FastRate_se(N);
%     
%     fraction(i) = trackedEnsemble(i).bs_SlowFraction_mean(N);
%     fraction_se(i) = trackedEnsemble(i).bs_SlowFraction_se(N);
    
    % use slow rate to calculate normalized RTs
    if i ~=1
        cal_norm_RTs (i) = 1/(slow(i)-slow(1));   % calculated normalized RTs from slow rate, this is different from RTs, but they are very close. 
        cal_norm_RTs_se (i) = (1/(slow(i)-slow(1)))-(1/(slow(i)+slow_se(i)-slow(1)));
        
        cal_norm_RTns (i) = 1/(fast(i)-slow(1));   % calculated normalized RTs from slow rate, this is different from RTs, but they are very close. 
        cal_norm_RTns_se (i) = (1/(fast(i)-slow(1)))-(1/(fast(i)+fast_se(i)-slow(1)));
        
    end
      
end

disp('---- plot bar graph calculated RT with correction ----');

plot_name = [figure_name, 'Cal_Norm_bar_t1'];
% get relevant data of selected samples
j = 0;
for i = plot_sample(2:end)
    j = j +1; 
    y5 (j) = cal_norm_RTs (i);
    y5_se (j) = cal_norm_RTs_se (i);
    y6 (j) = cal_norm_RTns (i);
    y6_se (j) = cal_norm_RTns_se (i);
    label(j) = legend_label(i);
   % b.CData(j,:) = color_plate{i}; % color bars 
end
x = categorical(label); % default is alphabetical
x = reordercats(x,cellstr(x)');  % change it to the original order

% plot the specific RT
figure; hold all;
% bar graph 
b = bar(x,y5,'FaceColor','flat');  % categorical(label) is the x label
ylabel('Residence time (sec)', 'FontSize', 20); title(sprintf(['RTs with correction \n Fitting threshold = 1 sec']), 'FontSize',20, 'FontWeight', 'Normal');
% overlay with error bar
errorbar(x,y5, y5_se, '.', 'color', 'k', 'linewidth', 2); hold off % errorbar(x,y,err)
% adjust format 
ax = gca; ax.FontSize = 20; ax.LineWidth = 2; ax.XAxis.TickLabelRotation = 90;  %ax.PlotBoxAspectRatio = [1 1.5 1];   %Relative length of each axis, specified as a three-element vector of the form [px py pz] 
h=gcf; h.OuterPosition = [ 100 100 width 600]; % [left bottom width height][ 100 100 366 600];
% save figure 
h=gcf;set(h,'Units','Inches'); pos = get(h,'Position'); set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])  %pos(3) and pos(4) are width and height of this windows
print(h,[output_path, plot_name, '_RTs', '.pdf'], '-dpdf', '-r0');

% plot the nonspecific RT
figure; hold all;
% bar graph 
b = bar(x,y6,'FaceColor','flat');  % categorical(label) is the x label
ylabel('Residence time (sec)', 'FontSize', 20); title(sprintf([ 'RTns with correction \n Fitting threshold = 1 sec']), 'FontSize',20, 'FontWeight', 'Normal');
% overlay with error bar
errorbar(x,y6, y6_se, '.', 'color', 'k', 'linewidth', 2); hold off % errorbar(x,y,err)
% adjust format
ax = gca; ax.FontSize = 20; ax.LineWidth = 2; ax.XAxis.TickLabelRotation = 90;  %ax.PlotBoxAspectRatio = [1 1.5 1];   %Relative length of each axis, specified as a three-element vector of the form [px py pz] 
h=gcf; h.OuterPosition = [ 100 100 width 600]; % [left bottom width height]
% save figure 
h=gcf;set(h,'Units','Inches'); pos = get(h,'Position'); set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])  %pos(3) and pos(4) are width and height of this windows
print(h,[output_path, plot_name, '_RTns', '.pdf'], '-dpdf', '-r0');

% plot

disp('---- script completed ----');
%% plot part 
% plot the specific RT
x = categorical(label(1:3)); % default is alphabetical
x = reordercats(x,cellstr(x)');  % change it to the original order
y5 = y5(1:3); y5_se = y5_se(1:3); 
figure; hold all;
% bar graph 
b = bar(x,y5,'FaceColor','flat');  % categorical(label) is the x label
ylabel('Residence time (sec)', 'FontSize', 20); title(sprintf(['RTs with correction \n Fitting threshold = 1 sec']), 'FontSize',20, 'FontWeight', 'Normal');
% overlay with error bar
errorbar(x,y5, y5_se, '.', 'color', 'k', 'linewidth', 2); hold off % errorbar(x,y,err)
% adjust format 
ax = gca; ax.FontSize = 20; ax.LineWidth = 2; ax.XAxis.TickLabelRotation = 90;  %ax.PlotBoxAspectRatio = [1 1.5 1];   %Relative length of each axis, specified as a three-element vector of the form [px py pz] 
h=gcf; h.OuterPosition = [ 100 100 300 600]; % [left bottom width height][ 100 100 366 600];
% save figure 
h=gcf;set(h,'Units','Inches'); pos = get(h,'Position'); set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])  %pos(3) and pos(4) are width and height of this windows
print(h,[output_path, plot_name, '_RTs_2', '.pdf'], '-dpdf', '-r0');
disp('---- script completed ----');
%% save variables for statistics significance calculation
StatT1 = struct;
StatT1.label  = legend_label;
StatT1.threshold = '1sec';
StatT1.RTs_mean = y1; 
StatT1.RTs_std = y1_se;
StatT1.RTns_mean = y2; 
StatT1.RTns_std = y2_se;
StatT1.RTs_CalNorm_mean = y5; 
StatT1.RTs_CalNorm_std = y5_se;
StatT1.RTns_CalNorm_mean = y6; 
StatT1.RTns_CalNorm_std = y6_se;

save([output_path, ori_figure_name, '_stat_sig.mat'], 'StatT1');