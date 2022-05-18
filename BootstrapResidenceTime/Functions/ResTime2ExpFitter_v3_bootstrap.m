function [ SlowFraction, FastRate, SlowRate, PrefixA, FractionBound, FittingTrajNum] = ResTime2ExpFitter_v3_bootstrap(ResidenceFrames, StartFrameForFit, Exposure)
%RESTIME2EXPFITTER this function fits 2-exponentials to SPT data 
%StartTimeForFit: only fit after this time threshold
%MinFrames: only consider trajectories with at least this many frames
%Zuhui Wang: move out ResidenceFrames variable.


FittingTrajNum = length(ResidenceFrames);  % added by WD to show trajectory # used for fitting
%Empirical CDF

%Bin into a full histogram:
HistVec = 0:1:max(ResidenceFrames)+1;
HistVecTime = Exposure.*[HistVec(2:end) HistVec(end)+1];
ResidenceProb = histc(ResidenceFrames, HistVec)./length(ResidenceFrames); %must use histc not histcounts!  by ZW
ResidenceCDF = zeros(1,length(ResidenceProb));
for i=2:length(ResidenceProb)
    ResidenceCDF(1,i) = sum(ResidenceProb(1:i));
end
%ResidenceProb = ResidenceProb(2:end); %since single count frames is the second element
FractionBound = 1-ResidenceCDF;


%%%%%%%%%%%%%%%%%%%%%%%% CURVE FITTING %%%%%%%%%%%%%%%%%%%%%%%%

%Fit two exponentials to all of the data
yProb = FractionBound(StartFrameForFit:end);
xTime = HistVecTime(StartFrameForFit:end);

f = fittype('A*(F*exp(-a*x) + (1-F)*exp(-b*x))');
% 
%[TwoExp_fit, TwoExp_param] = fit(xTime', yProb', f, 'Lower', [0.8 0 0.3 0], 'Upper', [25 1 5 1], 'StartPoint', [4 0.9 2 0.02]); 
% this is default value by Anders. I am changing it because some of my data
% hit the lower boundary 0.3
[TwoExp_fit, TwoExp_param] = fit(xTime', yProb', f, 'Lower', [0.8 0 0.2 0], 'Upper', [25 1 5 1], 'StartPoint', [4 0.9 2 0.02]); 
% [A, F, a, b], a is the decay rate of nonspecific fraction, b is the decay rate of SPECIFIC fraction.
%'Upper' ? Upper bounds on coefficients to be fitted; 'StartPoint' ? Initial values for the coefficients
TwoExp_CI = confint(TwoExp_fit);

SlowFraction = 1-TwoExp_fit.F;
FastRate = TwoExp_fit.a;
SlowRate = TwoExp_fit.b;
PrefixA = TwoExp_fit.A;

end


