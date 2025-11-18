% Stats
% This can be used after running the code in loadData.m
% This is done to compare the filters on a smaller dataset to determine the
% parameters to be used on the full dataset.

fs = 30;
refHRs = new_T.RefHR;
Qual = new_T.Quality;

[HRs, bfilts, methods, s_sqi] = calc_stats(new_T,refHRs,fs);


%% Determine new Signal Quality
HR_diffs = abs(refHRs - HRs);

% make new signal quality
% 1 for acceptable and 0 for poor
new_Qual = double(HR_diffs < 5);

%% Save as table

% filtered dataset 1.
varNames = ["Ref HR", "Calculated HR", "Best Filter", "Best Method", "Motion", "Signal Quality"];
F1 = table(refHRs, HRs, bfilts, methods, new_T.Motion, new_Qual,'VariableNames',varNames);

%% HR stats
RMSE = rmse(HRs,refHRs)
MAE = mean(abs(refHRs-HRs))

% finding the number of the HR differences greater than 5 bpm.
% calculating the percentage of HR differences that were out of range.
num_outliers = length(find((HR_diffs > 5))) 
pct_outliers = (num_outliers/(length(refHRs)))*100

% number of new quality annotations that agree with the original annotations
num_agreement = length(find(Qual == new_Qual));

%% HR Calculation Method Stats
avg_method = mean(methods)
method_cts = groupcounts(methods)
mode_method = mode(methods)
% med_method = median(methods)

%% Filtering Method Stats
avg_filt = mean(bfilts)
filt_cts = groupcounts(bfilts)
mode_filt = mode(bfilts)
% med_filt = median(bfilts)

%% Signal Quality Stats

G1 = groupcounts(F1,["Motion","Signal Quality"],"IncludeEmptyGroups",true)
G2 = groupcounts(F1,["Best Method","Signal Quality"],"IncludeEmptyGroups",true)
G3 = groupcounts(F1,["Best Filter","Signal Quality"],"IncludeEmptyGroups",true)

%% Correlations

% R is correlation coefficient and tells linear relationship btw data -1 to 1.
% positive for positive correlation. neg for negative correlation
% if p value is less than 0.05, then R is significantly different from
% zero.
[R2, p2] = corr(HRs(:),refHRs(:))

% quantifies how well the model fits the data
R_squared = R2^2


%% helper functions

function [HRs, bfilts, methods, s_sqi] = calc_stats(T,refHRs,fs)
% ensure T is a table

    HRs = zeros(51,1);
    bfilts = zeros(51,1);
    methods = zeros(51,1);
    s_sqi = zeros(51,1);
    
    for i=1:length(refHRs)
        sig = T.Signal(i,:)';
        tm = T.Time(i,:)';
    
        % preprocessing/smoothing
        PPG=sig-movmean(sig, 30);
        [HRs(i), bfilts(i), methods(i), ~] = filterComparison(refHRs(i),PPG,tm,fs,0,0,22);
    
        % stats
        s_sqi(i) = skew_sqi(sig);
    
        disp([num2str((i/length(refHRs))*100),'% done.'])
    end
end

% calculates the sqi over a window and takes the average
function sqi = skew_sqi(PPG)
    w = 60; % number of samples in the window. must result in integer for array sqis
    sqis = zeros(length(PPG)/w,1);
    for i=w:w:length(PPG)
        sqi_sum = 0;
        % special case for first index
        if i==w
            pmean = mean(PPG(1:w));
            pstd = std(PPG(1:w));
            for j=1:w
                x_point = ((PPG(j) - pmean)/pstd)^3;
                sqi_sum = sqi_sum + x_point;
            end
            sqis(1) = (1/length(PPG(1:w)))*sqi_sum;
        else
            pmean = mean(PPG(i-w:i));
            pstd = std(PPG(i-w:i));
            for j=i-w:i
                x_point = ((PPG(j) - pmean)/pstd)^3;
                sqi_sum = sqi_sum + x_point;
            end
            sqis(i/w) = (1/length(PPG(i-w:i)))*sqi_sum;
        end
    end
    sqi = mean(sqis); %return the average skewness
end