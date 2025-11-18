
% HR1 is finding HR with findpeaks. Method 1 in det_HR.
% HR2 is finding HR with zero padding fft. Method 5 in det_HR.
% Before running this, load in necessary tables.
% Used on the full dataset 3888, and the original dataset 48

%% For the Full DataSet 
L = 3888;

no_filter = loadTable(All,0,22,L);
no_filter_stats = table_Stats(no_filter,All);

butter_filter = loadTable(All,1,22,L);
butter_stats = table_Stats(butter_filter,All);

IIR_filter = loadTable(All,2,22,L);
IIR_stats = table_Stats(IIR_filter,All);

cheby_filter = loadTable(All,3,22,L);
cheby_stats = table_Stats(cheby_filter,All);

%% For the 48 Original Dataset
L = 48;

no_filter_48 = loadTable(All,0,22,L);
no_filter_48_stats = table_Stats(no_filter_48,All);

butter_filter_48 = loadTable(All,1,22,L);
butter_48_stats = table_Stats(butter_filter_48,All);

IIR_filter_48 = loadTable(All,2,22,L);
IIR_stats_48 = table_Stats(IIR_filter_48,All);

cheby_filter_48 = loadTable(All,3,22,L);
cheby_stats_48 = table_Stats(cheby_filter_48,All);


%% Motion Stats
groupcounts(T48,["Motion","Quality"],"IncludeEmptyGroups",true)
%% Location Stats
groupcounts(T48,["Location","Quality"],"IncludeEmptyGroups",true)

%% Plotting
newcolors = ["#721f81" "#f1605d"];

[no_err, no_mae] = hrMAE(no_filter_48);
[but_err, but_mae] = hrMAE(butter_filter_48);
[iir_err, iir_mae] = hrMAE(IIR_filter_48);
[cheb_err, cheb_mae] = hrMAE(cheby_filter_48);

X = categorical({'No filter','Butterworth','IIR','Chebyshev'});
X = reordercats(X,{'No filter','Butterworth','IIR','Chebyshev'});

Y = [no_mae; but_mae; iir_mae; cheb_mae];

% error bar
E = [no_err; but_err; iir_err; cheb_err];
ngroups = size(Y, 1);
nbars = size(Y, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));

figure;
colororder(newcolors)
bar(X,Y); hold on;
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, Y(:,i), E(:,i), 'o','LineStyle','none','Color','k');
end
ylim([0 20])
xlabel('Filter Type')
ylabel(["MAE (bpm)","|HR_{Ref} - HR_{Calculated}|"])
title('MAE for filters (48)')
legend('HR1','HR2')
hold off


%%
[no_err, no_mae] = maeMotion(no_filter_48);
[but_err, but_mae] = maeMotion(butter_filter_48);
[iir_err, iir_mae] = maeMotion(IIR_filter_48);
[cheb_err, cheb_mae] = maeMotion(cheby_filter_48);

X = categorical({'No filter','Butterworth','IIR','Chebyshev'});
X = reordercats(X,{'No filter','Butterworth','IIR','Chebyshev'});

Y = [no_mae(1,:); but_mae(1,:); iir_mae(1,:); cheb_mae(1,:)];

% error bar
E = [no_err(1,:); but_err(1,:); iir_err(1,:); cheb_err(1,:)];
ngroups = size(Y, 1);
nbars = size(Y, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));

figure;
tiledlayout

nexttile
colororder(newcolors)
bar(X, Y); hold on;
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, Y(:,i), E(:,i), 'o','LineStyle','none','Color','k');
end
xlabel('Filter Type')
ylabel(["MAE (bpm)","|HR_{Ref} - HR_{Calculated}|"])
title('MAE During Rest and Motion for HR1')
legend('Rest','Motion','Location','northwest')
ylim([0 35])
hold off

nexttile
Y = [no_mae(2,:); but_mae(2,:); iir_mae(2,:); cheb_mae(2,:)];
% error bar
E = [no_err(2,:); but_err(2,:); iir_err(2,:); cheb_err(2,:)];
ngroups = size(Y, 1);
nbars = size(Y, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));

colororder(newcolors)
bar(X, Y); hold on;
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, Y(:,i), E(:,i), 'o','LineStyle','none','Color','k');
end
xlabel('Filter Type')
ylabel(["MAE (bpm)","|HR_{Ref} - HR_{Calculated}|"])
title('MAE During Rest and Motion for HR2')
legend('Rest','Motion','Location','northwest')
ylim([0 35])
hold off

%%
function [err_val, mae_val] = maeMotion(T)
%calculate MAE and SEM for rest and activity. for plotting
% before using on table. turn all values that are not 0, to 1 for grouping
 
    nonzero_HR1 = T.HR1 ~= 0;
    nonzero_HR2 = T.HR2 ~= 0;

    [group, Motion] = findgroups(T.Motion);
    func = @(x,y) mean(abs(x-y));
    func2 = @(x,y) (std(abs(x-y)))/(sqrt(length(abs(x-y))));

    % calculate mae and sem for motion groups
    mae_1 = splitapply(func,T.("Ref HR")(nonzero_HR1),T.HR1(nonzero_HR1),group); % [rest , motion]
    mae_2 = splitapply(func,T.("Ref HR")(nonzero_HR2),T.HR2(nonzero_HR2),group);

    sem_1 = splitapply(func2,T.("Ref HR")(nonzero_HR1),T.HR1(nonzero_HR1),group);
    sem_2 = splitapply(func2,T.("Ref HR")(nonzero_HR2),T.HR2(nonzero_HR2),group);
    
    mae_val = [mae_1.'; mae_2.'];
    err_val = [sem_1.'; sem_2.'];
    % result = table(Motion,mae_1,mae_2,sem_1,sem_2); % turn into table for easier indexing

end

function [err_val,mae_val] = hrMAE(T)
%for plotting
    nonzero_HR1 = T.HR1 ~= 0;
    nonzero_HR2 = T.HR2 ~= 0;
    % calcualte error bar
    ab_er1 = abs(T.("Ref HR")(nonzero_HR1)-T.HR1(nonzero_HR1));
    ab_er2 = abs(T.("Ref HR")(nonzero_HR2)-T.HR2(nonzero_HR2));

    std1 = std(ab_er1);
    std2 = std(ab_er2);

    sem1 = std1/sqrt(length(ab_er1)); %standard error of MAE
    sem2 = std2/sqrt(length(ab_er2));
    err_val = [sem1 sem2];

    % calculate mae
    mae1 = mean(ab_er1);
    mae2 = mean(ab_er2);
    mae_val = [mae1 mae2];
end

function tableStats = table_Stats(T,All)
    n = size(T,1);

    % Print the quality stats
    HR1_qual_stat = groupcounts(T, "HR1 Sig Quality","IncludeEmptyGroups",true);
    HR2_qual_stat = groupcounts(T, "HR2 Sig Quality","IncludeEmptyGroups",true);
    disp(HR1_qual_stat)
    disp(HR2_qual_stat)
    
    % RMSE
    RMSE1 = rmse(T.HR1,T.("Ref HR"));
    RMSE2 = rmse(T.HR2,T.("Ref HR"));
    
    % MAE
    nonzero_HR1 = T.HR1 ~= 0;
    nonzero_HR2 = T.HR2 ~= 0;
    MAE1 = mean(abs(T.("Ref HR")(nonzero_HR1)-T.HR1(nonzero_HR1)));
    MAE2 = mean(abs(T.("Ref HR")(nonzero_HR2)-T.HR2(nonzero_HR2)));
    
    % Pearson correlation stats
    [R1, p_R1] = corr(T.HR1(:),T.("Ref HR")(:));
    [R2, p_R2] = corr(T.HR2(:),T.("Ref HR")(:));
    
    R_squared1 = R1^2;
    R_squared2 = R2^2;
    
    % quality agreement in original file and newly calculated quality
    num_qual_agreement1 = length(find(All.Quality(1:n) == T.("HR1 Sig Quality")));
    num_qual_agreement2 = length(find(All.Quality(1:n) == T.("HR2 Sig Quality")));
    
    varNames = ["Ttest Decision","Ttest P-value","RMSE","MAE", "Pearson R",...
        "Pearson R^2","Pearson P-value", "HR Quality Agreement"];

    tableStats = table([h1;h2],[p1;p2], [RMSE1;RMSE2], [MAE1;MAE2],[R1;R2], ...
        [R_squared1; R_squared2],[p_R1; p_R2],[num_qual_agreement1; num_qual_agreement2], ...
        'VariableNames',varNames,'RowNames',["HR1"; "HR2"]);
end

function name = loadTable(T,filt_type,wind,n)
% T is the table containing the information from the CSV combined with the
% signal files. Filt_type is an integer (0-3) to specify which filter to
% use. Wind is the window for movmax for flattening at the end. 
% n will be the number of files to use. starts from 1 to 3888
    HR1 = zeros(n,1);
    HR2 = zeros(n,1);
    refHRs = T.RefHR(1:n);

    for i=1:n
        sig = T.Signal(i,:)';
        % tm = T.Time(i,:)';
    
        % preprocessing/smoothing
        PPG=sig-movmean(sig, 30);
        %filter function
        [HR1(i),HR2(i)] = filterAll(PPG,filt_type,wind);
    
        disp([num2str((i/length(refHRs))*100),'% done.'])
    end
    
    % Determine new Signal Quality
    HR_diffs1 = abs(refHRs - HR1);
    HR_diffs2 = abs(refHRs - HR2);
    % 1 for acceptable and 0 for poor
    new_Qual1 = double(HR_diffs1 < 5);
    new_Qual2 = double(HR_diffs2 < 5);

    % new table
    varNames = ["Ref HR", "HR1", "HR2", "Location", "Motion","HR1 Sig Quality", "HR2 Sig Quality"];
    name = table(refHRs, HR1, HR2, T.Location(1:n), T.Motion(1:n), new_Qual1, new_Qual2,'VariableNames',varNames);

end

function [HR1,HR2] = filterAll(sig, filt_type, wind)
    fs = 30;

    switch filt_type
        case 0
            % no filter
            PPG_flat = sys_norm(sig,wind);
            [HR1, HR2] = getHR(PPG_flat,fs);
        case 1
            % butterworth filter
            [b, a] = butter(2, [0.5 10] / (fs / 2), 'bandpass');
            PPG_butter = filtfilt(b, a, sig);
            PPG_flat_butter = sys_norm(PPG_butter,wind);
            [HR1,HR2] = getHR(PPG_flat_butter,fs);
        case 2
            % IIR filter
            PPG_filt=highpass(sig,0.8,fs,'ImpulseResponse','iir','Steepness',0.5); % baseline wander
            PPG_filt=lowpass(PPG_filt,3,fs,'ImpulseResponse','iir','Steepness',0.8); % powerline noise
            PPG_flat_IIR = sys_norm(PPG_filt,wind);
            [HR1, HR2] = getHR(PPG_flat_IIR,fs);
        case 3
            % Chebyshev filter
            [b, a] = cheby2(4, 20, [0.5 10] / (fs / 2));
            PPG_cheby = filtfilt(b, a, sig);
            PPG_flat_cheby = sys_norm(PPG_cheby,wind);
            [HR1, HR2] = getHR(PPG_flat_cheby,fs);
    end

end

function PPG_flat = sys_norm(PPG_filt,wind)
    PPG_trend = movmax(abs(PPG_filt),wind); % track systolic peaks
    PPG_flat = PPG_filt./PPG_trend; % normalize and flatten
end

function [HR1,HR2] = getHR(sig,fs)

    [~, locs] = findpeaks(sig, 'MinPeakProminence',0.05, 'MinPeakDistance',8,'MinPeakHeight',0.45);
    interval = diff(locs/fs); %convert to seconds
    HR_1 = (60/mean(interval));
    if (isnan(HR_1)) | isempty(HR_1)
        HR1 = 0;
    else 
        HR1 = HR_1;
    end

    N_fft = length(sig);
    n = 2^nextpow2(N_fft);
    f2 = fs*(0:(n-1))/n;

    zero_fft = abs((fft(sig, n))/N_fft);

    % limit frequency range to 0.5-4 Hz
    mask = (0.5 < f2) & (f2 < 4);
    limit = zero_fft(mask);
    f_limit = f2(mask);
    
    % find peak within that range
    [peaks, loc] = findpeaks(limit);
    [~, idx] = max(peaks); %max value and its index
    peak_freq = f_limit(loc(idx));
    HR_2 = (peak_freq*60);
    if (isempty(HR_2)) | isnan(HR_2)
        HR2 = 0;
    else 
        HR2 = HR_2;
    end

end
