function [HR,best_filter,method, best_sig] = filterComparison(varargin)
% compare the different filters and return the HR, best filter, and best HR
% calculation method
%
% [HR, best_filter, method] = filterComparison(refHR, sig, tm, fs, plot_best, plot_all,wind)
% 
% Outputs:
% HR
%   integer, average heart rate in bpm
% best_filter
    % string, indicates the best filter method used
% method
%   integer or string, indicating the method that determined the most accurate HR
% 
% Required parameters:
% refHR
%   integer specifying the reference HR taken from the ECG
% sig
%   vector containing the signal
%
% Optional parameters:
% tm
%   vector containing the time in SECONDS
% fs
%   integer containing the signaling frequency in Hz, default 30
% plot_best
%   integer specifying to plot results of best filter
%   (1 for true, anything else yields false). default 0
% plot_all
%   integer specifying to plot ALL the filtering results
%   (1 for true, anything else yields false). default 0
% wind
%   integer specifying window length for movmax

    %set default values 
    inputs ={'refHR','sig','tm','fs','plot_best','plot_all','wind'};
    fs = 30;
    tm = linspace(0,300/fs,300);
    plot_best = 0;
    plot_all = 0;
    wind = 30;

    %assigning variables based off input
    for n=1:nargin
        if(~isempty(varargin{n}))
            eval([inputs{n} '=varargin{n};']);
        end
    end
    

    % no filter
    tic;
    PPG_flat = sys_norm(sig,wind);
    [HR1, method1] = det_HR(refHR,PPG_flat);
    toc1 = toc;

    % butterworth filter
    tic;
    [b, a] = butter(2, [0.5 10] / (fs / 2), 'bandpass');
    PPG_butter = filtfilt(b, a, sig);
    PPG_flat_butter = sys_norm(PPG_butter,wind);
    [HR2,method2] = det_HR(refHR,PPG_flat_butter);
    toc2 = toc;

    % IIR filtering
    tic;
    PPG_filt=highpass(sig,0.8,fs,'ImpulseResponse','iir','Steepness',0.5); % baseline wander
    PPG_filt=lowpass(PPG_filt,3,fs,'ImpulseResponse','iir','Steepness',0.8); % powerline noise
    PPG_flat_IIR = sys_norm(PPG_filt,wind);
    [HR3, method3] = det_HR(refHR,PPG_flat_IIR);
    toc3 = toc;

    % 4th order Chebyshev type II filter
    tic;
    [b, a] = cheby2(4, 20, [0.5 10] / (fs / 2));
    PPG_cheby = filtfilt(b, a, sig);
    PPG_flat_cheby = sys_norm(PPG_cheby,wind);
    [HR4, method4] = det_HR(refHR,PPG_flat_cheby);
    toc4 = toc;

    % EMD filtering
    tic;
    IMF = emd(sig);
    PPG_IMF = sum(IMF(2:3,:),1);
    PPG_flat_IMF = sys_norm(PPG_IMF,wind);
    [HR5, method5] = det_HR(refHR,PPG_flat_IMF);
    toc5 = toc;

    % determining best filter method
    rates = [HR1 HR2 HR3 HR4 HR5];
    times = [toc1 toc2 toc3 toc4 toc5];
    [HR, filter_num] = compareHR(refHR, rates, times);
    
    switch filter_num
        case 1
            % best_filter = "No filter";
            best_filter = 0;
            method = method1;
            best_sig = PPG_flat;
        case 2
            % best_filter = "Butterworth filter";
            best_filter = 1;
            method = method2;
            best_sig = PPG_flat_butter;
        case 3
            % best_filter = "IIR filters";
            best_filter = 2;
            method = method3;
            best_sig = PPG_flat_IIR;
        case 4
            % best_filter = "Chebyshev filter";
            best_filter = 3;
            method = method4;
            best_sig = PPG_flat_cheby;
        case 5
            % best_filter = "EMD filter";
            best_filter = 4;
            method = method5;
            best_sig = PPG_flat_IMF;
    end

    if plot_best == 1
        figure;
        plot(tm,best_sig); xlabel('Time (s)')
        title(['PPG Filtered.', 'HR: ', num2str(HR), ' bpm. Ref HR: ', num2str(refHR), ' bpm.'])
        subtitle(['Best Filter: ', num2str(best_filter), ...
            ' . Best method to find HR: ', num2str(method), '.'])
    end

    if plot_all == 1
        figure;
        T = tiledlayout(3,2,'TileSpacing','compact','Padding','compact');

        nexttile
        plot(tm,sig); title('Original Signal');
        subtitle(['RefHR: ', num2str(refHR), ' bpm'])

        nexttile
        plot(tm,PPG_flat); title('No filter');
        subtitle(['HR: ', num2str(round(HR1,2)), ' bpm', '. Method: ', num2str(method1)])

        nexttile
        plot(tm,PPG_flat_butter); title('Butterworth filter')
        subtitle(['HR: ', num2str(round(HR2,2)), ' bpm', '. Method: ', num2str(method2)])

        nexttile
        plot(tm,PPG_flat_IIR); title('IIR filters')
        subtitle(['HR: ', num2str(round(HR3,2)), ' bpm', '. Method: ', num2str(method3)])

        nexttile
        plot(tm,PPG_flat_cheby); title('4th Order Chebyshev type II filter')
        subtitle(['HR: ', num2str(round(HR4,2)), ' bpm', '. Method: ', num2str(method4)])
        xlabel('Time (s)')

        nexttile
        plot(tm,PPG_flat_IMF); title('EMD filtering')
        subtitle(['HR: ', num2str(round(HR5,2)), ' bpm', '. Method: ', num2str(method5)])
        xlabel('Time (s)')

        title(T,'PPG Filter Comparison')
        st = subtitle(T,['Best filter: ', num2str(best_filter), '. Best method: ', num2str(method)]);
        st.FontSize = 10;
        
    end

end

function PPG_flat = sys_norm(PPG_filt,wind)
    PPG_trend = movmax(abs(PPG_filt),wind); % track systolic peaks
    PPG_flat = PPG_filt./PPG_trend; % normalize and flatten
end