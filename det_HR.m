function [HR,method] = det_HR(varargin)
% Determine best method of finding HR and return it
% can be modified later to include statistics
%
% [HR, method] = det_HR(refHR, sig, tm, fs, plots, minProm, minDist, minHeight)
%
% Outputs:
% HR
%   integer, average heart rate in bpm
% method
%   integer, indicating the method that determined the most accurate HR
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
%   integer containing the signaling frequency in Hz
% plots
%   integer specifying to plot results (1 for true, anything else for false). default 0
% minProm
%   integer specifying MinPeakProminence in findpeaks
% minDist
%   integer specifying MinPeakDistance in findpeaks
% minHeight
%   integer specifying MinPeakHeight in findpeaks

    %set default values
    inputs ={'refHR','sig','tm','fs','plots','minProm','minDist','minHeight'};
    fs = 30;
    tm = linspace(0,300/fs,300);
    plots = 0;
    minProm = 0.05;
    minDist = 8;
    minHeight = 0.45;

    % assigning variables based off input
    for n=1:nargin
        if(~isempty(varargin{n}))
            eval([inputs{n} '=varargin{n};']);
        end
    end
    
    tic;
    if minHeight ~= 0
        [pks, locs] = findpeaks(sig, 'MinPeakProminence',minProm, 'MinPeakDistance',minDist,'MinPeakHeight',minHeight);
    else 
        [pks, locs] = findpeaks(sig, 'MinPeakProminence',minProm, 'MinPeakDistance',minDist);
    end
    toc0 = toc;
    tic;
    interval = diff(locs/fs); %convert to seconds
    HR1 = (60/mean(interval));% method 1. using STI interval method to calculate HR
    toc1 = toc + toc0;

    tic;
    HR2 = (length(pks)*6); % method 2. multiply number of systolic peaks by 6 to get BPM
    toc2 = toc + toc0;

    % method 3 - FFT
    [HR3, method_fft, f, sigfft, peak_freq, maxpeak, toc3] = detFFT(refHR, sig,fs);

    % put all HRs into array, find differences, find minimum diff and
    % return that and its index. provide times in case of duplicate
    rates = [HR1, HR2, HR3];
    times = [toc1 toc2 toc3];
    [HR, method] = compareHR(refHR,rates, times);

    if method == 3
        method = method_fft;
    end

    % plotting
    if plots == 1
        figure; 
        subplot(2,1,1);
        plot(tm,sig); hold on; plot(locs/fs, pks, 'ro')
        xlabel('Time (s)'); title(['PPG signal. Ref HR: ', num2str(refHR), ' bpm.'])
        subtitle(['HR 1: ', num2str(round(HR1)), ' bpm. HR2: ', num2str(round(HR2)),' bpm.'])
        subplot(2,1,2);
        plot(f(2:end),sigfft(2:end)); hold on;
        plot(peak_freq,maxpeak,'ro')
        xlabel('Frequency (Hz)'); 
        title('FFT of PPG Signal')
        subtitle(['HR3: ',num2str(round(HR3)), ' bpm.'])

    end

end

function [HR_fft, method_fft, f, sigfft, best_peak_freq, best_maxpeak, best_toc] = detFFT(refHR, sig,fs)
% determine the best method of calculating HR from FFT and return that
%
% Outputs:
% HR_fft
%   integer, average heart rate in bpm
% method_fft
%   string, indicating the method that determined the most accurate HR
% f
%   vector, containing frequencies. used for plotting
% sigfft
%   vector, containing the amplitudes from the FFT. used for plotting
% best_peak_freq
%   integer, containing the peak frequency. used for plotting
% best_maxpeak
%   integer, containing the amplitude of the peak frequency. used for plotting
% best_toc
%   double, containing runtime for the best method used
% 
% Parameters:
% refHR
%   integer specifying the reference HR taken from the ECG
% sig
%   vector containing the signal
% fs
%   integer containing the signaling frequency in Hz

    tic;
    N_fft = length(sig);
    toc_nfft = toc;

    tic;
    f1=((0:N_fft-1)./N_fft)*fs;
    toc0 = toc;

    % original FFT
    tic;
    og_fft = abs(fft(sig,N_fft));
    [HR1, peak_freq1, maxpeak1] = det_peak(f1,og_fft);
    toc1 = toc + toc_nfft + toc0;

    % Hann Window to increase resolution
    tic;
    hann = hanning(length(sig)).*sig;
    toc_hann = toc;

    tic;
    hann_fft = abs((fft(hann))/(N_fft/2));
    [HR2, peak_freq2, maxpeak2] = det_peak(f1,hann_fft);
    toc2 = toc + toc_nfft + toc0 + toc_hann;

    % Padd with Zeros to increase precision
    tic;
    n = 2^nextpow2(N_fft);
    f2 = fs*(0:(n-1))/n;
    toc_ze = toc;

    tic;
    zero_fft = abs((fft(sig, n))/N_fft); % Compute the FFT
    [HR3, peak_freq3, maxpeak3] = det_peak(f2,zero_fft);
    toc_3 = toc;
    toc3 = toc_3 + toc_nfft + toc_ze;

    % Hann Window and Zero padding
    tic;
    haze_fft = abs((fft((hann),n))/(N_fft/2));
    [HR4, peak_freq4, maxpeak4] = det_peak(f2,haze_fft);
    toc_4 = toc;
    toc4 = toc_4 + toc_ze + toc_hann + toc_nfft;

    % determine closest HR
    fft_rates = [HR1, HR2, HR3, HR4];
    fft_times = [toc1, toc2, toc3, toc4];
    [HR_fft, method_num] = compareHR(refHR, fft_rates, fft_times);

    switch method_num
        case 1
            % method_fft = "FFT - Regular";
            method_fft = 3;
            sigfft = og_fft;
            f = f1;
            best_peak_freq = peak_freq1;
            best_maxpeak = maxpeak1;
            best_toc = toc1;
        case 2
            % method_fft = "FFT - Hann Window";
            method_fft = 4;
            sigfft = hann_fft;
            f = f1;
            best_peak_freq = peak_freq2;
            best_maxpeak = maxpeak2;
            best_toc = toc2;
        case 3
            % method_fft = "FFT - Zero Padding";
            method_fft = 5;
            sigfft = zero_fft;
            f = f2;
            best_peak_freq = peak_freq3;
            best_maxpeak = maxpeak3;
            best_toc = toc3;
        case 4
            % method_fft = "FFT - Hann Window + Zero Padding";
            method_fft = 6;
            sigfft = haze_fft;
            f = f2;
            best_peak_freq = peak_freq4;
            best_maxpeak = maxpeak4;
            best_toc = toc4;
        otherwise
            warning('Error. Something went wrong in logic.')
    end

end

function [HR, peak_freq, maxpeak] = det_peak(f, sigfft)
% Outputs:
% HR
%   heart rate, in BPM
% peak_freq
%   integer, containing the peak frequency
% maxpeak
%   integer, containing the amplitude of the peak frequency
%
% Parameters:
% f
%   vector, containing frequencies.
% sigfft
%   vector, containing the amplitudes from the FFT.

    % limit frequency range to 0.5-4 Hz
    mask = (0.5 < f) & (f < 4);
    limit = sigfft(mask);
    f_limit = f(mask);
    
    % find peak within that range
    [peaks, loc] = findpeaks(limit);
    [maxpeak, idx] = max(peaks); %max value and its index
    peak_freq = f_limit(loc(idx));
    HR = (peak_freq*60);

end