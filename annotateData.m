function annotateData(sig,tm,fs)
% This function will recreate the annotating software similar to the image
% shown in Figure 3 of the paper
%
% Reference: 
% Nemcova A et al., Brno University of Technology Smartphone PPG Database (BUT PPG): 
% Annotated Dataset for PPG Quality Assessment and Heart Rate Estimation.
% BioMed Research International. 2021 Sep 7;2021. https://doi.org/10.1155/2021/3453007.
%
% Parameters:
% sig
%   vector containing the signal
% tm
%   vector containing the time in SECONDS
% fs
%   integer containing the signaling frequency in Hz

    % For some files, tm is 0 and must be manually calculated
    % disp(mean(tm)==0) %for debugging purposes
    if mean(tm) == 0
        tm = linspace(0,length(sig)/fs,length(sig));
    end

    % Wavelet Decomposition and Projection
    % Code based off the matlab wavedec function example 3
    
    level = 6; % number of decompositions
    wv = 'db4'; % wavelet to use
    [C, L] = wavedec(sig,level,wv);
    
    mra = zeros(level+1,numel(sig)); %preallocate
    % obtain projections onto wavelet detail subspaces
    for k=1:level
        mra(k,:) = wrcoef("d",C,L,wv,k);
    end
    mra(end,:) = wrcoef("a",C,L,wv,level); %approximation subspace

    
    % Determination of HR from FFT
    % calculate FFT
    N_fft = length(sig);
    sigfft = abs(fft(sig,N_fft));
    f=((0:N_fft-1)./N_fft)*fs;
    
    % limit frequency range to 0.5-4 Hz
    mask = (0.5 < f) & (f < 4);
    limit = sigfft(mask);
    f_limit = f(mask);
    
    % find peak within that range
    [peaks, loc] = findpeaks(limit);
    [maxpeak, idx] = max(peaks); %max value and its index
    peak_freq = f_limit(loc(idx));
    HR_fft = peak_freq*60;
    
    % disp(HR_fft) %for debugging

    % Plotting everything
    figure;
    T = tiledlayout(6,2,'TileSpacing','compact');
    
    nexttile
    plot(tm,sig); title('PPG signal');
    nexttile([2 1])
    plot(tm,sig);  title('PPG'); xlabel('Time (s)')
    
    t = (1:length(sig))'./fs;
    j = 1;
    k = 6;
    for i = 2:k
        [pks, locs] = findpeaks(mra(i,:), 'MinPeakProminence',0.01, 'MinPeakDistance',8);
        interval = diff(locs/fs); %convert to seconds
        med = median(60./interval);
        HR = 60/mean(interval);% using RR interval method to calculate HR
        
        j = j+2;
        nexttile(j) 
        plot(t,mra(i,:)); hold on;
        plot(locs/fs, pks, 'ro')
        th = title(['SWT Band ', num2str(i), '. Mean HR: ', num2str(round(HR)), ...
            ' bpm. Median HR: ', num2str(round(med)), ' bpm.']);
        th.FontSize = 7;
        %xlim([0 10]);
    end
    xlabel('Time (s)')
    
    nexttile(8,[3 1])
    % skipping first data point to avoid skew
    plot(f(2:end),sigfft(2:end)); hold on;
    plot(peak_freq,maxpeak,'ro')
    xlabel('Frequency (Hz)'); 
    title(['FFT of PPG Signal. HR: ',num2str(round(HR_fft)), ' bpm.'])

end






