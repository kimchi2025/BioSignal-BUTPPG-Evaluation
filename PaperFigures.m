% For plotting
% Sections are for plotting different figures

% Select a signal by commenting out the ones you do not want to look at. 
[sig,fs,tm] = rdsamp('butppg/112080/112080_PPG',1); %good quality one 
refHR = 71;
% [sig,fs,tm] = rdsamp('butppg/112002/112002_PPG',1); %bad quality one 
% refHR=64;
% [sig,fs,tm] = rdsamp('butppg/112096/112096_PPG',1); %good quality one 
% refHR=63;
% [sig,fs,tm] = rdsamp('butppg/100004/100004_PPG'); %bad quality 
% refHR=85;
% [sig,fs,tm] = rdsamp('butppg/101001/101001_PPG'); %good quality 
% refHR=67;

% creates a time vector if there is none and transposes signal vector to
% match dimensions with the newer data files
if mean(tm) == 0
    sig = sig';
    tm = linspace(0,length(sig)/fs,length(sig));
    tm = tm';
end

% preprocessing/smoothing
PPG=sig-movmean(sig, 30);

%% Repeat what authors did
% Figure 1 and in Figure 4
annotateData(sig,tm,fs);

%%
% testing different movmeans. not used in paper
figure;
T = tiledlayout(4,1,'TileSpacing','compact');

nexttile
plot(tm,sig); title("Original PPG")

nexttile
k = 25;
PPG=sig-movmean(sig, k);
plot(tm, PPG); title(['Movmean k = ', num2str(k)])

nexttile
k = 30;
PPG=sig-movmean(sig, k);
plot(tm,PPG); title(['Movmean k = ', num2str(k)])

nexttile
k = 40;
PPG=sig-movmean(sig, k);
plot(tm,PPG); title(['Movmean k = ', num2str(k)])

xlabel('time (secs)');

%% FFT Plots
% Not used in paper

N_fft = length(PPG);
f1=((0:N_fft-1)./N_fft)*fs;
n = 2^nextpow2(N_fft);
f2 = fs*(0:(n-1))/n;

figure;
T = tiledlayout(4,1,'TileSpacing','compact');

nexttile
% original FFT

og_fft = abs(fft(PPG,N_fft));
plot(f1,og_fft); title('Original Signal FFT')

nexttile
% Hann window
hann = hanning(length(PPG)).*PPG;
hann_fft = abs((fft(hann))/(N_fft/2));
plot(f1, hann_fft); title('Hann Window FFT')

nexttile
% Zero padding
zero_fft = abs((fft(PPG, n))/N_fft);
plot(f2,zero_fft); title('Zero-Padding FFT')

nexttile
% Hann window and zero padding
haze_fft = abs((fft((hann),n))/(N_fft/2));
plot(f2,haze_fft); title('Hann Window and Zero Padding FFT')

%% Comparing Filters and Normalized Data
% Figure 2 in Paper
figure;
T = tiledlayout(5,2,'TileSpacing','compact');
title(T,'PPG Signal Filters')
st = subtitle(T,'Filtering, then Normalizing and Flattening the Peaks');
st.FontSize = 10;

nexttile
plot(tm,sig); title("Original PPG")
nexttile
plot(tm,PPG); title('Signal Smoothed with Movmean')

% Butterworth
[b, a] = butter(2, [0.5 10] / (fs / 2), 'bandpass');
sig2 = filtfilt(b, a, PPG);
PPG_trend_butter=movmax(abs(sig2), 30); %relative systolic peaks
PPG_flat_butter=sig2./PPG_trend_butter;

nexttile
plot(tm,sig2); title("Butterworth Filter")
nexttile
plot(tm,PPG_flat_butter); title("Butterworth + Normalized")

% IIR filtering
PPG_filt=highpass(PPG,0.8,fs,'ImpulseResponse','iir','Steepness',0.5);%high pass filter to eliminate motion artifacts and other slow processes 
PPG_filt=lowpass(PPG_filt,3,fs,'ImpulseResponse','iir','Steepness',0.8);%low pass filter to eliminate motion artifacts and other slow processes 
PPG_trend=movmax(abs(PPG_filt), 30); %relative systolic peaks
PPG_flat=PPG_filt./PPG_trend;  %normalize and flatten

nexttile
plot(tm,PPG_filt); title("Highpass then Lowpass Filters IIR")
nexttile
plot(tm,PPG_flat); title('IIR + Normalized')

% Chebyshev
[b, a] = cheby2(4, 20, [0.5 10] / (fs / 2));
PPG_cheby = filtfilt(b, a, PPG);
PPG_trend_cheby=movmax(abs(PPG_cheby), 30); %relative systolic peaks
PPG_flat_cheby=PPG_cheby./PPG_trend_cheby;

nexttile
plot(tm,PPG_cheby); title('Chebyshev type II Filter')
nexttile
plot(tm,PPG_flat_cheby); title('Cheby + Normalized')

% IMF
IMF = emd(PPG);
PPG_IMF = sum(IMF(2:3,:),1);
PPG_trend=movmax(abs(PPG_IMF), 30); %relative systolic peaks
PPG_flat_IMF=PPG_IMF./PPG_trend; 

nexttile
plot(tm,PPG_IMF); title('EMD Filtering')
xlabel('Time (s)')
nexttile
plot(tm,PPG_flat_IMF); title('EMD + Normalized')
xlabel('Time (s)')

%% Testing best window to use
% Not used in paper. Was done for making the functions more precise.

% A few trials on 6 signals
% tallies(ties all get 1 point):
%   15: 2 points   22: 4 points   30: 2 points
% tallies (testing with 22,26,30):
%   22: 5 points   26: 3 points   30: 2 points
% could be worth doing this again with more files. but i manually count the
% points to account for ties. 22 won

window = [22 26 30];

HRs = zeros(1,length(window));
bfilts = repmat(" ",size(window));
methods = repmat(" ", size(window));
times = zeros(1,length(window));

for i=1:length(window)
    tic;
    [HRs(i), bfilt, method, PPG_filtered] = filterComparison(refHR, PPG, tm, fs, 0,0,window(i));
    times(i) = toc;
    bfilts(i) = string(bfilt);
    methods(i) = string(method);
end

[bHR, bMethod] = compareHR(refHR,HRs,times);
disp(['Ref HR: ', num2str(refHR)])
disp(['Best window: ',num2str(window(bMethod))])
disp(HRs)
disp(['Best HR: ', num2str(bHR)])
disp(['Best filter: ',num2str(bfilts(bMethod))])
disp(['Best HR claclulation method: ', num2str(methods(bMethod))])
disp(bfilts)
disp(" ")