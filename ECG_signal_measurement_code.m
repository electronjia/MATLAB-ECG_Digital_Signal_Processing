%  BE 3344       Lab 1 - ECG         03/28/2022

clear
close all
clc

%% Resistor used: 4.7K Ohm Resistor
%Import and plot raw ECG data
data = readmatrix('4.7k_ohmECG.CSV'); % imports oscilloscope obtained CSV data into a matrix
A = data(:,2); % extracts just the voltage values (col #2)
figure
subplot(2,1,1)
plot(A)
title('Raw ECG Signal - 4.7K Ohm Resistor')
xlabel('Time (ms)')
ylabel('Voltage (V)')

% Generate frequency magnitude and phase spectra using Fast Fourier Transform
fs = 50000; % sampling frequency (Hz)
L = length(A); % length of signal
y = fft(A); 
ds = abs(y/L); % double-sided amplitude spectrum
ss = ds(1:(L/2)+1);
ss(2:end-1) = 2*ss(2:end-1); % single-sided amplitude spectrum
f = (0:L/2)*fs/L; % frequency axis in Hz
subplot(2,1,2)
plot(f,ss)
xlabel('Frequency (Hz)')
ylabel('Amplitude')
title('Single-Sided Amplitude Spectrum - 4.7K Ohm Resistor')

% Finite Impulse Response (FIR) lowpass filter 
fn = fs/2; % Nyquist frequency
fc = 100; % cut-off frequency
n = 50; % filter order
b = fir1(n,fc/fn);
z = filtfilt(b,1,A);
figure
plot(z)
title('FIR Filtered (LP) ECG Signal - 4.7K Ohm Resistor')
xlabel("Time (ms)")
ylabel("Voltage (V)")


%{
% Finite Impulse Response (FIR) bandpass filter
fc_L = 50; % lower bound cut-off frequency (Hz)
fc_H = 400; % upper bound cut-off frequency (Hz)
b1 = fir1(n,[(fc_L/fn) (fc_H/fn)]);
z1 = filtfilt(b1,1,A); % zero-phase filtering
figure
plot(z1)
title('FIR Filtered (BP) ECG Signal - 4.7K Ohm Resistor')
%}

% Detect R peaks
[pk,x] = findpeaks(z,'MinPeakHeight',0.15); % finds all peaks w/ amplitudes >=0.16 V; stores their indices in vector x
hold on
plot(x,pk,'ro')


% Computing heart rate
% Method 1
hb = length(pk); % total # of R peaks = # of heart beats
ts = L/fs; % total # of samples/sampling rate (samples/sec) = time in seconds
tm = ts/60; % total time in minutes
bpm1 = hb/tm % heart rate (beats per minute)
% Method 2
ts = (x(2)-x(1))/fs; % time interval b/w 2 adjacent R peaks (sec)
tm = ts/60; % time interval in minutes
bpm2 = length(pk)/tm % heart rate (beats per minute)

% Method 3 
trr=x(2)-x(1)
bpm3=1/(trr*10^(-3))*60

% Frequency and phase response of lowpass FIR filter
figure
subplot(2,1,1)
freqz(b) % frequency response
title('Frequency Response - 4.7K Ohm Resistor')
subplot(2,1,2)
phasez(b) % phase response
title('Phase Response - 4.7K Ohm Resistor')

%% Resistor value: 10K Ohm Resistor

clear
clc
close all
% Import and plot raw ECG data
data = readmatrix('10k_ohmECG.CSV'); % imports data from Excel into a matrix
A = data(:,2); % extracts just the voltage values (col #2)
figure
subplot(2,1,1)
plot(A)
title('Raw ECG Signal - 10K Ohm Resistor')
xlabel('Time (ms)')
ylabel('Voltage (V)')

% Generate frequency magnitude spectrum using FFT
fs = 50000; % sampling frequency (Hz)
L = length(A); % length of signal
y = fft(A); 
ds = abs(y/L); % double-sided amplitude spectrum
ss = ds(1:(L/2)+1);
ss(2:end-1) = 2*ss(2:end-1); % single-sided amplitude spectrum
f = (0:L/2)*fs/L; % frequency axis in Hz
subplot(2,1,2)
plot(f,ss)
xlabel('Frequency (Hz)')
ylabel('Amplitude')
title('Single-Sided Amplitude Spectrum - 10K Ohm Resistor')

% FIR lowpass filter 
fn = fs/2; % Nyquist frequency
fc = 310; % cut-off frequency
n = 30; % filter order
b = fir1(n,fc/fn);
z = filtfilt(b,1,A); % zero-phase filtering
figure
plot(z)
title('FIR Filtered (LP) ECG Signal - 10K Ohm Resistor')
xlabel("Time (ms)")
ylabel("Voltage (V)")

% Detect R peaks
[pk,x] = findpeaks(z,'MinPeakHeight',0.16); % finds all R peaks and stores their indices in vector x
hold on
plot(x,pk,'ro')
% Method 3 
trr=x(2)-x(1)
bpm3=1/(trr*10^(-3))*60

% Frequency and phase response of lowpass FIR filter
figure
subplot(2,1,1)
freqz(b) % frequency response
title('Frequency Response - 10K Ohm Resistor')
subplot(2,1,2)
phasez(b) % phase response
title('Phase Response - 10K Ohm Resistor')

%% 150 Ohm Resistor

clear
clc
close all
% Import and plot raw ECG data
data = readmatrix('150_ohmECG.CSV');
A = data(:,2); 
figure
subplot(2,1,1)
plot(A)
title('Raw ECG Signal - 150 Ohm Resistor')
xlabel('Time (ms)')
ylabel('Voltage (V)')

% Generate frequency magnitude spectrum using FFT
fs = 50000; % sampling frequency (Hz)
L = length(A); % length of signal
y = fft(A); 
ds = abs(y/L); % double-sided amplitude spectrum
ss = ds(1:(L/2)+1);
ss(2:end-1) = 2*ss(2:end-1); % single-sided amplitude spectrum
f = (0:L/2)*fs/L; % frequency axis in Hz
subplot(2,1,2)
plot(f,ss)
xlabel('Frequency (Hz)')
ylabel('Amplitude')
title('Single-Sided Amplitude Spectrum - 150 Ohm Resistor')

% FIR lowpass filter 
fn = fs/2; % Nyquist frequency
fc = 250; % cut-off frequency
n = 30; % filter order
b = fir1(n,fc/fn);
z = filtfilt(b,1,A);
figure
plot(z)
title('FIR Filtered (LP) ECG Signal - 150 Ohm Resistor')
xlabel("Time (ms)")
ylabel("Voltage (V)")

% Detect R peaks
[pk,x] = findpeaks(z,'MinPeakHeight',0.2); % finds all R peaks and stores their indices in vector x
hold on
plot(x,pk,'ro')
% Method 3 
trr=x(2)-x(1)
bpm3=1/(trr*10^(-3))*60

% Frequency and phase response of lowpass FIR filter
figure
subplot(2,1,1)
freqz(b) % frequency response
title('Frequency Response - 150 Ohm Resistor')
subplot(2,1,2)
phasez(b) % phase response
title('Phase Response - 150 Ohm Resistor')

%% 470 Ohm Resistor
clear
clc
close all

% Import and plot raw ECG data
data = readmatrix('470_ohmECG.CSV');
A = data(:,2); 
figure
subplot(2,1,1)
plot(A)
title('Raw ECG Signal - 470 Ohm Resistor')
xlabel('Time (ms)')
ylabel('Voltage (V)')

% Generate frequency magnitude spectrum using FFT
fs = 50000; % sampling frequency (Hz)
L = length(A); % length of signal
y = fft(A); 
ds = abs(y/L); % double-sided amplitude spectrum
ss = ds(1:(L/2)+1);
ss(2:end-1) = 2*ss(2:end-1); % single-sided amplitude spectrum
f = (0:L/2)*fs/L; % frequency axis in Hz
subplot(2,1,2)
plot(f,ss)
xlabel('Frequency (Hz)')
ylabel('Amplitude')
title('Single-Sided Amplitude Spectrum - 470 Ohm Resistor')

% FIR lowpass filter 
fn = fs/2; % Nyquist frequency
fc = 250; % cut-off frequency
n = 30; % filter order
b = fir1(n,fc/fn);
z = filtfilt(b,1,A);
figure
plot(z)
title('FIR Filtered (LP) ECG Signal - 470 Ohm Resistor')
xlabel("Time (ms)")
ylabel("Voltage (V)")

% Detect R peaks
[pk,x] = findpeaks(z,'MinPeakHeight',0.15); % finds all R peaks and stores their indices in vector x
hold on
plot(x,pk,'ro')
% Method 3 
trr=x(2)-x(1)
bpm3=1/(trr*10^(-3))*60

% Frequency and phase response of lowpass FIR filter
figure
subplot(2,1,1)
freqz(b) % frequency response
title('Frequency Response - 470 Ohm Resistor')
subplot(2,1,2)
phasez(b) % phase response
title('Phase Response - 470 Ohm ')
