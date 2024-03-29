function [signal_filt] = bpfilt(signal, f1, f2, fs, isplot)
%% Bandpass filtering
%
% Syntax:
%   y = bpfilt(signal, f1, f2, [options])
%
% Description:
%   This function performs bandpass filtering of a time series 
%   with rectangle window.
%
% Input Arguments:
%   signal 	- a column vector of time series.
%   f1 		- the lower bound of frequencies (in Hz).
%   f2 		- the upper bound of frequencies (in Hz).
%
% Options:
%   fs      - the sampling frequency in Hz. Default is 1 Hz.
%   isplot  - whether to produce plots.
%
% Output Arguments:
%   y 		- the filtered time series.
%
% Examples:
%   fs = 100;
%   t  = 1:1/fs:10;
%   x  = sin(t);
%   y  = bpfilt(x,20,30);

%% getting options
if nargin < 4 || isempty(fs)
	fs = 1;
end

if nargin < 5 || isempty(isplot)
	isplot = 1;
end

%% define variables
if isrow(signal)
    signal = signal';
end
N  = length(signal);
dF = fs/N;
f  = (-fs/2:dF:fs/2-dF)';

%% Band-Pass Filter:
if isempty(f1) || f1==-Inf
    BPF = (abs(f) < f2);
elseif isempty(f2) || f2==Inf
    BPF = (f1 < abs(f));
else
    BPF = ((f1 < abs(f)) & (abs(f) < f2));
end


%% Power spectrum of the original signal
signal 	 = signal-mean(signal);
spektrum = fftshift(fft(signal));%/N;

if isplot
    figure;
    subplot(2,1,1);
    plot(f,abs(spektrum));
    title('Power spectrum of the original signal');
end

%% Power spectrum of the band-pass filtered signal
spektrum = BPF.*spektrum;
%power_filt = spektrum;
if isplot
    subplot(2,1,2);
    plot(f,abs(spektrum));
    title(sprintf('Power spectrum of the band-pass filtered signal in (%.3f, %.3f) Hz',f1,f2));
end

%% The band-pass filtered time series
signal_filt = ifft(ifftshift(spektrum)); %inverse ifft
signal_filt = real(signal_filt);

if isplot
    time = 1:length(signal_filt);
	
    figure;
    subplot(2,1,1);
    plot(time,signal);
    title('The original time series');
    subplot(2,1,2);
    plot(time,signal_filt);
    title(sprintf('The band-pass filtered time series in (%.3f, %.3f) Hz',f1,f2));
end