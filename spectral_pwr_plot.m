function [pxx,f] = spectral_pwr_plot(csc_mtx, varargin)
%function [pxx,f] = spectral_pwr_plot(csc_mtx, filterbounds, freq_rng, timebounds)
%line plot of power over all frequencies (freq shown in log)


%handle inputs
switch nargin>1
    case nargin == 2
        filterbounds = varargin{1};
        freq_rng = filterbounds;
        timebounds = [csc_mtx(1,1) csc_mtx(end,1)];
    case nargin == 3
        filterbounds = varargin{1};
        freq_rng = varargin{2};
        timebounds = [csc_mtx(1,1) csc_mtx(end,1)];
    case nargin == 4
        filterbounds = varargin{1};
        freq_rng = varargin{2};
        timebounds = varargin{3};
    otherwise
        filterbounds = [1 300];
        freq_rng = filterbounds;
        timebounds = [csc_mtx(1,1) csc_mtx(end,1)];
end
if isempty(freq_rng)
    freq_rng = [1 300];
end


%filterbounds = [4 300];
%freq_rng = filterbounds;


%sampling rate
%{
plausible_sample_rates = 32000./(1:1000);
samplerate = round(size(csc_mtx,1)/csc_mtx(end,1));
samplerate = plausible_sample_rates(abs(plausible_sample_rates-samplerate)...
    ==min(abs(plausible_sample_rates-samplerate)));
%}
samplerate = 1600;

%avoid nyquist-freq-related problems
if filterbounds(2) > floor(samplerate/2)
    filterbounds(2) = floor(samplerate/2) -1;
end
filterbounds = filterbounds

%timewindow
%number of points in 50ms
%tw = ceil(samplerate/20);
%number of points in 1000ms (1s)
tw = samplerate/20;

%timebound index
tbi = csc_mtx(:,1)>=timebounds(1) & csc_mtx(:,1)<=timebounds(2);

hold on
for icolumn = 2:size(csc_mtx,2)
    
    %filter
    csc_filt = bpfilt(csc_mtx(tbi, icolumn), filterbounds(1), filterbounds(2), samplerate, 0);

    %power and frequency
    %pwelch(csc_filt, tw, [], [], samplerate);
    [pxx,f] = pwelch(csc_filt, tw, [], [], samplerate, 'power');
    
    
    %limit frequency range for plotting
    f1 = find(abs(f-freq_rng(1))==min(abs(f-freq_rng(1))), 1, 'first');
    f2 = find(abs(f-freq_rng(2))==min(abs(f-freq_rng(2))), 1, 'last');

    %plot
    plot(f(f1:f2), 10*log10(pxx(f1:f2)))
    
end

%aesthetics
xlabel('Frequency (Hz)')
ylabel('Power (mV^2 / Hz)') %check these units!!!



%fourier transform


end