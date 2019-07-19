function [csc_bandpower, csc_filt, filt_times, window_duration] = bandpower_boxcar(csc_mtx, timebounds, column, filterbounds)
%filters to 150hz - 250hz (rippleband)

%some input options
low =  filterbounds(1); %hz
hi = filterbounds(2); %hz

%check timebounds
if isempty(timebounds)
    timebounds(1) = csc_mtx(1,1);
    timebounds(2) = csc_mtx(end,1);
end

%samplerate
samplerate =  1600;%size(csc_mtx,1)/csc_mtx(end,1);

%execute filter
csc_filt = bpfilt(csc_mtx(csc_mtx(:,1)>=timebounds(1) & csc_mtx(:,1)<=timebounds(2), column), low, hi, samplerate, 0);
filt_times = csc_mtx(csc_mtx(:,1)>=timebounds(1) & csc_mtx(:,1)<=timebounds(2),1);

%calculate rippleband power over time
window_duration = 0.05; %seconds
windowsize = ceil(samplerate*window_duration); %rows
if rem(windowsize,2)==1
    windowsize = windowsize+1;
end
increment = ceil(windowsize/10);
calc_vect = 1+windowsize/2:increment:length(csc_filt)-windowsize/2;

%preallocate
csc_bandpower = nan(length(csc_filt), 1);
power_time = nan(size(csc_bandpower));

%loop calculate bandpower of time
for i = calc_vect
    csc_bandpower(i) = bandpower(csc_filt(i-windowsize/2:i+windowsize/2));
    power_time(i) = mean(filt_times(i-windowsize/2:i+windowsize/2));
end


%interp at all time points
csc_bandpower = interp1(power_time(~isnan(power_time)),csc_bandpower(~isnan(csc_bandpower)),filt_times);

%zscore bandpower
%csc_bandpower = (csc_bandpower-nanmean(csc_bandpower))./nanstd(csc_bandpower);


















