function [bc_rates, times] = boxcar_rates(spiketimes, timerange, bc_duration, slide_duration, varargin)
%takes spike times and computes rates using a sliding boxcar

%input number of trials vector (size(spiketimes))
if ~isempty(varargin)
    num_trials = varargin{1};
else
    num_trials = ones(size(spiketimes));
end

%times of interest
times_lo = min(timerange):slide_duration:max(timerange)-bc_duration;
times_hi = times_lo + bc_duration;
times = mean([times_lo;times_hi]);

%preallocate
bc_counts = nan(length(spiketimes), length(times_lo));


%iterate through times
for itime = 1:length(times_lo)

    %each cluster
    for iclust = 1:length(spiketimes)
        
        %count spikes during times of interest
        bc_counts(iclust,itime) = sum(spiketimes{iclust}>=times_lo(itime) & spiketimes{iclust}<times_hi(itime));                   
    end
end

%counts to rates
bc_rates = bc_counts./bc_duration;

%divide by number of trials
bc_rates = bc_rates./repmat(num_trials, 1, size(bc_rates,2));

%plot normalized rates
%figure; imagesc(times, 1:length(spiketimes), norm_mtx(bc_rates')')

end