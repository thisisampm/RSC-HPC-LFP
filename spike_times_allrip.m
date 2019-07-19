function [spiketimes, num_trials] = spike_times_allrip(datamtx, clusters, timebounds, varargin)
%identify spike times around each ripple event


%input ripple column
ctrl = 0;
if ~isempty(varargin)
    ripple_col = varargin{1};
    if ripple_col == 0
        ctrl = 1;
        ripple_col = 12;
    end
else
    ripple_col = 12;
end

%preallocate
spiketimes = cell(size(clusters,1), 1);

%ripples
ripple_idx = datamtx(:,ripple_col)>0 & datamtx(:,11)==1 & ismember(datamtx(:,7), [1 5]);
ripple_events = unique(datamtx(ripple_idx, ripple_col))';
num_trials = length(ripple_events);
if num_trials == 0
    error('No ripple events detected')
end

%for each ripple event
for irip = ripple_events
    rip_start = min(datamtx(datamtx(:, ripple_col)==irip,1)); 
    if ctrl==1
        rip_start = rip_start+1;
    end
    rip_timebounds = [rip_start+timebounds(1) rip_start+timebounds(2)];
    
    %spike times for each cells
    [~, st_norm] = spike_times(datamtx, clusters, rip_timebounds);
    
    %load each cell
    for iclust = 1:size(clusters,1)
        spiketimes{iclust} = [spiketimes{iclust}; st_norm{iclust}];
    end
    
end

%set to timebounds
%load each cell
for iclust = 1:size(clusters,1)
    spiketimes{iclust} = (spiketimes{iclust}.*diff(timebounds))+timebounds(1);
    spiketimes{iclust} = sort(spiketimes{iclust});
end

%plot
%{
downsamp = 5;
figure; hold on
for iclust = 1:size(clusters,1)
    plot([spiketimes{iclust}(1:downsamp:end) spiketimes{iclust}(1:downsamp:end)], [iclust-0.5 iclust+0.5], 'k');
end
xlim(timebounds)
set(gca,'TickLength',[0, 0]); box off;
hold on; plot([0 0], ylim, 'r-')

%}
end