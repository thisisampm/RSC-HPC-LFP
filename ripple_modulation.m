function mod_scores = ripple_modulation(datamtx, csc_mtx, clusters, RSC_TTs, HPC_TTs, col)
% computes a score for each cell indicating how strongly the firing rate
% changes in response to the ripple event

%hard input
time_before = -0.50; %s
time_after = 0.60; %s
sesh_startend = datamtx([1 end], 1);
ripple_duration = 0.07; %s

%preallocate
cellevent_times = cell(length(clusters),1);

%region of interest column
%col % HPC=12, RSC=14, control (1s delay) = 0
if col == 0
    ctrl = 1;
    col = 12;
else
    ctrl = 0;
end

for irip = min(datamtx(datamtx(:,12)>0,col)):max(datamtx(:,col))
    timerng = [min(datamtx(datamtx(:,col)==irip,1))+time_before min(datamtx(datamtx(:,col)==irip,1))+time_after];
    if ctrl == 1
        timerng = timerng + 1;
    end

     if ~isempty(timerng) && timerng(1)>sesh_startend(1) && timerng(2)<sesh_startend(2)
        [~, ~, ~, cellevent_times_local] = lfp_plot(csc_mtx, timerng, RSC_TTs, HPC_TTs, datamtx, clusters, 0);
        
        %cell events
        for ic = 1:length(cellevent_times_local)
            %normalize times to ripple duration
            cellevent_times_local{ic} = (cellevent_times_local{ic}-timerng(1))./diff(timerng);
            %load
            cellevent_times{ic} = sort([cellevent_times{ic}; cellevent_times_local{ic}]);
        end
    end
end

%binning for rates
tw_length = time_after - time_before; %seconds
bin_duration = 0.005; %seconds
bins = ceil(tw_length/bin_duration);
xedges = linspace(0,1+realmin,bins+1);

%compute bin firing rates
rm = nan(size(clusters,1), bins);
count = 0;
for icluster = 1:size(clusters,1)
    count = count+1;
    cellevent_times_plot = cellevent_times{icluster};
    spikecounts = histcounts(cellevent_times_plot, xedges);
    spikecounts = smooth(spikecounts, 3);
    rm(count,:) = spikecounts; %not normalized
end

%find bin immediately after 0s
xedges_time = linspace(time_before,time_after,bins+1);
rip_bin_idx = abs(xedges_time)==min(abs(xedges_time));

%slide ripple-length TW across ratewindow
ripple_binwidth = floor(ripple_duration./bin_duration);
if ripple_binwidth==0; ripple_binwidth = 1; end
    ripple_rates = nan(size(rm,1), size(rm,2)-ripple_binwidth);
for ibin = 1:size(rm,2)-ripple_binwidth
    ripple_rates(:, ibin) = sum(rm(:,ibin:ibin+ripple_binwidth), 2);
end

%zscore all bins
ripple_rates = zscore_mtx(ripple_rates')';

%code clusters in terms of regions
region_ids = zeros(size(clusters(:,1)));
%RSCtts = 1
region_ids(ismember(floor(clusters(:,1)), RSC_TTs))=1;
%HPCtts = 2
region_ids(ismember(floor(clusters(:,1)), HPC_TTs))=2;

%ouput zscores of ripple time window
mod_scores = [region_ids ripple_rates(:, rip_bin_idx)];





end