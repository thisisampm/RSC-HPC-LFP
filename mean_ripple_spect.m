function [spectrogram_mtx, cell_rms, cellevent_times, frequencies, ripple_power] = mean_ripple_spect(datamtx, csc_mtx, clusters, time_beforeafter, RSC_TTs, HPC_TTs, varargin)
%plot the mean spectrogram for all ripples (aligned to ripple start)


if nargin==7
   ploton = varargin{1};
else
   ploton = 1;
end
rsc_hpc_tts{1} = RSC_TTs; rsc_hpc_tts{2} = HPC_TTs;
sesh_startend = datamtx([1 end], 1);
spectrogram_mtx = [];
all_ripple_power_rsc = [];
all_ripple_power_hpc = [];
cellevent_times = cell(1, length(clusters));
time_before = time_beforeafter(1);
time_after = time_beforeafter(2);

for irip = unique(datamtx(datamtx(:,12)>0,12))'
    timerng = [min(datamtx(datamtx(:,12)==irip,1))+time_before min(datamtx(datamtx(:,12)==irip,1))+time_after];

    if ~isempty(timerng) && timerng(1)>sesh_startend(1) && timerng(2)<sesh_startend(2)
        [spectrogram_mtx_local, frequencies, ~, cellevent_times_local] = lfp_plot(csc_mtx, timerng, RSC_TTs, HPC_TTs, datamtx, clusters, 0); 

        spectrogram_mtx = cat(4, spectrogram_mtx, spectrogram_mtx_local);

        %ripple power
        ripple_power_rsc = datamtx(datamtx(:,1)>=timerng(1) & datamtx(:,1)<=timerng(2) & datamtx(:,8)==0, 9)';
        ripple_power_hpc = datamtx(datamtx(:,1)>=timerng(1) & datamtx(:,1)<=timerng(2) & datamtx(:,8)==0, 10)';

        ripple_power_rsc = interp1(linspace(1, size(spectrogram_mtx_local,2), length(ripple_power_rsc)), ripple_power_rsc, 1:size(spectrogram_mtx_local,2));
        ripple_power_hpc = interp1(linspace(1, size(spectrogram_mtx_local,2), length(ripple_power_hpc)), ripple_power_hpc, 1:size(spectrogram_mtx_local,2));

        all_ripple_power_rsc = [all_ripple_power_rsc; ripple_power_rsc]; 
        all_ripple_power_hpc = [all_ripple_power_hpc; ripple_power_hpc]; 

        %cell events
        for ic = 1:length(cellevent_times_local)
            %normalize times to ripple duration
            cellevent_times_local{ic} = (cellevent_times_local{ic}-timerng(1))./diff(timerng);
            %load
            cellevent_times{ic} = sort([cellevent_times{ic}; cellevent_times_local{ic}]);
        end
    end
end
ripple_power = [mean(all_ripple_power_rsc); mean(all_ripple_power_hpc)];

size(spectrogram_mtx)

%average
spectrogram_mtx = mean(spectrogram_mtx,4);

if ploton==1
    %spectrogram subplots
    figure
    for i = 1:size(spectrogram_mtx,3)
        subplot(size(spectrogram_mtx,3),1,i)

        %prepare output
        s_mtx = spectrogram_mtx(:,:,i).*frequencies;
        s_mtx = 10*log10(s_mtx);
        %s_mtx = s_mtx.*frequencies;

        imagesc(s_mtx)
        set(gca,'Ydir','normal')
        caxis([60 95])
        colorbar
        ylim([1 find(frequencies>300,1,'first')])

        %labels
        newyticklabels = 0:50:300;
        newyticks = interp1(frequencies, 1:size(spectrogram_mtx,1), newyticklabels);
        yticks(newyticks)
        yticklabels(newyticklabels)
        ylabel('Frequency (Hz)')
        holdxlim = xlim;
        zero_sec_bins = interp1(linspace(time_beforeafter(1), time_beforeafter(2), size(s_mtx,2)), 1:size(s_mtx,2), 0);
        xticks([1 zero_sec_bins size(s_mtx,2)]); 
        xticklabels([time_beforeafter(1) 0 time_beforeafter(2)])
        xlabel('Time (s)')
        hold on; plot([1 1].*zero_sec_bins, ylim, 'r-')

        %aesthetics
        set(gca,'TickLength',[0, 0]); box off

    end
end

    %binning for rates
    tw_length = time_after - time_before; %seconds
    bin_duration = 0.01; %seconds
    bins = ceil(tw_length/bin_duration);
    xedges = linspace(0-realmin,1+realmin,bins+1);

    if ploton==1
        figure
    end
    cell_events_out = cell(1, size(spectrogram_mtx,3));
    cell_rms = cell(1, size(spectrogram_mtx,3));
    
    size(spectrogram_mtx)
    
    for isubp = 1:size(spectrogram_mtx,3)
        if ploton==1
            subplot(2,1,isubp); hold on
        end
        
        rsc_hpc_tts
        
        isubp
        rsc_hpc_tts{isubp}
        
        %plot cell events
        cluster_assigns = find(ismember(floor(clusters), rsc_hpc_tts{isubp}))';
        rm = nan(length(cluster_assigns), bins);
        count = 0;
        for icluster = cluster_assigns
            count = count+1;
            cellevent_times_plot = cellevent_times{icluster};
            spikecounts = histcounts(cellevent_times_plot, xedges);
            %spikecounts = smooth(spikecounts, 3);
            norm_rates = (spikecounts-min(spikecounts))./max((spikecounts-min(spikecounts)));
            rm(count,:) = norm_rates;
            %rm(count,:) = spikecounts; %not normalized
        end

        if ploton==1
        if ~isempty(rm)

            %sort rm
            %[~,maxcols] = max(rm,[],2);
            %[~,sort_idx] = sort(maxcols);    
            %sort_idx = flipud(sort_idx);

            %plot
            %imagesc(rm(sort_idx,:)); caxis([0 1])
            imagesc(rm); caxis([0 1])

            %label
            ylim([0.5 length(cluster_assigns)+0.5])
            xlim([0.5 size(rm,2)+0.5])
            set(gca,'TickLength',[0, 0]); box off
            set(gca,'Ydir','reverse')
            zero_sec_bins = interp1(linspace(time_beforeafter(1), time_beforeafter(2), size(rm,2)), 1:size(rm,2), 0);
            xticks([1 zero_sec_bins size(rm,2)]); 
            xticklabels([time_beforeafter(1) 0 time_beforeafter(2)])
                hold on; plot([1 1].*zero_sec_bins, ylim, 'r-')
            
            yticks(1:length(cluster_assigns))
            relevant_cells = clusters(cluster_assigns);
            yticklabels(relevant_cells)
            xlabel('Time (s)')
            ylabel('Cell IDs')
            
        end
        end

        %load out
        cell_events_out{isubp} = rm;        
        cell_rms{isubp} = rm;
    end
    
    
end

    

