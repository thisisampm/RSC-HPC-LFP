function [spectrogram_mtx, frequencies, times, cellevent_times] = lfp_plot(csc_mtx, timebounds, RSC_TTs, HPC_TTs, varargin)
%plots lfp signal during timebounds [lo hi] found in timestamps
%columns indicate which columns of csc_mtx should be plotted
%lfp cell contains hpc open, hpc ripple, rsc open, rsc ripple
%
%lfp_plot(csc_mtx, [2 4 11 13], [1606 1609]+20)

%spectrogram_mtx preallocate
spectrogram_mtx = [];

%can plot spikes
if nargin == 6
    datamtx = varargin{1};
    clusters = varargin{2};
    plot_on = 1;
elseif nargin == 7
    datamtx = varargin{1};
    clusters = varargin{2};
    plot_on = varargin{3};
end

unique_chans = size(csc_mtx,2)-1;
subplot_idx = (repmat(1:unique_chans, 3, 1) + [0;2;4]);
subplot_idx_end = subplot_idx(end,:);
subplot_idx = subplot_idx(1:end-1,:);

spectrogram_idx = csc_mtx(:,1)>timebounds(1) & csc_mtx(:,1)<timebounds(2); 
spectr_timewindow = 0.1; %s
timebound_rng = range(timebounds);

%figure
if plot_on == 1
    figure;
end

%spectogram inputs
plausible_sample_rates = 32000./(1:1000);
samplerate = round(size(csc_mtx,1)/csc_mtx(end,1));
samplerate = plausible_sample_rates(abs(plausible_sample_rates-samplerate)==min(abs(plausible_sample_rates-samplerate)));
if plot_on ==0
    binsize = round(((spectr_timewindow/timebound_rng)*(timebound_rng*samplerate))/10);
elseif plot_on==1
    binsize = round(((spectr_timewindow/timebound_rng)*(timebound_rng*samplerate))/4);
end
num_y_bins = floor(samplerate/4);

cols = [nan 14 12];
count = 1;
for i = 1:unique_chans

    %data
    time_idx = csc_mtx(:,1)>timebounds(1) & csc_mtx(:,1)<timebounds(2);
    times = csc_mtx(time_idx, 1);
    csc = csc_mtx(time_idx, i+1);
      
    %lineplot
    if plot_on == 1
        subplot(unique_chans+1, 2, subplot_idx(count))
        hold on
        count = count+1;
        plot(times, csc)
        set(gca, 'TickLength', [0 0])
        xlim(timebounds)
        ylim([-5 5].*10^4)
        colorbar
        xtick_hold = xticks;
        xticklab_hold = xticklabels;
        
        
        %red hpc ripple line
        rip_idx = csc_dmidx(csc_mtx, datamtx, datamtx(:,cols(i+1))>0);
        %times_select = find(rip_idx(time_idx)==1);
        times_deselect = find(rip_idx(time_idx)~=1);
        xplot_rip = times; xplot_rip(times_deselect) = nan;
        yplot_rip = csc; yplot_rip(times_deselect) = nan;
        
        plot(xplot_rip, yplot_rip, 'r')
        
    end
       
    %spectogram plot
    if plot_on == 1
        subplot(unique_chans+1, 2, subplot_idx(count)) 
        count = count+1;
    end

    %spectrogram
    if plot_on == 1
        [~, frequencies, times, spect_mtx_chan] = spectrogram(csc_mtx(spectrogram_idx,i+1), binsize, floor(binsize/3), num_y_bins, samplerate, 'yaxis', 'power');

        s_mtx = spect_mtx_chan.*frequencies;
        s_mtx = 10*log10(s_mtx);
        imagesc(s_mtx)
        set(gca,'Ydir','normal')
        caxis([60 95])
        colorbar
        ylim([1 find(frequencies>300,1,'first')])
        set(gca, 'TickLength', [0 0])
        hold_ylim = ylim;
        xtick_hold = (xtick_hold-timebounds(1))./(timebounds(2)-timebounds(1));
        xticks(xtick_hold.*size(s_mtx,2))
        xticklabels(xticklab_hold)
        set(gcf, 'Position', [157         197        1210         677])
    
    elseif plot_on == 0

        [~, frequencies, times, spect_mtx_chan] = spectrogram(csc_mtx(spectrogram_idx,i+1), binsize, floor(binsize/3), num_y_bins, samplerate, 'yaxis', 'power');
        
        
    end
    
    spectrogram_mtx = cat(3, spectrogram_mtx, spect_mtx_chan);
        
end

if plot_on == 1
   
    %plot spike events
    subplot(unique_chans+1, 2, subplot_idx_end(1)); hold on
    cellevent_plot(datamtx, clusters(ismember(floor(clusters), RSC_TTs),1), timebounds, 1);
    colorbar
    
    subplot(unique_chans+1, 2, subplot_idx_end(2)); hold on
    cellevent_plot(datamtx, clusters(ismember(floor(clusters), HPC_TTs),1), timebounds, 1);
    colorbar

else
    cellevent_times = cellevent_plot(datamtx, clusters, timebounds, 0);
end


if plot_on == 1
   figure
   subplot(1,2,1); hold on
   timerng_idx = datamtx(:,1)>timebounds(1) & datamtx(:,1)<timebounds(2);
   plot(datamtx(timerng_idx,2), datamtx(timerng_idx,3), '.', 'markersize', 20)
   maze_outline
   title position
   
   subplot(1,2,2); hold on
   plot(datamtx(timerng_idx,1), datamtx(timerng_idx,5))
   plot(xlim, repmat(max(datamtx(datamtx(:,11)==1,5)), 1, 2), 'k--')
   xlim(timebounds)
   ylim([0 .2])
   set(gca,'TickLength',[0, 0]); box off
   title velocity
   
   set(gcf, 'Position', [1431         451         944         337])
end

if plot_on ==1
    [~, rsc_ripplefilt] = bandpower_boxcar(csc_mtx, [0 inf], 2, [150 250]);
    [~, hpc_ripplefilt] = bandpower_boxcar(csc_mtx, [0 inf], 3, [150 250]);
    ripplefilts = [csc_mtx(:,1) rsc_ripplefilt hpc_ripplefilt];
    times = csc_mtx(time_idx, 1);
    
    figure; hold on
    csc = ripplefilts(time_idx, 2);
    plot(times, csc);
    rip_idx = csc_dmidx(csc_mtx, datamtx, datamtx(:,14)>0);
    times_deselect = find(rip_idx(time_idx)~=1);
    xplot_rip = times; xplot_rip(times_deselect) = nan;
    yplot_rip = csc; yplot_rip(times_deselect) = nan;
    plot(xplot_rip, yplot_rip, 'r')
    set(gca,'TickLength',[0, 0]); box off;
    
    figure; hold on  
    csc = ripplefilts(time_idx, 3);
    plot(times, csc);
    rip_idx = csc_dmidx(csc_mtx, datamtx, datamtx(:,12)>0);
    times_deselect = find(rip_idx(time_idx)~=1);
    xplot_rip = times; xplot_rip(times_deselect) = nan;
    yplot_rip = csc; yplot_rip(times_deselect) = nan;
    plot(xplot_rip, yplot_rip, 'r')
    set(gca,'TickLength',[0, 0]); box off;
end

end