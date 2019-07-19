function [all_hpc_cwt, all_rsc_cwt, hpc_freq, rsc_freq, tr] = scalogram_timeseries(datamtx, csc_mtx, datamtx_base, varargin)
%plots average scalograms from all ripple events

%inputs
if nargin==4
    ploton = varargin{1};
else
    ploton = 1;
end

%sample rate
samprate = cscmtx_samprate(csc_mtx);

%find csc_mtx rows corresponding to datamtx
lo_time = datamtx(1,1);
hi_time = datamtx(end,1);
lo_row = find(abs(csc_mtx(:,1) - lo_time) == min(abs(csc_mtx(:,1) - lo_time)), 1, 'first');
hi_row = find(abs(csc_mtx(:,1) - hi_time) == min(abs(csc_mtx(:,1) - hi_time)), 1, 'last');

%regional csc
hpc_rip = csc_mtx(lo_row:hi_row, 3);
rsc_rip = csc_mtx(lo_row:hi_row, 2);

%compute continuous waveform transformations (cwt)
[hpc_cwt, ~, ~, hpc_freq] = cwt_ampm(hpc_rip, samprate);
[rsc_cwt, ~, ~, rsc_freq] = cwt_ampm(rsc_rip, samprate);

%absolute values are customary
hpc_cwt = abs(hpc_cwt);
rsc_cwt = abs(rsc_cwt);

%baseline
[all_hpc_cwt_base, all_rsc_cwt_base] = scalogram_baseline(datamtx_base, csc_mtx, hi_time-lo_time, samprate);

%standardize by baseline
hpc_cwt = hpc_cwt ./ std(all_hpc_cwt_base, [], 3);
rsc_cwt = rsc_cwt ./ std(all_rsc_cwt_base, [], 3);


%remove rows outside of neuralynx filter
min_over300 = min(hpc_freq(hpc_freq>300));
hpc_cwt = hpc_cwt(fliplr(hpc_freq)>=4 & fliplr(hpc_freq)<=min_over300,:);
rsc_cwt = rsc_cwt(fliplr(rsc_freq)>=4 & fliplr(rsc_freq)<=min_over300,:);
hpc_freq = hpc_freq(fliplr(hpc_freq)>=4 & fliplr(hpc_freq)<=min_over300);
rsc_freq = rsc_freq(fliplr(rsc_freq)>=4 & fliplr(rsc_freq)<=min_over300);

%plot
if ploton==1   
    figure; imagesc(hpc_cwt)

        %time labels
        xticks(linspace(1, size(hpc_cwt,2), 11))
        xticklabels(linspace(lo_time, hi_time, 11))

        %freq labels
        Yticks = 2.^(round(log2(min(hpc_freq))):round(log2(max(hpc_freq))));    
        yticks(fliplr(interp1(hpc_freq, 1:size(hpc_cwt,1), Yticks, 'linear', 'extrap')))
        yticklabels(fliplr(Yticks))

        %aesthetics
        title('HPC')
        set(gca,'TickLength',[0, 0]); box off; axis square;

    figure; imagesc(rsc_cwt)

        %time labels
        xticks(linspace(1, size(rsc_cwt,2), 11))
        xticklabels(linspace(lo_time, hi_time, 11))

        %freq labels
        Yticks = 2.^(round(log2(min(rsc_freq))):round(log2(max(rsc_freq))));
        yticks(fliplr(interp1(rsc_freq, 1:size(rsc_cwt,1), Yticks, 'linear', 'extrap')))
        yticklabels(fliplr(Yticks))

        %aesthetics
        title('RSC')
        set(gca,'TickLength',[0, 0]); box off; axis square;
end

