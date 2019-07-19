function [all_hpc_cwt, all_rsc_cwt, hpc_freq, rsc_freq, tr] = scalogram_ripples(datamtx, csc_mtx, ripple_column, varargin)
%plots average scalograms from all ripple events

%inputs
if nargin==4
    ploton = varargin{1};
    datamtx_base = datamtx;
elseif nargin==5
    ploton = varargin{1};
    datamtx_base = varargin{2};
else
    ploton = 1;
    datamtx_base = datamtx;
end


%hard inputs
stages = [1 5]; %1 3 and/or 5
tr = [-0.5 0.6]; %time range
%ripple_column = 12; %hpc=12, rsc=14, control=0;

%sample rate
samprate = cscmtx_samprate(csc_mtx);

%preallocate
all_hpc_cwt = [];
all_rsc_cwt = [];
all_hpc_cwt_base = [];
all_rsc_cwt_base = [];

%for all ripples
if ripple_column~=0
    rc = ripple_column;
else %find control "ripple event" times (delayed)
    delay = 1;%s
    rc = 12;
end
all_rips = unique(datamtx(ismember(datamtx(:,7), stages) & datamtx(:,rc)>0  & datamtx(:,11)==1, rc));

for irip = 1:length(all_rips)
    current_rip = all_rips(irip);
    
    %ripple times (to plot)
    lo_time = min(datamtx(datamtx(:,rc)==current_rip,1)) + tr(1);
    
    if ripple_column == 0
        lo_time = lo_time + delay;
    end
    
    lo_row = find(abs(csc_mtx(:,1) - lo_time) == min(abs(csc_mtx(:,1) - lo_time)), 1, 'first');
    hi_row = lo_row + floor((diff(tr)*samprate));
    
    %quit if we're maxed out the csc_mtx
    if hi_row > size(csc_mtx,1)
        break
    end
    
    %regional csc
    hpc_rip = csc_mtx(lo_row:hi_row, 3);
    rsc_rip = csc_mtx(lo_row:hi_row, 2);
    
    %compute continuous waveform transformations (cwt)
    [hpc_cwt, ~, ~, hpc_freq] = cwt_ampm(hpc_rip, samprate);
    [rsc_cwt, ~, ~, rsc_freq] = cwt_ampm(rsc_rip, samprate);
    
    %absolute values are customary
    hpc_cwt = abs(hpc_cwt);
    rsc_cwt = abs(rsc_cwt);
    
    %load cwts
    all_hpc_cwt = cat(3, all_hpc_cwt, hpc_cwt);
    all_rsc_cwt = cat(3, all_rsc_cwt, rsc_cwt);
end

%baseline
[all_hpc_cwt_base, all_rsc_cwt_base] = scalogram_baseline(datamtx_base, csc_mtx, tr, samprate);

%standardize by baseline
all_hpc_cwt = (mean(all_hpc_cwt,3) - mean(all_hpc_cwt_base,3)) ./ std(all_hpc_cwt_base, [], 3);
all_rsc_cwt = (mean(all_rsc_cwt,3) - mean(all_rsc_cwt_base,3)) ./ std(all_rsc_cwt_base, [], 3);

%zscore
%all_hpc_cwt = (all_hpc_cwt - mean(all_hpc_cwt(:)))./std(all_hpc_cwt(:));
%all_rsc_cwt = (all_rsc_cwt - mean(all_rsc_cwt(:)))./std(all_rsc_cwt(:));

%remove rows outside of neuralynx filter
min_over300 = min(hpc_freq(hpc_freq>300));
all_hpc_cwt = all_hpc_cwt(fliplr(hpc_freq)>=4 & fliplr(hpc_freq)<=min_over300,:);
all_rsc_cwt = all_rsc_cwt(fliplr(rsc_freq)>=4 & fliplr(rsc_freq)<=min_over300,:);
hpc_freq = hpc_freq(fliplr(hpc_freq)>=4 & fliplr(hpc_freq)<=min_over300);
rsc_freq = rsc_freq(fliplr(rsc_freq)>=4 & fliplr(rsc_freq)<=min_over300);

%plot
if ploton==1   
    figure; imagesc(all_hpc_cwt)

        %time labels
        xticks(linspace(1, size(all_hpc_cwt,2), 11))
        xticklabels(linspace(tr(1), tr(2), 11))

        %freq labels
        Yticks = 2.^(round(log2(min(hpc_freq))):round(log2(max(hpc_freq))));    
        yticks(fliplr(interp1(hpc_freq, 1:size(all_hpc_cwt,1), Yticks, 'linear', 'extrap')))
        yticklabels(fliplr(Yticks))

        %aesthetics
        title('HPC')
        set(gca,'TickLength',[0, 0]); box off; axis square;

    figure; imagesc(all_rsc_cwt)

        %time labels
        xticks(linspace(1, size(all_rsc_cwt,2), 11))
        xticklabels(linspace(tr(1), tr(2), 11))

        %freq labels
        Yticks = 2.^(round(log2(min(rsc_freq))):round(log2(max(rsc_freq))));
        yticks(fliplr(interp1(rsc_freq, 1:size(all_rsc_cwt,1), Yticks, 'linear', 'extrap')))
        yticklabels(fliplr(Yticks))

        %aesthetics
        title('RSC')
        set(gca,'TickLength',[0, 0]); box off; axis square;
end

