function [all_hpc_cwt_base, all_rsc_cwt_base] = scalogram_baseline(datamtx, csc_mtx, timerng, samprate)
%finds the average scalogram for all periods in csc_mtx of size timerng
%applys filtering criteria described below

%inputs
if length(timerng)==2
    timerng = diff(timerng);
end

%filtering:
%no rip, still, bowl stages
filt_idx = datamtx(:,12)==0 & datamtx(:,11)==1 & ismember(datamtx(:,7), [1 5]);

%remaining csc_mtx
csc_mtx_filtered = csc_mtx(csc_dmidx(csc_mtx, datamtx, filt_idx),:);

%fold csc_mtx into pages size timerng
rows_in_timerng = floor(timerng*samprate)+1;
csc_mtx_filtered = csc_mtx_filtered(1: end - rem(size(csc_mtx_filtered,1), rows_in_timerng), :);
num_pages = size(csc_mtx_filtered,1)/rows_in_timerng;
csc_mtx_rshp = nan(rows_in_timerng, size(csc_mtx_filtered,2), num_pages);
for icol = 1:size(csc_mtx_filtered,2)
   csc_mtx_rshp(:,icol,:) = reshape(csc_mtx_filtered(:,icol,:), rows_in_timerng, 1, num_pages);
end

%compute all scalograms
all_hpc_cwt_base = [];
all_rsc_cwt_base = [];
for ipage = 1:size(csc_mtx_rshp,3)
    
    %regional csc
    hpc_rip = csc_mtx_rshp(:, 3, ipage);
    rsc_rip = csc_mtx_rshp(:, 2, ipage);
    
    %compute continuous waveform transformations (cwt)
    [hpc_cwt_base] = cwt_ampm(hpc_rip, samprate);
    [rsc_cwt_base] = cwt_ampm(rsc_rip, samprate);
    
    %absolute values are customary
    hpc_cwt_base = abs(hpc_cwt_base);
    rsc_cwt_base = abs(rsc_cwt_base);

    %load cwts
    all_hpc_cwt_base = cat(3, all_hpc_cwt_base, hpc_cwt_base);
    all_rsc_cwt_base = cat(3, all_rsc_cwt_base, rsc_cwt_base);
end



    