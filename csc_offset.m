function [offset, rhold] = csc_offset(csc_mtx, bandpower_bounds)
%plots correlations between columns 2 and 3 or csc_mtx across a range of
%offsets. The more negative the peak, the more RSC is delayed compared to HPC.

%compute bandpower
[bp_rsc] = bandpower_boxcar(csc_mtx, [0 inf], 2, bandpower_bounds);
[bp_hpc] = bandpower_boxcar(csc_mtx, [0 inf], 3, bandpower_bounds);

%compute samplerate
samprate = 1600; %cscmtx_samprate(csc_mtx);

%compute correlations
%

%first sliding hpc backwards (test if rsc proceeds hpc)
rhold1 = nan(samprate/4,1); %.25s
for i = 1:length(rhold1)
    
    %shift vectors
    cbp_rsc = bp_rsc(1:end+1-i);
    cbp_hpc = bp_hpc(i:end); 
        
    %correlate and load
    nnan_idx = ~isnan(cbp_hpc) & ~isnan(cbp_rsc); 
    rhold1(i) = corr(cbp_hpc(nnan_idx), cbp_rsc(nnan_idx)); 
    
end

%then sliding rsc backwards (test if hpc proceeds rsc)
rhold2 = nan(samprate/4,1); %.25s
for i = 1:length(rhold2)
    
    %shift vectors
    cbp_rsc = bp_rsc(i:end);
    cbp_hpc = bp_hpc(1:end+1-i); 
    
    %correlate and load
    nnan_idx = ~isnan(cbp_hpc) & ~isnan(cbp_rsc); 
    rhold2(i) = corr(cbp_hpc(nnan_idx), cbp_rsc(nnan_idx)); 
    
end

%combine vectors (with first reversed)
rhold = [rhold1(end:-1:1); rhold2(2:end)]; %first item of both vectors is the same

%compute offset(offset of the rsc signal relative to hpc)
offset = (find(rhold==max(rhold))-(length(rhold)/2))/samprate;

%plot
%{
figure; hold on;
plot(rhold)
plot([(length(rhold)/2)+.5 (length(rhold)/2)+.5], ylim, 'k-')
plot([1 1].*find(rhold==max(rhold)), ylim, 'r-')
plot(xlim, [0 0], 'k--')
set(gca,'TickLength',[0, 0]); box off;
xticks(0:80:length(rhold))
xticklabels((-(length(rhold)/2):80:(length(rhold)/2))./samprate)
%}



end