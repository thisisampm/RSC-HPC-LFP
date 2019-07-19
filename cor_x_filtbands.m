function bin_cors = cor_x_filtbands(csc_mtx)
%calculates the correlation between RSC and HPC ba[150 ndpower at bins within
%the range of bandbounds

%preallocate
bin_edges = [4 10 20 40 80 150 250 400];
bins = length(bin_edges)-1;
bin_cors = nan(bins,2);


%parse band bounds into bins
%bin_edges = linspace(band_bounds(1), band_bounds(2), bins+1);

rsc_bandpower = bandpower_boxcar(csc_mtx, [0 inf], 2, [4 400]);
hpc_bandpower = bandpower_boxcar(csc_mtx, [0 inf], 3, [4 400]);
nnan_idx = ~isnan(rsc_bandpower) & ~isnan(hpc_bandpower);
bin_cors = corr(rsc_bandpower(nnan_idx), hpc_bandpower(nnan_idx));

%iterate through each bin
%{
for ibin = 1:bins

    %calculate bandpower
    rsc_bandpower = bandpower_boxcar(csc_mtx, [0 inf], 2, [bin_edges(ibin) bin_edges(ibin+1)]);
    hpc_bandpower = bandpower_boxcar(csc_mtx, [0 inf], 3, [bin_edges(ibin) bin_edges(ibin+1)]);

    %correlate
    nnan_idx = ~isnan(rsc_bandpower) & ~isnan(hpc_bandpower);
    [r, p] = corr(rsc_bandpower(nnan_idx), hpc_bandpower(nnan_idx));
    
    %load output
    bin_cors(ibin, :) = [r, p];
       
end




%plot barplot of correlations
figure; hold on
bar(bin_cors(:,1));
ylim([0 1])
set(gca,'TickLength',[0, 0]); box off;

b = (bin_edges(1:end-1) + bin_edges(2:end))./2;
xticks(1:bins)
xticklabels(b)
%}
end