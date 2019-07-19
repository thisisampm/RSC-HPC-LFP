function spike_raster(spiketimes, timerange, varargin)
%plots vertical tick lines at each time in each cell of spike_times alone
%a x axis defined by timerange

%input subsampling
if ~isempty(varargin)
    subsamp_div = varargin{1};
else
    subsamp_div = 1;
end


figure; hold on; 
for i = 1:length(spiketimes) 
    times = spiketimes{i}; 
    sub_times = times(1:subsamp_div:end);
    plot([sub_times sub_times]', [repmat(i-.5, size(sub_times)) repmat(i+.5, size(sub_times))]', 'k'); 
end 
xlim(timerange)
ylim([0.5 length(spiketimes)+0.5])
set(gca,'TickLength',[0, 0]); box off;
set(gca,'Ydir','reverse')

end