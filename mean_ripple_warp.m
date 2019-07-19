function [rm, cellevent_times, spike_event_counts] = mean_ripple_warp(datamtx, clusters, col, varargin)
%plot the mean spectrogram for all ripples (aligned to ripple start)

%col is 12 for hpc ripple, 14 for rsc ripple

%input check
if nargin==4
   ploton = varargin{1};
else
   ploton = 1;
end

%speed prep
num_cells = size(clusters,1);



%calculate times between ripples
a = datamtx(datamtx(:,8)==0,12);
a(a==0)=nan;
b = diff(find([1,diff(a)',1]));
b = b((b~=1));



%calculate spike firing times within each ripple event
rip_count = 0;
cellevent_times = cell(length(clusters),1);
for irip = min(datamtx(datamtx(:,col)>0,col)):max(datamtx(:,col))
    rip_count = rip_count+1;
    
    %ripple time bounds
    %timerng = [min(datamtx(datamtx(:,col)==irip,1))-0.05 max(datamtx(datamtx(:,col)==irip,1))+0.05] ;
    timerng = [min(datamtx(datamtx(:,col)==irip,1)) max(datamtx(datamtx(:,col)==irip,1))] ;
        
    if ~isempty(timerng)
        ripple_duration = diff(timerng);    
        
        %cell events
        [~, normalized_spiketimes] = spike_times(datamtx, clusters, timerng, 0);
        
        %load
        for ic = 1:num_cells
            cellevent_times{ic} = sort([cellevent_times{ic}; normalized_spiketimes{ic}]);
        end
    end
end

%number of spikes for each cell
spike_event_counts = nan(num_cells,1);
for ic = 1:num_cells
    spike_event_counts(ic) = length(cellevent_times{ic});
end

%calculate normalized rates for each cell
rm = boxcar_rates(cellevent_times, [0 1], 0.1, 0.01, repmat(rip_count, size(cellevent_times)));
rm = norm_mtx(rm')';


%plot
if ploton==1 && ~isempty(rm)

    figure
    imagesc(rm)

    %label
    ylim([0.5 size(rm,1)+0.5])
    xlim([0.5 size(rm,2)+0.5])
    set(gca,'TickLength',[0, 0]); box off
    xticks(linspace(1,size(rm,2),11))
    xticklabels(linspace(0,1,11))
    yticks(1:length(clusters))
    yticklabels(clusters)
    xlabel('Time during ripple (norm)')
    ylabel('Cell IDs')
    colorbar; caxis([0 1])

end        



    
    
end

    

