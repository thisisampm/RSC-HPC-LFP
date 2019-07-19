function cellevent_times_out = cellevent_plot(datamtx, clusters, timerange, varargin)
%calculates spike times within a timerange. can plot

colororder = get(gca,'ColorOrder');

if nargin == 4
    ploton = varargin{1};
else
    ploton = 1;
end

if ploton == 1
    %figure; hold on
end

cellevent_times_out = cell(1, length(clusters));
for i = 1:size(clusters,1)
    cellevent_times = datamtx(datamtx(:,1)>=timerange(1) & datamtx(:,1)<=timerange(2) & datamtx(:,8) == clusters(i), 1);
    cellevent_times_out{i} = cellevent_times;
    
    %downsample
    if ploton==1
        %
        if length(cellevent_times)>100
            %cellevent_times = cellevent_times(1:floor(length(cellevent_times)/20):length(cellevent_times));
        end
        %}
        %cellevent_times
        plot(([1 1].*cellevent_times)', repmat([i-.5 i+.5], size(cellevent_times,1), 1)', '-', 'color', colororder(i,:))
        
    end
end

if ploton == 1
    ylim([0.25 size(clusters,1)+0.75])
    xlim([timerange(1) timerange(2)])
    set(gca,'TickLength',[0, 0]); box off
    yticks(1:size(clusters,1))
    yticklabels(clusters)
    xlabel('Time (s)')
    ylabel('Cell IDs')
end
end