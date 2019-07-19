function dotplot(datamtx, clusters)
%DOTPLOT plots colored dots at each location of a spike event for each
% cluster id in the row vector 'clusters'. All positions are indicated by a
% gray line.

%load held figure
figure; hold on

%legend array
legend_in{1} = 'position';

%plot all spatial positions
plot3(datamtx(:,2), datamtx(:,3), datamtx(:,1), '-', 'color', [1 1 1].*0.8)


%trials index
trl_idx = ~isnan(datamtx(:,13));

%plot cluster events
count = 0;
for iclust = clusters(clusters>0)
    count = count+1;
    
    iclust_idx = datamtx(:,8)==iclust;
    iclust_idx = iclust_idx & trl_idx; %comment out
    plot3(datamtx(iclust_idx,2), datamtx(iclust_idx,3), datamtx(iclust_idx,1), '.')
    legend_in{count+1} = num2str(iclust);
end

%set 2D view
view([0 90])

%axis
%zhold = zlim;
%axis([0 700 0 700 zhold(1)-50 zhold(2)])

%axis labels
xlabel('X position')
ylabel('Y position')
zlabel('Time (s)')

%legend
legend(legend_in)

%asthetics
set(gca, 'TickLength', [0 0])
maze_outline

end

