function plot_trials(datamtx)


figure; hold on
plot3(datamtx(datamtx(:,7)==3,2), datamtx(datamtx(:,7)==3,3), datamtx(datamtx(:,7)==3,1), 'color', [.9 .9 .9])
for i = 1:length(unique(datamtx(~isnan(datamtx(:,13)),13)))
    plot3(datamtx(datamtx(:,13)==i,2), datamtx(datamtx(:,13)==i,3), datamtx(datamtx(:,13)==i,1))
end
maze_outline

end