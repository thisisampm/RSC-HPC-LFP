function HDposplot(datamtx)
%plot positions with HD arrows
%downsample by at least 500

%plot all spatial positions
figure; hold on; view([0 90])
plot3(datamtx(:,2), datamtx(:,3), datamtx(:,1), 'color', [1 1 1].*0.9)

%downsample
downsamp = 500;
datamtx = datamtx(datamtx(:,5)==0,:); %only vid samples
datamtx = datamtx(1:downsamp:end,:);

%plot (subsampled) HD arrows
quiver3(datamtx(:,2), datamtx(:,3), datamtx(:,1),...
    cosd(-datamtx(:,4)+90), sind(-datamtx(:,4)+90), zeros(size(datamtx(:,1)))...
    , 0.001, 'color', [0    0.4470    0.7410]); %last input determins arrow size/length

%axis labels
xlabel xpos
ylabel ypos
zlabel('time (s)')

%add black rectangle to show maze
maze_outline

%legend
legend({'position', 'heading'})

end