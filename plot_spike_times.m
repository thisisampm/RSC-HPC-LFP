function plot_spike_times(all_FRs, all_spiketimes, times, region_ID_idx)
%plotting function for ALL_ripple_FRs

%downsamp for each region, rsc and hpc
reg_downsamps = [15 5];

%iterate plotting rasters and heatmaps
for regionid = 1:2
    if sum(region_ID_idx==regionid)==0
        continue
    end
    all_spikes_local = all_spiketimes(region_ID_idx==regionid);
    all_FRs_local = all_FRs(region_ID_idx==regionid, :);
    spike_raster(all_spikes_local, [-.5 .6], reg_downsamps(regionid)); 
    title(['Region ' num2str(regionid)]) 
    figure; imagesc(times, 1:size(all_FRs_local,1), norm_mtx(all_FRs_local')'); 
    title(['Region ' num2str(regionid)]) 
    set(gca,'TickLength',[0, 0]); box off;
    yticks(1:sum(region_ID_idx==regionid))
    %yticklabels(clusters(region_ID_idx==regionid,1))
    set(gca,'TickLength',[0, 0]);
end


%line plot of mean activity +/-SE
rsc_rm_mean = mean(zscore_mtx(all_FRs(region_ID_idx==1,:)'), 2);
rsc_rm_se = std(zscore_mtx(all_FRs(region_ID_idx==1,:)'), [], 2)./sqrt(sum(region_ID_idx==1));
hpc_rm_mean = mean(zscore_mtx(all_FRs(region_ID_idx==2,:)'), 2);
hpc_rm_se = std(zscore_mtx(all_FRs(region_ID_idx==2,:)'), [], 2)./sqrt(sum(region_ID_idx==2));
figure; hold on;
plot(hpc_rm_mean, 'color', [1 1 1].*0.4, 'linewidth', 3)
plot(hpc_rm_mean+hpc_rm_se, 'color', [1 1 1].*0.4)
plot(hpc_rm_mean-hpc_rm_se, 'color', [1 1 1].*0.4)
plot(rsc_rm_mean, 'color', [1 1 1].*0.8, 'linewidth', 3)
plot(rsc_rm_mean+rsc_rm_se, 'color', [1 1 1].*0.8)
plot(rsc_rm_mean-rsc_rm_se, 'color', [1 1 1].*0.8)

xlim([1 length(hpc_rm_mean)])
xticks(1:length(hpc_rm_mean)/11:length(hpc_rm_mean))
xticklabels(times(1:length(hpc_rm_mean)/11:length(hpc_rm_mean)))
hold on; plot(xlim, [0 0], 'k--')
set(gca,'TickLength',[0, 0]); box off;
xlabel('Time (s)')
ylabel('Firing Rate (z)')
