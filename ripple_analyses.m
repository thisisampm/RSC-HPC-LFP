%plot joint distribution of RSC and HPC ripple power
%[all_div, all_obs, all_exp] = ALL_rip_expected;

    %plot mean lfp and cell activity in window around ripples
    % OLD OLD OLD
    %[ALL_spectrogram_mtx, ALL_cell_rms, ALL_cell_event_times, all_ripple_power_hold] = ALL_mean_ripple_spect;

%plot scalograms (~normalized spectrograms) and rippleband power
%[ALL_scalograms_hpc_col12, ALL_scalograms_rsc_col12, hpc_freq_col12, rsc_freq_col12, timerng_col12] = ALL_scalogram_ripples(12);
%[ALL_scalograms_hpc_col14, ALL_scalograms_rsc_col14, hpc_freq_col14, rsc_freq_col14, timerng_col14] = ALL_scalogram_ripples(14);
%[ALL_scalograms_hpc_col10, ALL_scalograms_rsc_col10, hpc_freq_col10, rsc_freq_col10, timerng_col10] = ALL_scalogram_ripples(0);

%plot rsc csc temporal offset relative to hpc (is rsc slightly behind hpc?)
%check code - it output each session, i think
%[all_offsets, all_rholds] = ALL_csc_offset([150 250], 1); %bowl rest


%plot all cell event times
%region_ID_idx = ALL_regionIDs;
%[all_FRs, all_spiketimes, times] = ALL_ripple_FRs(12); %hpc
%plot_spike_times(all_FRs, all_spiketimes, times, region_ID_idx)
%[all_FRs, all_spiketimes, times] = ALL_ripple_FRs(14); %rsc
%plot_spike_times(all_FRs, all_spiketimes, times, region_ID_idx)
%[all_FRs, all_spiketimes, times] = ALL_ripple_FRs(0); %ctrl 1s-hpc-rip-del
%plot_spike_times(all_FRs, all_spiketimes, times, region_ID_idx)


%compute and plot ripple modulation scores for every cell
%{
region_ID_idx = ALL_regionIDs;
all_mod_scores_HPCrip = ALL_ripple_modulation(12); %hpc
    subplot(2,1,1);histogram(all_mod_scores_HPCrip(all_mod_scores_HPCrip(:,1)==1,2),-5:.25:5, 'normalization', 'probability'); 
    title rsc; set(gca,'TickLength',[0, 0]); box off; ylabel('proportion of cells'); xlabel('ripple modulation')
    subplot(2,1,2); histogram(all_mod_scores_HPCrip(all_mod_scores_HPCrip(:,1)==2,2),-5:.25:5, 'normalization', 'probability'); 
    title hpcRip; set(gca,'TickLength',[0, 0]); box off; ylabel('proportion of cells'); xlabel('ripple modulation')
all_mod_scores_RSCrip = ALL_ripple_modulation(14); %rsc
    subplot(2,1,1);  histogram(all_mod_scores_RSCrip(all_mod_scores_RSCrip(:,1)==1,2),-5:.25:5, 'normalization', 'probability'); 
    title rsc; set(gca,'TickLength',[0, 0]); box off; ylabel('proportion of cells'); xlabel('ripple modulation')
    subplot(2,1,2); histogram(all_mod_scores_RSCrip(all_mod_scores_RSCrip(:,1)==2,2),-5:.25:5, 'normalization', 'probability'); 
    title rscRip; set(gca,'TickLength',[0, 0]); box off; ylabel('proportion of cells'); xlabel('ripple modulation')
all_mod_scores_ctl = ALL_ripple_modulation(0); %ctrl 1s-hpc-rip-del
    subplot(2,1,1);histogram(all_mod_scores_ctl(all_mod_scores_ctl(:,1)==1,2),-5:.25:5, 'normalization', 'probability'); 
    title rsc; set(gca,'TickLength',[0, 0]); box off; ylabel('proportion of cells'); xlabel('ripple modulation')
    subplot(2,1,2); histogram(all_mod_scores_ctl(all_mod_scores_ctl(:,1)==2,2),-5:.25:5, 'normalization', 'probability'); 
    title control; set(gca,'TickLength',[0, 0]); box off; ylabel('proportion of cells'); xlabel('ripple modulation')
%}
    
    
%ripple density before and after
%{
all_ripples_per_sec = ALL_ripple_density;
for i = 1:4; cell_e{i} = all_ripples_per_sec(:,i); end
errorbar_plot(cell_e(1:2)) 
[~, b, ~, d] = ttest(cell_e{1}, cell_e{2});
title(['HPC, r=' num2str(b) ', p=' num2str(abs(d.tstat))])
ylabel('Ripples per second')
xticklabels({'Before', 'After'})
errorbar_plot(cell_e(3:4))
[~, b, ~, d] = ttest(cell_e{3}, cell_e{4});
title(['RSC, r=' num2str(b) ', p=' num2str(abs(d.tstat))])
ylabel('Ripples per second')
xticklabels({'Before', 'After'})

%}
    
%spatial characteristics


%plot spiking activity on the maze against spiking activity during ripples
sesh_stages = [1 5];
[rm_out, all_replay_scores_r, all_replay_scores_dist, ALL_rm_space, ALL_rm_ripple] = ALL_ripple_v_space(sesh_stages, all_mod_scores);

%correlate firing activity in space and ripple
%
figure; hold on;
%imagesc(smooth2a(histcounts2(rm_out{1,1}(:),rm_out{2,1}(:), 20),1));
%set(gca,'TickLength',[0, 0]); box off; axis square; set(gca,'Ydir','normal')
[hpc_r1, hpc_p1] = fit_line(rm_out{1,1}, rm_out{2,1}); 
axis square; axis([-.05 1.05 -0.05 1.05])
set(gca,'TickLength',[0, 0]); box off
title hpc
figure; hold on;
%imagesc(smooth2a(histcounts2(rm_out{1,2}(:),rm_out{2,2}(:), 20),1));
%set(gca,'TickLength',[0, 0]); box off; axis square; set(gca,'Ydir','normal')
[rsc_r1, rsc_p1] = fit_line(rm_out{1,2}, rm_out{2,2}); 
axis square; axis([-.05 1.05 -0.05 1.05])
set(gca,'TickLength',[0, 0]); box off
title rsc
%}

%plot replay fidelity by spatial coding strength
%
figure; [r_rsc, p_rsc] = fit_line(all_mis(all_region_ids==1,1), all_replay_scores_r(all_region_ids==1,1)); %rsc
xlim([-1 1])
set(gca,'TickLength',[0, 0]); box off; axis square;
%ylim([-5 5])
hold on; plot(xlim, [0 0], 'k--')
hold on; plot([0 0], ylim, 'k--')
title(['rsc  ' num2str([r_rsc, p_rsc])])

figure; [r_hpc, p_hpc] = fit_line(all_mis(all_region_ids==2,1), all_replay_scores_r(all_region_ids==2,1)); %hpc
xlim([-1 1])
set(gca,'TickLength',[0, 0]); box off; axis square;
%ylim([-5 5])
hold on; plot(xlim, [0 0], 'k--')
hold on; plot([0 0], ylim, 'k--')
title(['hpc  ' num2str([r_hpc, p_hpc])])
%}

%plot ripple modulation by spatial coding strength
%
[all_mis, all_clusts, all_region_ids] = ALL_mis;

figure; [r_rsc, p_rsc] = fit_line(all_mis(all_region_ids==1,1), all_mod_scores(all_region_ids==1,2)); %rsc
xlim([-1 1])
set(gca,'TickLength',[0, 0]); box off; axis square;
ylim([-5 5])
hold on; plot(xlim, [0 0], 'k--')
hold on; plot([0 0], ylim, 'k--')
title(['rsc  ' num2str([r_rsc, p_rsc]) ' SpaceCoding by RippleMod'])

figure; [r_hpc, p_hpc] = fit_line(all_mis(all_region_ids==2,1), all_mod_scores(all_region_ids==2,2)); %hpc
xlim([-1 1])
set(gca,'TickLength',[0, 0]); box off; axis square;
ylim([-5 5])
hold on; plot(xlim, [0 0], 'k--')
hold on; plot([0 0], ylim, 'k--')
title(['hpc  ' num2str([r_hpc, p_hpc]) ' SpaceCoding by RippleMod'])
%}

%plot replay fidelity by ripple modulation
%
figure; [r_rsc, p_rsc] = fit_line(all_mod_scores(all_region_ids==1,2), all_replay_scores_r(all_region_ids==1,1)); %rsc
ylim([-1 1])
set(gca,'TickLength',[0, 0]); box off; axis square;
%ylim([-5 5])
hold on; plot(xlim, [0 0], 'k--')
hold on; plot([0 0], ylim, 'k--')
title(['rsc  ' num2str([r_rsc, p_rsc])])

figure; [r_hpc, p_hpc] = fit_line(all_mod_scores(all_region_ids==2,2), all_replay_scores_r(all_region_ids==2,1)); %hpc
ylim([-1 1])
set(gca,'TickLength',[0, 0]); box off; axis square;
%ylim([-5 5])
hold on; plot(xlim, [0 0], 'k--')
hold on; plot([0 0], ylim, 'k--')
title(['hpc  ' num2str([r_hpc, p_hpc])])
%}


