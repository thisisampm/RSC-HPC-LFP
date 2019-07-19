function [ALL_scalograms_hpc, ALL_scalograms_rsc, hpc_freq, rsc_freq, timerng] = ALL_scalogram_ripples(ripple_column)
% combine ripple and cell events details across days
% plot
%ripple column --->hpc=12, rsc=14, random=0

%should be input
folderpath = 'C:\Users\ampm1\Documents\MATLAB\tt_ephys\neurodata_new';




%preallocate
%
ALL_scalograms_hpc = [];
ALL_scalograms_rsc = [];


%subject folders
file_list_subjects = dir(folderpath);
file_list_subjects(1:2) = [];
for isubject = 1:length(file_list_subjects)
    current_subj = file_list_subjects(isubject).name;
    
    %session folders
    file_list_sessions = dir([folderpath '\' current_subj]);
    file_list_sessions(1:2) = [];
    for isession = 1:length(file_list_sessions)
        current_sesh = file_list_sessions(isession).name

        %load data file
        load([folderpath '\' current_subj '\' current_sesh])

        %compute
        [all_hpc_cwt, all_rsc_cwt, hpc_freq, rsc_freq, timerng] = scalogram_ripples(datamtx, csc_mtx, ripple_column, 0);

        %load ripple info
        ALL_scalograms_hpc = cat(3, ALL_scalograms_hpc, all_hpc_cwt);   
        ALL_scalograms_rsc = cat(3, ALL_scalograms_rsc, all_rsc_cwt);
        
    end
end


%averages
all_hpc_cwt = mean(ALL_scalograms_hpc,3);
all_rsc_cwt = mean(ALL_scalograms_rsc,3);

%zscore
%all_hpc_cwt = (all_hpc_cwt-mean(all_hpc_cwt(:)))./std(all_hpc_cwt(:));
%all_rsc_cwt = (all_rsc_cwt-mean(all_rsc_cwt(:)))./std(all_rsc_cwt(:));


%PLOT
figure; imagesc(all_hpc_cwt)

    %time labels
    xticks(linspace(1, size(all_hpc_cwt,2), 11))
    xticklabels(linspace(timerng(1), timerng(2), 11))

    %freq labels
    yticks_freq = 2.^(round(log2(min(hpc_freq))):round(log2(max(hpc_freq))));
    new_yticks = fliplr(interp1(hpc_freq, 1:size(all_hpc_cwt,1), yticks_freq,'spline','extrap'));
    yticks(new_yticks)
    yticklabels(fliplr(yticks_freq))

    %aesthetics
    title('HPC')
    set(gca,'TickLength',[0, 0]); box off;
    ylabel('Frequency (Hz)')
    xlabel('Time (s)')
    colorbar; caxis_hold = caxis;

figure; imagesc(all_rsc_cwt)

    %time labels
    xticks(linspace(1, size(all_rsc_cwt,2), 11))
    xticklabels(linspace(timerng(1), timerng(2), 11))

    %freq labels
    yticks_freq = 2.^(round(log2(min(rsc_freq))):round(log2(max(rsc_freq))));
    new_yticks = fliplr(interp1(rsc_freq, 1:size(all_rsc_cwt,1), yticks_freq,'spline','extrap'));
    yticks(new_yticks)
    yticklabels(fliplr(yticks_freq))

    %aesthetics
    title('RSC')
    set(gca,'TickLength',[0, 0]); box off;
    ylabel('Frequency (Hz)')
    xlabel('Time (s)')
    colorbar; caxis(caxis_hold)
    
    
figure; hold on
set(gca,'TickLength',[0, 0]); box off;
ripple_rng_idx = hpc_freq>=100 & hpc_freq<=250;

plot(mean(mean(ALL_scalograms_hpc(ripple_rng_idx,:,:),3)))
    plot(mean(mean(ALL_scalograms_hpc(ripple_rng_idx,:,:),3))+ mean(std(ALL_scalograms_hpc(ripple_rng_idx,:,:),[],3))./sqrt(size(ALL_scalograms_hpc,3)))
    plot(mean(mean(ALL_scalograms_hpc(ripple_rng_idx,:,:),3))- mean(std(ALL_scalograms_hpc(ripple_rng_idx,:,:),[],3))./sqrt(size(ALL_scalograms_hpc,3)))

plot(mean(mean(ALL_scalograms_rsc(ripple_rng_idx,:,:),3)))
    plot(mean(mean(ALL_scalograms_rsc(ripple_rng_idx,:,:),3))+ mean(std(ALL_scalograms_rsc(ripple_rng_idx,:,:),[],3))./sqrt(size(ALL_scalograms_rsc,3)))
    plot(mean(mean(ALL_scalograms_rsc(ripple_rng_idx,:,:),3))- mean(std(ALL_scalograms_rsc(ripple_rng_idx,:,:),[],3))./sqrt(size(ALL_scalograms_rsc,3)))

        
        