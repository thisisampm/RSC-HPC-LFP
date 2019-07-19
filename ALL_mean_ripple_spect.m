function [ALL_spectrogram_mtx, ALL_cell_rms, ALL_cell_event_times, all_ripple_power_hold] = ALL_mean_ripple_spect
% combine ripple and cell events details across days
% plot

%should be input
folderpath = 'C:\Users\ampm1\Documents\MATLAB\tt_ephys\neurodata_new';

%preallocate
%
ALL_spectrogram_mtx = [];
all_ripple_power = []; %zed
all_ripple_power_hold = []; %non zed
ALL_cell_rms = cell(2,1);
ALL_cell_event_times = [];
time_beforeafter = [-0.5 0.6];



%subject folders
file_list_subjects = dir(folderpath);
file_list_subjects(1:2) = [];
for isubject = 1:length(file_list_subjects)
    current_subj = file_list_subjects(isubject).name;
    
    %session folders
    file_list_sessions = dir([folderpath '\' current_subj]);
    file_list_sessions(1:2) = [];
    for isession = 1: length(file_list_sessions)
        current_sesh = file_list_sessions(isession).name
        
        %load data file
        load([folderpath '\' current_subj '\' current_sesh])

        %compute ripple info
        [spectrogram_mtx, cell_rms, cellevent_times, frequencies, ripple_power] = ...
            mean_ripple_spect(datamtx(ismember(datamtx(:,7), [1 5]) & datamtx(:,11)==1,:), ...
            csc_mtx, clusters, time_beforeafter, RSC_TTs, HPC_TTs, 1);

        %load ripple info
        ALL_spectrogram_mtx = cat(4, ALL_spectrogram_mtx, spectrogram_mtx);   
        ALL_cell_rms{1} = [ALL_cell_rms{1}; cell_rms{1}];
        ALL_cell_rms{2} = [ALL_cell_rms{2}; cell_rms{2}];
        all_ripple_power_hold = cat(3,all_ripple_power, ripple_power);
        all_ripple_power = cat(3,all_ripple_power, zscore_mtx(ripple_power')');
        ALL_cell_event_times = [ALL_cell_event_times cellevent_times];
        
    end
end

%averages
spectrogram_mtx = mean(ALL_spectrogram_mtx,4);
%all_ripple_power = mean(all_ripple_power,3);


%PLOT
%
%spectrogram subplots
%
figure
for i = 1:size(spectrogram_mtx,3)
    subplot(2,1,i)

    %prepare output
    s_mtx = spectrogram_mtx(:,:,i).*frequencies;
    s_mtx = 10*log10(s_mtx);

    imagesc(s_mtx)
    set(gca,'Ydir','normal')
    caxis([60 95])
    colorbar
    ylim([1 find(frequencies>300,1,'first')])

    %labels
    newyticklabels = 0:50:300;
    newyticks = interp1(frequencies, 1:size(spectrogram_mtx,1), newyticklabels);
    yticks(newyticks)
    yticklabels(newyticklabels)
    ylabel('Frequency (Hz)')
    zero_sec_bins = interp1(linspace(time_beforeafter(1), time_beforeafter(2), size(s_mtx,2)), 1:size(s_mtx,2), 0);
    xticks([1 zero_sec_bins size(s_mtx,2)]); 
    xticklabels([time_beforeafter(1) 0 time_beforeafter(2)])
    xlabel('Time (s)')

    %aesthetics
    set(gca,'TickLength',[0, 0]); box off
end

%mean activity
rm_all = [];
rm_std = [];

%cellevent subplots
figure
cell_events_out = cell(1, size(spectrogram_mtx,3));
sort_idx_hold = cell(size(cell_events_out));
for isubplot = 1:size(spectrogram_mtx,3)
    subplot(2,1,isubplot); hold on

    rm = ALL_cell_rms{isubplot};
    if ~isempty(rm)

        %load rm all
        %
        %[~,maxcols] = max(rm,[],2);
        %[~,sort_idx] = sort(maxcols);    
        %sort_idx_hold{isubplot} = flipud(sort_idx);
        rm_std = [rm_std; std(rm)];
        rm_all = [rm_all; mean(rm)];
        %rm_all = [rm_all; sum(rm)./size(rm,1)];
        %}

        %plot
        %imagesc(rm(sort_idx_hold{isubplot},:)); caxis([0 1]); colorbar
        imagesc(rm); colorbar
        set(gca,'Ydir','reverse')
        
        %label
        ylim([0.5 size(rm,1)+0.5])
        xlim([0.5 size(rm,2)+0.5])
        set(gca,'TickLength',[0, 0]); box off
        
        
        zero_sec_bins = interp1(linspace(time_beforeafter(1), time_beforeafter(2), size(rm,2)), 1:size(rm,2), 0);
        xticks([1 zero_sec_bins size(rm,2)]); 
        xticklabels([time_beforeafter(1) 0 time_beforeafter(2)])

        xlabel('Time (s)')
        ylabel('Cells')
    end
end

figure;

%mean ripple-band power
subplot(2,1,1); hold on
title('ripple-band power')
plot(mean(all_ripple_power(1,:,:),3), 'linewidth', 2, 'color', [0    0.4470    0.7410])
    plot(mean(all_ripple_power(1,:,:),3)-std(all_ripple_power(1,:,:),[],3)./sqrt(size(all_ripple_power(1,:,:),3)), 'linewidth', 1, 'color', [0    0.4470    0.7410])
    plot(mean(all_ripple_power(1,:,:),3)+std(all_ripple_power(1,:,:),[],3)./sqrt(size(all_ripple_power(1,:,:),3)), 'linewidth', 1, 'color', [0    0.4470    0.7410])
plot(mean(all_ripple_power(2,:,:),3), 'linewidth', 2, 'color', [0.8500    0.3250    0.0980])
    plot(mean(all_ripple_power(2,:,:),3)-std(all_ripple_power(2,:,:),[],3)./sqrt(size(all_ripple_power(2,:,:),3)), 'linewidth', 1, 'color', [0.8500    0.3250    0.0980])
    plot(mean(all_ripple_power(2,:,:),3)+std(all_ripple_power(2,:,:),[],3)./sqrt(size(all_ripple_power(2,:,:),3)), 'linewidth', 1, 'color', [0.8500    0.3250    0.0980])

set(gca,'TickLength',[0, 0]); box off
xlim([0.5 size(s_mtx,2)+0.5])
zero_sec_bins = interp1(linspace(time_beforeafter(1), time_beforeafter(2), size(s_mtx,2)), 1:size(s_mtx,2), 0);
    xticks([1 zero_sec_bins size(s_mtx,2)]); 
    xticklabels([time_beforeafter(1) 0 time_beforeafter(2)])
ylabel('Power (Z)')
colorbar
legend({'rsc', 'hpc'})

%mean cell activity
subplot(2,1,2); hold on
plot(rm_all(1,:), 'color', [0    0.4470    0.7410], 'linewidth', 2)
    plot(rm_all(1,:)+rm_std(1,:)./sqrt(size(ALL_cell_rms{1},1)), 'color', [0    0.4470    0.7410])
    plot(rm_all(1,:)-rm_std(1,:)./sqrt(size(ALL_cell_rms{1},1)), 'color', [0    0.4470    0.7410])
plot(rm_all(2,:), 'color', [0.8500    0.3250    0.0980], 'linewidth', 2)
    plot(rm_all(2,:)+rm_std(2,:)./sqrt(size(ALL_cell_rms{2},1)), 'color', [0.8500    0.3250    0.0980])
    plot(rm_all(2,:)-rm_std(2,:)./sqrt(size(ALL_cell_rms{2},1)), 'color', [0.8500    0.3250    0.0980])
set(gca,'TickLength',[0, 0]); box off
xlim([0.5 size(rm,2)+0.5])
zero_sec_bins = interp1(linspace(time_beforeafter(1), time_beforeafter(2), size(rm,2)), 1:size(rm,2), 0);
        xticks([1 zero_sec_bins size(rm,2)]); 
        xticklabels([time_beforeafter(1) 0 time_beforeafter(2)])
ylabel('Hz (Norm)')
colorbar

