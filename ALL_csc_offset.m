function [all_offsets, all_rholds] = ALL_csc_offset(bandpower_bounds, beh_idx)
%should be input
folderpath = 'C:\Users\ampm1\Documents\MATLAB\tt_ephys\neurodata_new';

%preallocate
%
all_offsets = [];
all_rholds = [];

 
%subject folders
file_list_subjects = dir(folderpath);
file_list_subjects(1:2) = [];
for isubject = 1:length(file_list_subjects)
    current_subj = file_list_subjects(isubject).name
    
    %session folders
    file_list_sessions = dir([folderpath '\' current_subj]);
    file_list_sessions(1:2) = [];
    for isession = 1 : length(file_list_sessions)
        current_sesh = file_list_sessions(isession).name

        %load data file
        datamtx = [];
        csc_mtx = [];
        load([folderpath '\' current_subj '\' current_sesh])

        
        %set behavior index
        if beh_idx == 1 %bowl rest
            dmidx = datamtx(:,11)==1 & ismember(datamtx(:,7), [1 5]);
        elseif beh_idx == 2 %maze run
            dmidx = datamtx(:,13)>0 & datamtx(:,11)==0;
        end
        
        [offset, rhold] = csc_offset(csc_mtx(csc_dmidx(csc_mtx, datamtx, dmidx), :), bandpower_bounds);
        
        %bug bandaid
        if length(offset)>1
            offset = offset(offset~=0);
        end
        
        %load
        all_offsets = [all_offsets; offset];
        all_rholds = [all_rholds rhold];
        
    end
    
end

%plot
figure; hold on;
mean_rhold = mean(all_rholds,2);
plot(mean_rhold)
plot([(length(mean_rhold)/2)+.5 (length(mean_rhold)/2)+.5], ylim, 'k-')
plot([1 1].*find(mean_rhold==max(mean_rhold)), ylim, 'r-')
plot(xlim, [0 0], 'k--')
set(gca,'TickLength',[0, 0]); box off;
xticks(0:80:length(mean_rhold))
xticklabels((-(length(mean_rhold)/2):80:(length(mean_rhold)/2))./1600)
