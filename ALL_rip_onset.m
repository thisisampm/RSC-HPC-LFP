function [all_rip_times, time_diff] = ALL_rip_onset
% combine ripple and cell events details across days
% plot

%should be input
folderpath = 'C:\Users\ampm1\Documents\MATLAB\tt_ephys\neurodata\LinTrack_prelim\';

%preallocate
%
all_rip_times = [];


%subject folders
file_list_subjects = dir(folderpath);
file_list_subjects(1:2) = [];
for isubject = 1:length(file_list_subjects)
    current_subj = file_list_subjects(isubject).name;
    
    %session folders
    file_list_sessions = dir([folderpath '\' current_subj]);
    file_list_sessions(1:2) = [];
    for isession = 1 : length(file_list_sessions)
        current_sesh = file_list_sessions(isession).name

        %load data file
        load([folderpath '\' current_subj '\' current_sesh])
        
        %compute
        [rip_nums, rip_times] = simultaneous_ripples(datamtx);

        %
        all_rip_times = [all_rip_times; [ rip_times   repmat(isession, size(rip_times(:,1))) ] ];

        
    end
end

%when does rsc ripple start relative to hpc ripple onset
time_diff = all_rip_times(:,1) - all_rip_times(:,2);

figure; plot(all_rip_times(:,3), time_diff, 'o')
figure; histogram(time_diff, [-.3:.01:.3])


