function all_FRs = ALL_firingrates
% combine ripple and cell events details across days
% plot

%should be input
folderpath = 'C:\Users\ampm1\Documents\MATLAB\tt_ephys\neurodata\LinTrack_prelim\';

%preallocate
%
all_FRs = [];


%subject folders
file_list_subjects = dir(folderpath);
file_list_subjects(1:2) = [];
for isubject = 1:length(file_list_subjects)
    current_subj = file_list_subjects(isubject).name;
    
    %session folders
    file_list_sessions = dir([folderpath '\' current_subj]);
    file_list_sessions(1:2) = [];
    for isession = 1 : length(file_list_sessions)
        current_sesh = file_list_sessions(isession).name;

        %load data file
        load([folderpath '\' current_subj '\' current_sesh])

        stage_idx = ismember(datamtx(:,7), [1 2 3]);
        FRs = histcounts(datamtx(stage_idx,8), [clusters(1:end, 1); clusters(end,1)+1]);
        time = sum(datamtx(stage_idx,8)==0)./100; %s
        all_FRs = [all_FRs; (FRs/time)'];

        
    end
end

