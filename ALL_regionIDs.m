function region_ID_idx = ALL_regionIDs
% combine ripple and cell events details across days
% plot

%should be input
folderpath = 'C:\Users\ampm1\Documents\MATLAB\tt_ephys\neurodata_new';

%preallocate
%
region_ID_idx = [];


%subject folders
file_list_subjects = dir(folderpath);
file_list_subjects(1:2) = [];
for isubject = 1:2%length(file_list_subjects)
    current_subj = file_list_subjects(isubject).name;
    
    %session folders
    file_list_sessions = dir([folderpath '\' current_subj]);
    file_list_sessions(1:2) = [];
    for isession = 1 : length(file_list_sessions)
        current_sesh = file_list_sessions(isession).name;

        %load data file
        load([folderpath '\' current_subj '\' current_sesh])


        region_ids = zeros(size(clusters(:,1)));
        %RSCtts = 1
        region_ids(ismember(floor(clusters(:,1)), RSC_TTs))=1;
        %HPCtts = 2
        region_ids(ismember(floor(clusters(:,1)), HPC_TTs))=2;

        region_ID_idx = [region_ID_idx; region_ids];
    end
end

