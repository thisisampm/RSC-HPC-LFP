function all_clusts = ALL_clusters
% combine ripple and cell events details across days
% plot

%should be input
folderpath = 'C:\Users\ampm1\Documents\MATLAB\tt_ephys\neurodata_new\';

%preallocate
%
all_clusts = [];


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


[repmat([isubject, isession], size(clusters,1),1) clusters(:,1)]

        all_clusts = [all_clusts; [repmat([isubject, isession], size(clusters,1),1) clusters(:,1)]];

        
    end
end

