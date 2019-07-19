function [all_mis, all_clusters, all_regions_ids, all_sesh] = ALL_mis
% combine ripple and cell events details across days
% plot

%should be input
folderpath = 'C:\Users\ampm1\Documents\MATLAB\tt_ephys\neurodata\LinTrack_prelim\';

%preallocate
%
all_mis = [];
all_clusters = [];
all_regions_ids = [];
all_sesh = [];


%subject folders
file_list_subjects = dir(folderpath);
file_list_subjects(1:2) = [];

ams_idx = [0 0];

for isubject = 1:length(file_list_subjects)
    current_subj = file_list_subjects(isubject).name;
    
    %session folders
    file_list_sessions = dir([folderpath '\' current_subj]);
    file_list_sessions(1:2) = [];
    for isession = 1 : length(file_list_sessions)
        current_sesh = file_list_sessions(isession).name;

        %load data file
        load([folderpath '\' current_subj '\' current_sesh])
       
        %shuffle event cell ids
        %{
        cellevents = datamtx(datamtx(:,8)>0,8);
        cellevents = cellevents(randperm(length(cellevents)));
        datamtx(datamtx(:,8)>0,8) = cellevents;
        %}

        if length(unique(datamtx(datamtx(:,13)>0 & ~isnan(datamtx(:,13)),13))) < 4
            mis = nan;
        else
            [mis] = miller_info_score(datamtx, clusters, 20, 0);
        end
        all_mis = [all_mis; mis];
        all_clusters = [all_clusters; clusters(:,1)];
        
        %code clusters in terms of regions
        region_ids = zeros(size(clusters(:,1)));
        %RSCtts = 1
        region_ids(ismember(floor(clusters(:,1)), RSC_TTs))=1;
        %HPCtts = 2
        region_ids(ismember(floor(clusters(:,1)), HPC_TTs))=2;
        all_regions_ids = [all_regions_ids; region_ids];
        
        all_sesh = [all_sesh; [repmat(isubject, size(clusters(:,1))) repmat(isession, size(clusters(:,1)))]];
    
        ams_idx = [ams_idx(2) ams_idx(2)];
    
    end
end


%plot


