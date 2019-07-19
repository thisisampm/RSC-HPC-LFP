function [rm_out] = ALL_rate_xdim
% combine duration-normalized ripple and cell events details across days
% combine spatial cell events details across days
% plot

%preallocate
%
ALL_rm_space = [];
ALL_clusters_ = [];
ALL_cluster_idx = [];
ALL_region_IDs_allclusters = [];
ALL_region_IDs = [];
ALL_spike_event_counts_maze = [];

%iterate through folders and files
folderpath = 'C:\Users\ampm1\Documents\MATLAB\tt_ephys\neurodata\LinTrack_prelim\';
%subject folders
file_list_subjects = dir(folderpath);
file_list_subjects(1:2) = [];
ripmod_index = [0 0];
for isubject = 1:length(file_list_subjects)
    current_subj = file_list_subjects(isubject).name;
    
    %session folders
    file_list_sessions = dir([folderpath '\' current_subj]);
    file_list_sessions(1:2) = [];
    for isession = 1 : length(file_list_sessions)
        current_sesh = file_list_sessions(isession).name

        %load data file
        clusters = [];
        datamtx = [];
        csc_mtx = [];
        spike_waveforms = [];
        RSC_TTs = [];
        HPC_TTs = [];
        O_TTs = [];
        load([folderpath '\' current_subj '\' current_sesh])
        
        
        %code clusters in terms of regions
        region_ids = zeros(size(clusters(:,1)));
        %RSCtts = 1
        region_ids(ismember(floor(clusters(:,1)), RSC_TTs))=1;
        %HPCtts = 2
        region_ids(ismember(floor(clusters(:,1)), HPC_TTs))=2;
        
        %screen for too-few trials
        if length(unique(datamtx(~isnan(datamtx(:,13)),13)))<6
            ALL_clusters_ = [ALL_clusters_; clusters];
            ALL_cluster_idx = [ALL_cluster_idx; false(size(clusters(:,1)))];
            ALL_region_IDs_allclusters = [ALL_region_IDs_allclusters; region_ids];
            display('too few trials; skip session')
            continue
        else
            ALL_clusters_ = [ALL_clusters_; clusters];
            ALL_cluster_idx = [ALL_cluster_idx; true(size(clusters(:,1)))];
            ALL_region_IDs_allclusters = [ALL_region_IDs_allclusters; region_ids];
            ALL_region_IDs = [ALL_region_IDs; region_ids];
        end
        

       
        %compute spatial rates
        [~,rm_space, ~, spike_event_counts_maze] = rate_xdim(datamtx(datamtx(:,7)==3 & datamtx(:,11)==0,:), clusters, 0);  
        %[~,rm_space, ~, spike_event_counts_maze] = rate_xdim_dir(datamtx(datamtx(:,7)==3 & datamtx(:,11)==0,:), clusters, 0);  
        rm_space = norm_rows(rm_space);
        
        %load info
        ALL_rm_space = [ALL_rm_space; rm_space];
        ALL_spike_event_counts_maze = [ALL_spike_event_counts_maze; spike_event_counts_maze];


    end
end
bins = size(rm_space,2);


%apply clusters idx
ALL_clusters_included = ALL_clusters_(logical(ALL_cluster_idx));


%compute regional rate matrices
%spike count minimum index
sc_min_maze = 15;
sc_min_maze_idx = ALL_spike_event_counts_maze >= sc_min_maze;
rm_out = cell(1,2); 
clust_out = cell(1,2);
figure
for rsc_hpc = 1:2
    rm = ALL_rm_space;
    rm = rm(ALL_region_IDs==rsc_hpc & sc_min_maze_idx,:);
    clusters_local = ALL_clusters_included(ALL_region_IDs==rsc_hpc & sc_min_maze_idx,1);
    rm_out{rsc_hpc} = rm;
    clust_out{rsc_hpc} = clusters_local;
end

%sort region rms by spatial peaks
%
for rsc_hpc  = 1:2
    [~,sort_idx] = max(rm_out{rsc_hpc},[],2); 
    [~,sort_idx] = sort(sort_idx);
    rm_out{rsc_hpc} = rm_out{rsc_hpc}(sort_idx,:);
    clust_out{rsc_hpc} = clust_out{rsc_hpc}(sort_idx,:);
end


%plot rate matrices
for rsc_hpc = 1:2
    subplot(2,1,rsc_hpc)
    if ~isempty(rm_out{rsc_hpc})

        %PLOT
        %
        imagesc(rm_out{rsc_hpc});

        %aesthetics
        set(gca, 'TickLength', [0 0])
        xticks(0-realmin:bins/10:bins)
        xticklabels(0:.1:1)
        yticks(1:size(rm_out{rsc_hpc},1))
        yticklabels(clust_out{rsc_hpc})
        ylabel('Cell IDs')

    end
end

function mtx_out = norm_rows(mtx_in)
    mtx_out = mtx_in - min(mtx_in,[],2);
    mtx_out = mtx_out./max(mtx_out,[],2);
end



end