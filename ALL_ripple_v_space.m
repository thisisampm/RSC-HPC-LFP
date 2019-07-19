function [rm_out, replay_scores_r, ALL_rm_space, ALL_rm_ripple, ALL_region_IDs] = ALL_ripple_v_space(sesh_stages, col, varargin)
% combine duration-normalized ripple and cell events details across days
% combine spatial cell events details across days
% plot

%inputs
if nargin > 2
    ripmod = varargin{1}(:,2) >= 1; %ripmod score and cuttoff
end

%firing rate peak or weighted center
prefered_pixle = 2; %peak=1, wc = 2;

%preallocate
%
ALL_rm_ripple = [];
ALL_rm_space = [];
ALL_clusters = [];
ALL_cluster_idx = [];
ALL_region_IDs_allclusters = [];
ALL_region_IDs = [];
ALL_spike_event_counts_ripple = [];
ALL_spike_event_counts_maze = [];


%iterate through folders and files
folderpath = 'C:\Users\ampm1\Documents\MATLAB\tt_ephys\neurodata_new';
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
        %region_ids(ismember(floor(clusters(:,1)), HPC_TTs))=2;
        region_ids(ismember(floor(clusters(:,1)), 6:8))=2;
        
        %screen for too-few trials
        if length(unique(datamtx(~isnan(datamtx(:,13)),13)))<6
            ALL_clusters = [ALL_clusters; clusters];
            ALL_cluster_idx = [ALL_cluster_idx; false(size(clusters(:,1)))];
            ALL_region_IDs_allclusters = [ALL_region_IDs_allclusters; region_ids];
            display('too few trials; skip session')
            continue
        else
            ALL_clusters = [ALL_clusters; clusters];
            ALL_cluster_idx = [ALL_cluster_idx; true(size(clusters(:,1)))];
            ALL_region_IDs_allclusters = [ALL_region_IDs_allclusters; region_ids];
            ALL_region_IDs = [ALL_region_IDs; region_ids];
        end
        
        
        
        
        %compute ripple rates
        [rm_ripple, ~, spike_event_counts_ripple] = mean_ripple_warp(datamtx(datamtx(:,11)==1 & ismember(datamtx(:,7), sesh_stages),:), clusters, col, 0);
        rm_ripple = rm_ripple - min(rm_ripple,[],2);
        rm_ripple = rm_ripple./max(rm_ripple,[],2);
        
        %compute spatial rates
        [~,rm_space, ~, spike_event_counts_maze] = rate_xdim(datamtx(datamtx(:,7)==3 & datamtx(:,11)==0,:), clusters, 0);  
        %[~,rm_space, ~, spike_event_counts_maze] = rate_xdim_dir(datamtx(datamtx(:,7)==3 & datamtx(:,11)==0,:), clusters, 0);  
        rm_space = norm_rows(rm_space);
        
        %load info
        ALL_rm_ripple = [ALL_rm_ripple; rm_ripple];
        ALL_rm_space = [ALL_rm_space; rm_space];
        ALL_spike_event_counts_ripple = [ALL_spike_event_counts_ripple; spike_event_counts_ripple];
        ALL_spike_event_counts_maze = [ALL_spike_event_counts_maze; spike_event_counts_maze];
       

    end
end
bins = size(rm_ripple,2);



%apply clusters idx
ALL_clusters_included = ALL_clusters(logical(ALL_cluster_idx));
if exist('ripmod', 'var')
    ripmod_included = ripmod(logical(ALL_cluster_idx));
else
    ripmod_included = true(sum(ALL_cluster_idx),1);
end

%compute regional rate matrices
%spike count minimum index
sc_min_ripple = 15;
sc_min_maze = 15;
sc_min_ripple_idx = ALL_spike_event_counts_ripple >= sc_min_ripple;
sc_min_maze_idx = ALL_spike_event_counts_maze >= sc_min_maze;
rm_out = cell(2,2); 
clust_out = cell(2,2);
subplot_idx = [1 2; 3 4];
figure
for space_rip = 1:2
    for rsc_hpc = 1:2
        
        if space_rip == 1
            rm = ALL_rm_space;
        elseif space_rip == 2
            rm = ALL_rm_ripple;
        end

        rm = rm(ripmod_included & ALL_region_IDs==rsc_hpc & sc_min_ripple_idx & sc_min_maze_idx,:);
        clusters_local = ALL_clusters_included(ripmod_included & ALL_region_IDs==rsc_hpc & sc_min_ripple_idx & sc_min_maze_idx,1);
        rm_out{space_rip, rsc_hpc} = rm;
        clust_out{space_rip, rsc_hpc} = clusters_local;
    end
end

%sort region rms by spatial peaks
%
for rsc_hpc  = 1:2
    [~,sort_idx] = max(rm_out{1,rsc_hpc},[],2); 
    [~,sort_idx] = sort(sort_idx);
    
    rm_out{1,rsc_hpc} = rm_out{1,rsc_hpc}(sort_idx,:);
    rm_out{2,rsc_hpc} = rm_out{2,rsc_hpc}(sort_idx,:);
    
    clust_out{1,rsc_hpc} = clust_out{1,rsc_hpc}(sort_idx,:);
    clust_out{2,rsc_hpc} = clust_out{2,rsc_hpc}(sort_idx,:);
end
%}

%compute replay score
replay_scores_r = cell(1,2);
for space_rip = 1:2
    for rsc_hpc = 1:2
        replay_scores_r{rsc_hpc} = nan(size(rm_out{1,rsc_hpc}, 1), 2);
        max_idx = false(size(replay_scores_r{rsc_hpc},1),2);
        for ic = 1:size(rm_out{1,rsc_hpc}, 1)
            replay_scores_r{rsc_hpc}(ic,:) = [corr(rm_out{1,rsc_hpc}(ic,:)', rm_out{2,rsc_hpc}(ic,:)')...
                corr(rm_out{1,rsc_hpc}(ic,:)', fliplr(rm_out{2,rsc_hpc}(ic,:))')];
            max_idx(ic, replay_scores_r{rsc_hpc}(ic, :) == max(replay_scores_r{rsc_hpc}(ic, :))) = true;
            %max_idx(ic, 1) = true;
        end
        replay_scores_r{rsc_hpc} = max(replay_scores_r{rsc_hpc},[],2);

        %SUPERVISE ORIENTAITON OF REPLAY RELATIVE TO TRACK
        %flip ripple rm if corr is better opposite way
        rm_out{2,rsc_hpc}(max_idx(:,2),:) = fliplr(rm_out{2,rsc_hpc}(max_idx(:,2),:));
    end
end



%plot rate matrices
for space_rip = 1:2
    for rsc_hpc = 1:2
        subplot(2,2,subplot_idx(space_rip, setdiff(1:2, rsc_hpc)))
        if ~isempty(rm_out{space_rip,rsc_hpc})

            %PLOT
            %
            %imagesc(norm_rows(smooth2a(rm_out{space_rip,hpc_rsc}, 1)));
            imagesc(rm_out{space_rip,rsc_hpc});

            %aesthetics
            set(gca, 'TickLength', [0 0])
            xticks(0-realmin:bins/10:bins)
            xticklabels(0:.1:1)
            yticks(1:size(rm_out{space_rip,rsc_hpc},1))
            yticklabels(clust_out{space_rip,rsc_hpc})
            ylabel('Cell IDs')

            if space_rip == 1
                xlabel('Spatial position')
            elseif space_rip == 2
                xlabel('Ripple position')
            end
        end
    end
end



function mtx_out = norm_rows(mtx_in)
    mtx_out = mtx_in - min(mtx_in,[],2);
    mtx_out = mtx_out./max(mtx_out,[],2);
end


%dot plots
%correlate peaks in space and ripple
%

%hpc
    %space
    subplot(2,2,1); hold on
    hpc_space = rm_out{1,2};
    hpc_space_max = nan(size(hpc_space,1),1);
    for i = 1:size(hpc_space,1) 
        if prefered_pixle == 1
            hpc_space_max(i) = find(hpc_space(i,:) == max(hpc_space(i,:)), 1);
        elseif prefered_pixle == 2
            hpc_space_max(i) = weighted_center(hpc_space(i,:));
        end
        plot(hpc_space_max(i), i, 'ro')        
    end
    
    %ripple
    subplot(2,2,3); hold on
    hpc_rip = rm_out{2,2};
    hpc_rip_max = nan(size(hpc_rip,1),1);
    for i = 1:size(hpc_rip,1) 

        %simple peak (max)
        if prefered_pixle == 1
            hold_max = find(hpc_rip(i,:) == max(hpc_rip(i,:)), 1);
        elseif prefered_pixle == 2
            hold_max = weighted_center(hpc_rip(i,:));
        end

        if ~isempty(hold_max)
            hpc_rip_max(i) = hold_max; 
        end
        plot(hpc_rip_max(i), i, 'ro')

    end
    
%rsc
    %space
    subplot(2,2,2); hold on
    rsc_space = rm_out{1,1};
    rsc_space_max = nan(size(rsc_space,1),1);
    for i = 1:size(rsc_space,1)
        if prefered_pixle == 1
        	rsc_space_max(i) = find(rsc_space(i,:) == max(rsc_space(i,:)), 1); 
        elseif prefered_pixle == 2
            rsc_space_max(i) = weighted_center(rsc_space(i,:));
        end
        plot(rsc_space_max(i), i, 'ro')
    end
    
    %ripple
    subplot(2,2,4); hold on
    rsc_rip = rm_out{2,1};
    rsc_rip_max = nan(size(rsc_rip,1),1);
    for i = 1:size(rsc_rip,1) 
        %simple peak (max)
        if prefered_pixle == 1
        	hold_max = find(rsc_rip(i,:) == max(rsc_rip(i,:)), 1); 
        elseif prefered_pixle == 2
            hold_max = weighted_center(rsc_rip(i,:));
        end
        
        if ~isempty(hold_max)
            rsc_rip_max(i) = hold_max; 
        end
        plot(rsc_rip_max(i), i, 'ro')
    end
    

    figure;
    [hpc_r, hpc_p] = fit_line(hpc_space_max, hpc_rip_max)
    title hpc
    xlabel space
    xticks(linspace(1,size(rsc_space,2),11))
    xticklabels(linspace(0,1,11))
    ylabel ripple
    yticks(linspace(1,size(rsc_space,2),11))
    yticklabels(linspace(0,1,11))
    set(gca,'TickLength',[0, 0]); box off
    axis([-1 size(hpc_space,2)+1 -1 size(hpc_space,2)+1])
    axis square
    
    figure;
    [rsc_r, rsc_p] = fit_line(rsc_space_max, rsc_rip_max)
    title rsc
    xlabel space
    xticks(linspace(1,size(rsc_space,2),11))
    xticklabels(linspace(0,1,11))
    ylabel ripple
    yticks(linspace(1,size(rsc_space,2),11))
    yticklabels(linspace(0,1,11))
    set(gca,'TickLength',[0, 0]); box off
    axis([-1 size(rsc_space,2)+1 -1 size(rsc_space,2)+1])
    axis square
    
    
%r score histograms
figure; histogram(replay_scores_r{1}, -1:.1:1, 'normalization', 'probability')
set(gca,'TickLength',[0, 0]); box off;
title rsc
xlim([-1 1])
hold on; plot([0 0], ylim, 'k--')
ylabel('proportion of cells')
xlabel('replay score (r)')

figure; histogram(replay_scores_r{2}, -1:.1:1, 'normalization', 'probability')
set(gca,'TickLength',[0, 0]); box off;
title hpc
xlim([-1 1])
hold on; plot([0 0], ylim, 'k--')
ylabel('proportion of cells')
xlabel('replay score (r)')


%set output to size(ALL_clusters)
ALL_cluster_idx = logical(ALL_cluster_idx);
ripmod_included = set_output(ripmod_included, ALL_cluster_idx)==1;
sc_min_ripple_idx = set_output(sc_min_ripple_idx, ALL_cluster_idx)==1;
sc_min_maze_idx = set_output(sc_min_maze_idx, ALL_cluster_idx)==1;

size(ripmod_included)
size(sc_min_ripple_idx)

rm_out{1,1} = set_output(rm_out{1,1}, ALL_cluster_idx(ripmod_included & ALL_region_IDs_allclusters==1 & sc_min_ripple_idx & sc_min_maze_idx));
rm_out{2,1} = set_output(rm_out{2,1}, ALL_cluster_idx(ripmod_included & ALL_region_IDs_allclusters==1 & sc_min_ripple_idx & sc_min_maze_idx));
rm_out{1,2} = set_output(rm_out{1,2}, ALL_cluster_idx(ripmod_included & ALL_region_IDs_allclusters==2 & sc_min_ripple_idx & sc_min_maze_idx));
rm_out{2,2} = set_output(rm_out{2,2}, ALL_cluster_idx(ripmod_included & ALL_region_IDs_allclusters==2 & sc_min_ripple_idx & sc_min_maze_idx));
replay_scores_r{1} = set_output(replay_scores_r{1}, ALL_cluster_idx(ripmod_included & ALL_region_IDs_allclusters==1 & sc_min_ripple_idx & sc_min_maze_idx));
replay_scores_r{2} = set_output(replay_scores_r{2}, ALL_cluster_idx(ripmod_included & ALL_region_IDs_allclusters==2 & sc_min_ripple_idx & sc_min_maze_idx));
ALL_rm_space = set_output(ALL_rm_space, ALL_cluster_idx);
ALL_rm_ripple = set_output(ALL_rm_ripple, ALL_cluster_idx);

%ALL_region_IDs_allclusters = set_output(ALL_region_IDs_allclusters, ALL_cluster_idx);
%if exist('ripmod', 'var')
%    ripmod = set_output(ripmod, ALL_cluster_idx);
%end

function output = set_output(input, index)
    temp_input = nan(length(index), size(input,2));    
    temp_input(index,:) = input;
    output = temp_input;
end


end