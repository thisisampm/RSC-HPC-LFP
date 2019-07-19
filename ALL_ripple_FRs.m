function [all_FRs, all_spiketimes, times] = ALL_ripple_FRs(varargin)
% combine ripple and cell events details across days
% plot

%input ripple column (12 is HPC, 14 is RSC)
if ~isempty(varargin)
    ripple_col = varargin{1};
else
    ripple_col = 12;
end

%should be input
folderpath = 'C:\Users\ampm1\Documents\MATLAB\tt_ephys\neurodata_new';

%preallocate
%
all_FRs = [];
all_spiketimes = [];


%subject folders
file_list_subjects = dir(folderpath);
file_list_subjects(1:2) = [];
for isubject = 1:2%:length(file_list_subjects)
    current_subj = file_list_subjects(isubject).name
    
    %session folders
    file_list_sessions = dir([folderpath '\' current_subj]);
    file_list_sessions(1:2) = [];
    for isession = 1 : length(file_list_sessions)
        current_sesh = file_list_sessions(isession).name

        %load data file
        load([folderpath '\' current_subj '\' current_sesh])

        %compute
        [spiketimes, num_trials] = spike_times_allrip(datamtx, clusters, [-.5 .6], ripple_col);
        [bc_rates, times] = boxcar_rates(spiketimes, [-.5 .6], .01, .001, repmat(num_trials, size(spiketimes)));
        all_FRs = [all_FRs; bc_rates];
        
        %load output
        all_spiketimes = [all_spiketimes; spiketimes];

        
    end
end

%plot
%{
figure; imagesc(times, 1:size(all_FRs,1), norm_mtx(all_FRs')')
yticks(1:size(clusters,1))
yticklabels(clusters(:,1))
set(gca,'TickLength',[0, 0]);
%}