function all_mod_scores = ALL_ripple_modulation(col)
% combine ripple and cell events details across days
% plot


%region of interest column
%col % HPC=12, RSC=14, control= 0

display('ALL_ripple_modulation')

%should be input
folderpath = 'C:\Users\ampm1\Documents\MATLAB\tt_ephys\neurodata_new';

%preallocate
%
all_mod_scores = [];

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
        load([folderpath '\' current_subj '\' current_sesh])
        
        %ripple modulation
        %
        %original
        mod_scores = ripple_modulation(datamtx(ismember(datamtx(:,7), [5]) & datamtx(:,11)==1,:), csc_mtx, clusters, RSC_TTs, HPC_TTs, col);
        %
        %buszaki
        %mod_scores = 
        
        %load 
        all_mod_scores = [all_mod_scores; mod_scores];
        
    end
end


%plot


