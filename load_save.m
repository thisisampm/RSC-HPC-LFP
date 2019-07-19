%iterate through folders and files
load_folderpath = 'D:\HPC-RSC_replay';
save_folderpath = 'C:\Users\ampm1\Documents\MATLAB\tt_ephys\neurodata\LinTrack_prelim';
%subject folders
file_list_subjects = dir(load_folderpath);
file_list_subjects(1:2) = [];
for isubject = 1:length(file_list_subjects)
    current_subj = file_list_subjects(isubject).name;
    
    %session folders
    file_list_sessions = dir([load_folderpath '\' current_subj]);
    file_list_sessions(1:2) = [];
    
    for isession = 1:length(file_list_sessions)
        current_sesh_load = file_list_sessions(isession).name;
        current_sesh_save = file_list_sessions(isession).name(1:10)
        
        full_load_file = [load_folderpath '\' current_subj '\' current_sesh_load];
        full_save_file = [save_folderpath '\' current_subj '\' current_sesh_save];

        %load data file
        [datamtx, csc_mtx, clusters, spike_waveforms] = loadnlx(full_load_file); 

        RSC_TTs = 1:5; HPC_TTs = 7:8; O_TTs = 6;
        
        %save
        save(full_save_file, 'clusters', 'csc_mtx', 'datamtx', 'spike_waveforms', 'RSC_TTs', 'HPC_TTs', 'O_TTs');
    
        %clear
        clearvars('clusters', 'csc_mtx', 'datamtx', 'spike_waveforms', 'RSC_TTs', 'HPC_TTs', 'O_TTs');
        
    end
end
clear