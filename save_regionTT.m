function save_regionTT
%add files to each nuerodata folder

%should be input
folderpath = 'C:\Users\ampm1\Documents\MATLAB\tt_ephys\neurodata\LinTrack_prelim\';



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


        RSC_TTs = 1:5;
        HPC_TTs = 7:8;
        O_TTs = 6;
        
        save([folderpath '\' current_subj '\' current_sesh], 'clusters', 'csc_mtx', 'datamtx', 'spike_waveforms', 'RSC_TTs', 'HPC_TTs', 'O_TTs')
        clearvars('clusters', 'csc_mtx', 'datamtx', 'spike_waveforms', 'RSC_TTs', 'HPC_TTs', 'O_TTs');

        
    end
end

