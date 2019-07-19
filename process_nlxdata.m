function process_nlxdata
% identify neuralynx session folders within the original project data file
% process each session folder
% Create a neurodata folder, with subject sub-folders, then fill 
% subfolders with date.mat files


% original project data file
path_origin = 'E:\Projects\InProgress\HPC-RSC-RippleMod\data\dual-RSC-HPC';

% matlab project data file
path_destination = 'C:\Users\ampm1\Documents\MATLAB\tt_ephys\';

%get all the things in origin project folder...
file_list_subjects = dir(path_origin);
file_list_subjects(1:2) = [];
file_names_subjects = {file_list_subjects([file_list_subjects(:).isdir])};

%iterate through subjects
for subject = 1:size(file_names_subjects{:},1)

    %print update
    subj_name = file_names_subjects{:}(subject,1).name
    
    %get all the things in subject folder...
    file_list_sessions = dir([path_origin '\' file_names_subjects{:}(subject,1).name]);
    file_list_sessions(1:2) = [];
    file_names_sessions = {file_list_sessions([file_list_sessions(:).isdir])};
    
    %iterate through session folders
    for session = 1:size(file_names_sessions{:},1)
        
        %current session file name
        session_file = file_names_sessions{:}(session,1).name
        
        %check if destination has a neurodata folder, if not, create it.
        if ~exist([path_destination, 'neurodata_new'], 'dir')
            mkdir([path_destination, 'neurodata_new'])
        end

        %check if destination has corresponding subject folder, if not, create it.
        if ~exist([path_destination, 'neurodata_new\', subj_name], 'dir')
            mkdir([path_destination, 'neurodata_new\', subj_name])
        end

        %session number
        if session < 10
            session_number = ['0' num2str(session)];
        else
            session_number = num2str(session);
        end

        %check if destination has corresponding session file, if so, skip
        if exist([path_destination, 'neurodata_new\', subj_name, '\', session_number, '.mat'], 'file')
            disp('skip')
            continue
        end

        %load session
        if subject==1
            %load
            [datamtx, csc_mtx, clusters, spike_waveforms, TT_regions] = loadnlx_shortmaze([path_origin, '\', subj_name, '\', session_file]);

        else
            %load
            [datamtx, csc_mtx, clusters, spike_waveforms, TT_regions] = loadnlx_longmaze([path_origin, '\', subj_name, '\', session_file]);
        end

        %tt regions
        RSC_TTs = TT_regions{1};
        HPC_TTs = TT_regions{2};
        O_TTs = TT_regions{3};

        %save session in correct destination folder
        save([path_destination, 'neurodata_new\', subj_name, '\', session_number], 'datamtx', 'csc_mtx', 'clusters', 'RSC_TTs', 'HPC_TTs', 'O_TTs', 'spike_waveforms')

        %clear variables
        clear datamtx
        clear csc_mtx
        clear clusters
        clear RSC_TTs
        clear HPC_TTs
        clear spike_waveforms
              
    end
end

end
