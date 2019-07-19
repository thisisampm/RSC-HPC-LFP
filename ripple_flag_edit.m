function ripple_flag_edit
%loads and resaves session files with new ripple flags

%should be input
folderpath = 'C:\Users\ampm1\Documents\MATLAB\tt_ephys\neurodata\LinTrack_prelim\';

%region ripple power column (9=rsc, 10=hpc)
col_oi = 9;

%save column (14 = rsc, 12 = hpc);
save_col = 14;


%subject folders
file_list_subjects = dir(folderpath);
file_list_subjects(1:2) = [];
for isubject = 1:length(file_list_subjects)
    current_subj = file_list_subjects(isubject).name;
    
    %session folders
    file_list_sessions = dir([folderpath '\' current_subj]);
    file_list_sessions(1:2) = [];
    for isession = 1 : length(file_list_sessions)
        current_sesh = file_list_sessions(isession).name

        %load data file
        load([folderpath '\' current_subj '\' current_sesh])

            %Flag HPC ripple windows
            %
            power_threshold = 3; %z
            windowsize = 0.05; %s

            %stillness criteria
            still_mean = nanmean(datamtx(datamtx(:,11)==1 & ismember(datamtx(:,7), [1 5]), col_oi));
            still_std = nanstd(datamtx(datamtx(:,11)==1 & ismember(datamtx(:,7), [1 5]), col_oi));

            %preallocate
            ripple_flags = datamtx(:,col_oi)>(still_mean + power_threshold*still_std);
            ripple_flags_hold = zeros(size(datamtx(:,col_oi)));

            %flag ripple windows
            for ivs = 1:length(datamtx(:,col_oi))
                if ripple_flags(ivs)==1
                    %current timewindow idx
                    tw_idx = datamtx(:,1)>datamtx(ivs,1)-windowsize/2 & datamtx(:,1)<datamtx(ivs,1)+windowsize/2;
                    ripple_flags_hold(tw_idx) = 1;
                end
            end

            %count ripple windows
            count = 2; 
            for ivs = 1:length(ripple_flags_hold)
                if ripple_flags_hold(ivs)==1
                    reach = 0;
                    while ripple_flags_hold(ivs+reach)==1
                        ripple_flags_hold(ivs+reach) = count;
                        reach = reach+1;
                    end
                    count = count+1;
                end
            end
            ripple_flags_hold(ripple_flags_hold>0) = ripple_flags_hold(ripple_flags_hold>0)-1;
            

            %merge in datamtx
            datamtx(:,save_col) = ripple_flags_hold;

        save([folderpath '\' current_subj '\' current_sesh], 'clusters', 'csc_mtx', 'datamtx', 'spike_waveforms', 'RSC_TTs', 'HPC_TTs', 'O_TTs')
        clearvars('clusters', 'csc_mtx', 'datamtx', 'spike_waveforms', 'RSC_TTs', 'HPC_TTs', 'O_TTs');

    end
end

