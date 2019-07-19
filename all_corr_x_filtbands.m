function [all_bincors_1, all_bincors_2] = all_corr_x_filtbands

%should be input
folderpath = 'C:\Users\ampm1\Documents\MATLAB\tt_ephys\neurodata_new\';

%preallocate
%
all_bincors_1 = [];
all_bincors_2 = [];

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
        csc_mtx = [];
        datamtx = [];
        load([folderpath '\' current_subj '\' current_sesh])
        
        %compute
        %
        %{
        %time check BOWLS 1&2
        rest_time = length(datamtx(datamtx(:,11)==1 & ismember(datamtx(:,7), [1 5]) & datamtx(:,8)==0,1))/100;
        if rest_time < 60 %s
            bincors_1 = nan;
        else
            bincors_1 = cor_x_filtbands(csc_mtx(csc_dmidx(csc_mtx, datamtx, datamtx(:,11)==1 & ismember(datamtx(:,7), [1 5])),:));
        end
        
        %time check MAZE RUNS
        
        run_time = length(datamtx(datamtx(:,11)==0 & datamtx(:,13)>0 & datamtx(:,8)==0,1))/100;
        if run_time < 60 %s
            bincors_2 = nan;
        else
            bincors_2 = cor_x_filtbands(csc_mtx(csc_dmidx(csc_mtx, datamtx, datamtx(:,11)==0 & datamtx(:,13)>0),:));
        end
        %}
        
        
        %time check BOWL 1
        rest_time = length(datamtx(datamtx(:,11)==1 & ismember(datamtx(:,7), [1]) & datamtx(:,8)==0,1))/100;
        if rest_time < 60 %s
            bincors_1 = nan;
        else
            bincors_1 = cor_x_filtbands(csc_mtx(csc_dmidx(csc_mtx, datamtx, datamtx(:,11)==1 & ismember(datamtx(:,7), [1])),:));
        end
        
        
        %time check BOWL 2
        %
        rest_time2 = length(datamtx(datamtx(:,11)==1 & ismember(datamtx(:,7), [5]) & datamtx(:,8)==0,1))/100;
        if rest_time2 < 60 %s
            bincors_2 = nan;
        else
            bincors_2 = cor_x_filtbands(csc_mtx(csc_dmidx(csc_mtx, datamtx, datamtx(:,11)==1 & ismember(datamtx(:,7), [5])),:));
        end
        %}
        
        %load
        all_bincors_1 = [all_bincors_1 bincors_1(:,1)];
        all_bincors_2 = [all_bincors_2 bincors_2(:,1)];
        
    end
end




end