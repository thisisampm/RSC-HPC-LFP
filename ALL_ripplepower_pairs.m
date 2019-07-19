function ALL_ripplepower_pairs
% combine ripple and cell events details across days
% plot

%should be input
folderpath = 'C:\Users\ampm1\Documents\MATLAB\tt_ephys\neurodata\LinTrack_prelim\';

%preallocate
%
all_ripple_powers = cell(1,3);


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

        %load ripple band power
        count = 0;
        for istage = [1 3 5]
            count = count+1;
            
            %zscore session ripples during stillness
            all_sesh_rsc_rip_power = datamtx(datamtx(:,11)==1 & ismember(datamtx(:,7),[1 5]), 9);
            mean_rippwr = nanmean(all_sesh_rsc_rip_power);
            std_rippwr = nanstd(all_sesh_rsc_rip_power);
            
            stage_rips_withhpc = datamtx(datamtx(:,11)==1 & datamtx(:,12)>0 & datamtx(:,7)==istage, 9);
            stage_rips_withhpc = (stage_rips_withhpc-mean_rippwr)./std_rippwr;
            
            %load
            all_ripple_powers{count} = [all_ripple_powers{count}; stage_rips_withhpc];
        
        end
    end
end


%plot
figure; hold on
count = 0;
for istage = [1 3 5]
    count = count+1;
    subplot(3,1,count)

    %prepare output
    histogram(all_ripple_powers{count}, -5:0.1:15, 'normalization', 'probability')

    %labels
    xlabel('RSC ripple power (zscored)')
    ylabel(['Prop. HPC-rips Stage' num2str(istage)])

    %aesthetics
    set(gca,'TickLength',[0, 0]); box off
end

