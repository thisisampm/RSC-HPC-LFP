function [all_ripples_per_sec, all_joint_dens, all_joint_props] = ALL_ripple_density

%should be input
folderpath = 'C:\Users\ampm1\Documents\MATLAB\tt_ephys\neurodata_new';

%preallocate
all_ripples_per_sec = [];
all_joint_dens = [];
all_joint_props = [];
all_ind_pre = [];
all_joint_pre = [];
all_ind_post = [];
all_joint_post = [];

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

        %compute ripple density
        [ripples_per_sec] = ripple_density(datamtx);
        
        %compute joint ripple number and proportion
        [jct_pre, prop_pre, hpc_venn_pre, rsc_venn_pre, joint_venn_pre] = joint_rip_prop(datamtx(datamtx(:,7)==1, :));
        [jct_post, prop_post, hpc_venn_post, rsc_venn_post, joint_venn_post] = joint_rip_prop(datamtx(datamtx(:,7)==5, :));
        
        %still time
        still_time_pre = length(datamtx(datamtx(:,11)==1 & datamtx(:,8)==0 & datamtx(:,7)==1,1))/100;
        still_time_post = length(datamtx(datamtx(:,11)==1 & datamtx(:,8)==0 & datamtx(:,7)==5,1))/100;
        
        %concatonate
        all_ripples_per_sec = [all_ripples_per_sec; ripples_per_sec];

        all_joint_dens = [all_joint_dens; [jct_pre./still_time_pre  jct_post./still_time_post]];
        all_joint_props = [all_joint_props; [prop_pre prop_post]];

        
        %venn diagram raw numbers       
        all_ind_pre = [all_ind_pre; [hpc_venn_pre rsc_venn_pre]];
        all_joint_pre = [all_joint_pre; joint_venn_pre];
        all_ind_post = [all_ind_post; [hpc_venn_post rsc_venn_post]];
        all_joint_post = [all_joint_post; joint_venn_post];
        
    end
end


%plot density
for i = 1:4; cell_e_density{i} = all_ripples_per_sec(:,i); end

%{
errorbar_plot(cell_e_density(1:2))
[~, b, ~, d] = ttest(all_ripples_per_sec(:,1), all_ripples_per_sec(:,2))
title(['Density HPC, p=' num2str(b) ', t=' num2str(abs(d.tstat))])
xticklabels({'Before', 'After'})

errorbar_plot(cell_e_density(3:4))
[~, b, ~, d] = ttest(all_ripples_per_sec(:,3), all_ripples_per_sec(:,4))
title(['Density RSC, p=' num2str(b) ', t=' num2str(abs(d.tstat))])
xticklabels({'Before', 'After'})
%}

figure; [r_dense, p_dense] = fit_line(cell_e_density{2}-cell_e_density{1}, cell_e_density{4}-cell_e_density{3})
xlabel('HPC after - before')
ylabel('RSC after - before')
set(gca,'TickLength',[0, 0]); box off;
title(['Change in ripple density, r=' num2str(r_dense) ', p=' num2str(p_dense)])

%plot joint_num
for i = 1:2; cell_e_jointdense{i} = all_joint_dens(:,i); end
errorbar_plot(cell_e_jointdense([1 2]))
[~, b, ~, d] = ttest(all_joint_dens(:,1), all_joint_dens(:,2))
title(['Jointdense, p=' num2str(b) ', t=' num2str(abs(d.tstat))])
xticklabels({'Before', 'After'})

%plot joint_prop
for i = 1:2; cell_e_jointprop{i} = all_joint_props(:,i); end
errorbar_plot(cell_e_jointprop([1 2]))
[~, b, ~, d] = ttest(all_joint_props(:,1), all_joint_props(:,2))
title(['Jointprop, p=' num2str(b) ', t=' num2str(abs(d.tstat))])
xticklabels({'Before', 'After'})

%venn diagram
figure;
subplot(1,2,1)
venn(nanmean(all_ripples_per_sec(:, [1 3])), nanmean(all_joint_dens(:,1)))
set(gca,'TickLength',[0, 0]); box off;
axis equal
ahold = axis;
axis(ahold)
subplot(1,2,2)
venn(nanmean(all_ripples_per_sec(:, [2 4])), nanmean(all_joint_dens(:,2)))
set(gca,'TickLength',[0, 0]); box off;
axis(ahold)

%correlation between regional ripple density and joint ripple density
mean_region_ripdense = mean([[all_ripples_per_sec(:,2);all_ripples_per_sec(:,1)] [all_ripples_per_sec(:,4);all_ripples_per_sec(:,3)]],2);
joint_ripdense = [all_joint_dens(:,2);all_joint_dens(:,1)] ;
figure; [r, p] = fit_line(mean_region_ripdense, joint_ripdense)
