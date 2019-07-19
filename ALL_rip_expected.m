function [all_div, all_obs, all_exp] = ALL_rip_expected
% combine ripple and cell events details across days
% plot

%should be input
folderpath = 'C:\Users\ampm1\Documents\MATLAB\tt_ephys\neurodata_new\';

%preallocate
%
all_obs = [];
all_exp = [];
all_div = [];


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
        
        %bin ripple power
        binsize = 5; %10ms/bin]
        bindx = datamtx(:,8)==0 & datamtx(:,11)==1 & ismember(datamtx(:,7), [1 5]);
        RSC_rbp = datamtx(bindx,9);
            RSC_rbp = RSC_rbp(1:end-rem(size(RSC_rbp,1), binsize));
            RSC_rbp = reshape(RSC_rbp, binsize, size(RSC_rbp,1)/binsize);
            RSC_rbp = mean(RSC_rbp);
            RSC_rbp = zscore_mtx(RSC_rbp);
        HPC_rbp = datamtx(bindx,10);
            HPC_rbp = HPC_rbp(1:end-rem(size(HPC_rbp,1), binsize));
            HPC_rbp = reshape(HPC_rbp, binsize, size(HPC_rbp,1)/binsize);
            HPC_rbp = mean(HPC_rbp);
            HPC_rbp = zscore_mtx(HPC_rbp);
        
        %compute
        rng1 = [-5 30]; rng2 = rng1;
        heatmap_bins = 72;
        %RSC_rbp_still = zscore_mtx(datamtx(datamtx(:,11)==1 & ismember(datamtx(:,7), [1 5]), 9));
        %HPC_rbp_still = zscore_mtx(datamtx(datamtx(:,11)==1 & ismember(datamtx(:,7), [1 5]), 10));
        [div_mtx, observed, expected] = obs_expect_heatmap(RSC_rbp, HPC_rbp, rng1, rng2, heatmap_bins, 0);

        %load ouput
        all_obs = cat(3,all_obs, observed);
        all_exp = cat(3,all_exp, expected);
        all_div = cat(3,all_div, div_mtx);

    end
end

%compute overall
sum_all_obs = sum(all_obs,3); 
sum_all_exp = sum(all_exp,3); 
%remove zero probabilities
sum_all_exp(sum_all_exp<=0) = realmin;
%normalize
sum_all_obs = sum_all_obs./nansum(sum_all_obs(:));
sum_all_exp = sum_all_exp./nansum(sum_all_exp(:));
%set to percent
sum_all_obs = sum_all_obs.*100;
sum_all_exp = sum_all_exp.*100;
%divide
sum_all_div = sum_all_obs./sum_all_exp;

%plot
figure; 
subplot(1,3,1); hold on;
imagesc(sum_all_obs); colorbar; axis square; title observed
xticks(interp1(linspace(rng1(1), rng1(2), 7),linspace(1, size(sum_all_obs,2), 7),floor(linspace(rng1(1), rng1(2), 7))))
xticklabels(floor(linspace(rng1(1), rng1(2), 7)))
xlabel('HPC ripple power (Z)')
yticks(sort(size(sum_all_obs,2) +1 - interp1(linspace(rng1(1), rng1(2), 7),linspace(1, size(sum_all_obs,2), 7), floor(linspace(rng1(1), rng1(2), 7)))))
yticklabels(fliplr(floor(linspace(rng1(1), rng1(2), 7))))
ylabel('RSC ripple power (Z)')
set(gca,'ColorScale','log'); caxis([0.0001 100])

    %add dotted lines at zeros
    zero_tick_x = interp1(linspace(rng1(1), rng1(2), 7),linspace(1, size(sum_all_obs,2), 7), 0);
    zero_tick_y = size(sum_all_obs,2)-zero_tick_x+1;
    plot([1 1].*zero_tick_x, ylim, 'k--')
    plot(xlim, [1 1].*zero_tick_y, 'k--')
    axis([.5 heatmap_bins+.5 .5 heatmap_bins+.5])
    set(gca,'Ydir','reverse')
    
subplot(1,3,2); hold on
imagesc(sum_all_exp); colorbar; axis square; title expected
xticks(interp1(linspace(rng1(1), rng1(2), 7),linspace(1, size(sum_all_obs,2), 7),floor(linspace(rng1(1), rng1(2), 7))))
xticklabels(floor(linspace(rng1(1), rng1(2), 7)))
xlabel('HPC ripple power (Z)')
yticks(sort(size(sum_all_obs,2) +1 - interp1(linspace(rng1(1), rng1(2), 7),linspace(1, size(sum_all_obs,2), 7), floor(linspace(rng1(1), rng1(2), 7)))))
yticklabels(fliplr(floor(linspace(rng1(1), rng1(2), 7))))
ylabel('RSC ripple power (Z)')
set(gca,'ColorScale','log'); caxis([0.0001 100])

    %add dotted lines at zeros
    zero_tick_x = interp1(linspace(rng1(1), rng1(2), 7),linspace(1, size(sum_all_obs,2), 7), 0);
    zero_tick_y = size(sum_all_obs,2)-zero_tick_x+1;
    plot([1 1].*zero_tick_x, ylim, 'k--')
    plot(xlim, [1 1].*zero_tick_y, 'k--')
    axis([.5 heatmap_bins+.5 .5 heatmap_bins+.5])
    set(gca,'Ydir','reverse')

subplot(1,3,3); hold on
imagesc(sum_all_div); colorbar; axis square; title obs/exp
xticks(interp1(linspace(rng1(1), rng1(2), 7),linspace(1, size(sum_all_obs,2), 7),floor(linspace(rng1(1), rng1(2), 7))))
xticklabels(floor(linspace(rng1(1), rng1(2), 7)))
xlabel('HPC ripple power (Z)')
yticks(sort(size(sum_all_obs,2) +1 - interp1(linspace(rng1(1), rng1(2), 7),linspace(1, size(sum_all_obs,2), 7), floor(linspace(rng1(1), rng1(2), 7)))))
yticklabels(fliplr(floor(linspace(rng1(1), rng1(2), 7))))
ylabel('RSC ripple power (Z)')
set(gca,'ColorScale','log'); caxis([0.0001 100])

    %add dotted lines at zeros
    zero_tick_x = interp1(linspace(rng1(1), rng1(2), 7),linspace(1, size(sum_all_obs,2), 7), 0);
    zero_tick_y = size(sum_all_obs,2)-zero_tick_x+1;
    plot([1 1].*zero_tick_x, ylim, 'k--')
    plot(xlim, [1 1].*zero_tick_y, 'k--')
    axis([.5 heatmap_bins+.5 .5 heatmap_bins+.5])
    set(gca,'Ydir','reverse')

