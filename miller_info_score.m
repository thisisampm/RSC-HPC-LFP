function [mis, all_trial_rates] = miller_info_score(datamtx, clusters, bins, varargin)
% computes a score for each cell indicating how strongly the firing rate
% changes in response to the ripple event


if nargin >3
    ploton=varargin{1};
else
    ploton=1;
end

%trials
uniq_trls = unique(datamtx(datamtx(:,7)==3 & datamtx(:,11)==0 & ~isnan(datamtx(:,13)),13))';

%all trial rates
%all_trial_rates = ...
%    nan(length(uniq_trls), bins, size(clusters,1));

all_trial_rates = [];

%for each cell
mis = nan(size(clusters,1),3);
for ic = 1:size(clusters,1)
    clust = clusters(ic,1);
    
    
    %load trial firing rates
    dwell_min = 0.10; %s
    [~, ~, trl_bin_rates, trl_dwell_times] = rate_xdim(datamtx(datamtx(:,7)==3 & datamtx(:,11)==0, :), clust, bins, 0);
    
    %trl_rm_hold = nan(size(trl_bin_rates));
    %trl_rm_hold(trl_dwell_times>dwell_min) = trl_bin_rates(trl_dwell_times>dwell_min);
    trl_rm_hold = trl_bin_rates; trl_rm_hold(trl_dwell_times<dwell_min) = nan;
    %all_trial_rates(:,:,ic) = trl_rm_hold;
    all_trial_rates = cat(2, all_trial_rates, trl_rm_hold);

        
    %calculate mis independently for odd and even trials (running directions)
    odd_even_mis = nan(1,2);
    for ioe = 1:length(odd_even_mis)
        if ioe==1
            oe_idx = rem(uniq_trls,2)==1;
        else
            oe_idx = ~rem(uniq_trls,2)==1;
        end
        local_atr = all_trial_rates(oe_idx,:,ic);
        
        %impose trial minimum
        if size(local_atr,1) < 3
            continue
        end
        
        
        %calculate bin-by-bin variance on each trial (good variance)
        bin_var = nan(length(uniq_trls), 1);
        for it = 1:size(local_atr,1)
            bin_var(it) = nanvar(local_atr(it,:));
        end

        %calculate trial-by-trial variance at each bin (bad variance)
        trial_var = nan(bins, 1);
        for ib = 1:bins
            
            if sum(~isnan(local_atr(:,ib)))<3
                continue
            end
            
            trial_var(ib) = nanvar(local_atr(:,ib));
        end 
        
        odd_even_mis(ioe) = (nanmean(bin_var)-nanmean(trial_var))/...
            (nanmean(bin_var)+nanmean(trial_var));
        
    end

        %load mis
        if sum(isnan(odd_even_mis))<2
            mis(ic, :) = [max(odd_even_mis(~isnan(odd_even_mis))) odd_even_mis];
        end
end

%plot rates
if ploton==1
    for ic = 1:size(all_trial_rates,3)
        figure; hold on; 
        for it = 1:size(all_trial_rates,1) 
            if rem(it,2)==0
                plot(all_trial_rates(it,:,ic), 'color', [     0    0.4470    0.7410]) 
            else
                plot(all_trial_rates(it,:,ic), 'color', [0.8500    0.3250    0.0980])
            end
        end
        axis([0 bins+1 0 inf])
        set(gca,'TickLength',[0, 0]); box off; 
        title(['Cell ' num2str(clusters(ic,1)) ', MIS ' num2str(mis(ic,:))])
    end
end


end