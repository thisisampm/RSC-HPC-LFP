function [bin_rates, bin_rates_norm, trl_bin_rates, spike_counts] = rate_xdim(datamtx, clusters, varargin)
%plots rate over the xdim of the maze

if nargin == 3
    ploton = varargin{1};
else
    ploton = 1;
end

%position bounds
pos_bnds = [0 1];

%all trials
all_trials = unique(datamtx(~isnan(datamtx(:,13)),13))';
num_trials = length(all_trials);

%rate for each cluster
bin_rates = [];

%spike count matrix (cells,trials)
spike_counts = nan(size(clusters,1), num_trials);
   
%rate on each trial
trl_bin_rates = []; %preallocate trial rows
for itrl = 1:num_trials
    current_trial = all_trials(itrl);

    %local datamtx
    trl_dm = datamtx(datamtx(:,13)==current_trial, :);

    %spike events
    spike_pos = spike_positions(trl_dm, clusters, pos_bnds, 0);
    for isp = 1:length(spike_pos)
        spike_pos{isp} = (spike_pos{isp}-pos_bnds(1))./diff(pos_bnds);
        
        %load spike counts
        spike_counts(isp,itrl) = length(spike_pos{isp});
    end
    
    %time (vid) events
    dwell_pos = dwell_positions(trl_dm, pos_bnds, 0);
    for idp = 1:length(dwell_pos)
        dwell_pos{idp} = (dwell_pos{idp}-pos_bnds(1))./diff(pos_bnds);
    end
    
    %bin events
    pos_bnds = [0 1];
    spike_pos_binned = boxcar_rates(spike_pos, pos_bnds, 0.1, 0.01);
    dwell_pos_binned = boxcar_rates(dwell_pos, pos_bnds, 0.1, 0.01)./100;

    %rate in each bin
    trl_bin_rates = cat(3, trl_bin_rates, spike_pos_binned./dwell_pos_binned);

end

%sum spike counts over all trials
spike_counts = sum(spike_counts,2);

%normalize rates
bin_rates = nanmean(trl_bin_rates,3);  
bin_rates_norm = norm_mtx(bin_rates')';
%bin_rates_norm = zscore_mtx(bin_rates')';

if ploton == 1
    
    %single plot
    if length(clusters) == 1

        %plot each trial
        figure; hold on
        new_x = 1:size(trl_bin_rates(1,:,itrl),2);
        
        for itrl = 1:num_trials

            %plot each direction as different color
            if rem(itrl,2)==0
                plot(new_x, trl_bin_rates(1,:,itrl), 'color', [0    0.4470    0.7410]);
            else
                plot(new_x, trl_bin_rates(1,:,itrl), 'color', [0.8500    0.3250    0.0980]);
            end
        end
        
        %plot trial means
        for idir = 1:2
            
            %directional rates
            if idir==1
                dir_rates = trl_bin_rates(1, :, logical(rem(1:num_trials,2)));
            elseif idir==2
                dir_rates = trl_bin_rates(1, :, ~logical(rem(1:num_trials,2)));
            end
            dir_rates_mean = nanmean(dir_rates,3);
            dir_rates_std = nanstd(dir_rates, [], 3);
            
            %smooth
            dir_rates_mean = smooth(dir_rates_mean, 5);
            dir_rates_std = smooth(dir_rates_std, 5);
            
            % +/- 1 se
            dir_rates_se_hi = dir_rates_mean + dir_rates_std./sqrt(size(dir_rates,3));
            dir_rates_se_lo = dir_rates_mean - dir_rates_std./sqrt(size(dir_rates,3));

            %plot each direction as different color
            if idir==2
                plot(new_x, dir_rates_mean, 'color', [0    0.4470    0.7410], 'linewidth', 3);
                plot(new_x, dir_rates_se_hi, 'color', [0    0.4470    0.7410], 'linewidth', 2);
                plot(new_x, dir_rates_se_lo, 'color', [0    0.4470    0.7410], 'linewidth', 2);
            else
                plot(new_x, dir_rates_mean, 'color', [0.8500    0.3250    0.0980], 'linewidth', 3);
                plot(new_x, dir_rates_se_hi, 'color', [0.8500    0.3250    0.0980], 'linewidth', 2);
                plot(new_x, dir_rates_se_lo, 'color', [0.8500    0.3250    0.0980], 'linewidth', 2);
            end
        end

        
        %aesthetics
        xlabel('X position')
        ylabel('Firing Rate (Hz)')
        set(gca, 'TickLength', [0 0])
        set(gcf, 'Position', [600 450 800 300])
        title(['Cluster ' num2str(clusters)])
        xlim([1 size(dir_rates,2)])
        xticks(linspace(1,size(dir_rates,2),11))
        xticklabels(linspace(0,1,11))
    
    %multicell rate matrix
    else
        
        figure; 
        imagesc(bin_rates_norm);
        
        %aesthetics
        set(gca, 'TickLength', [0 0])
        xticks(linspace(1,size(bin_rates_norm,2),11))
        xticklabels(linspace(0,1,11))
        yticks(1:length(clusters))
        yticklabels(clusters)
        xlabel('X position')
        ylabel('Cell IDs')
        colorbar; caxis([0 1])

    end
end


end