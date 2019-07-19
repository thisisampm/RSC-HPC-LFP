function [shuf_out, obs_out] = shuffle_all_replay_scores_r(rm_out, shuffs)

%preallocate
obs_out = nan(1, 2);
shuf_out = nan(shuffs, 2);


figure
for rsc_hpc = 1:2

    %obs
    rms = rm_out{1,rsc_hpc}; 
    rmr = rm_out{2,rsc_hpc}; 
    allrs = nan(size(rms(:,1))); 
    for ic = 1:size(rms,1)
        allrs(ic) = max([corr(rms(ic,:)', rmr(ic,:)') corr(rms(ic,:)', fliplr(rmr(ic,:))')]); 
        %allrs(ic) = corr(rms(ic,:)', rmr(ic,:)');
        %allrs(ic) = corr(rms(ic,:)', fliplr(rmr(ic,:))');
    end
    obs_out(rsc_hpc) = mean(allrs);

    %shuffs
    for ishuf = 1:shuffs

        %shuff
        rmr = rmr(randperm(size(rmr,1)),:);

        %calculate
        allrs = nan(size(rms(:,1))); 
        for ic = 1:size(rms,1)
            allrs(ic) = max([corr(rms(ic,:)', rmr(ic,:)') corr(rms(ic,:)', fliplr(rmr(ic,:))')]); 
            %allrs(ic) = corr(rms(ic,:)', rmr(ic,:)');
            %allrs(ic) = corr(rms(ic,:)', fliplr(rmr(ic,:))');
        end

        shuf_out(ishuf, rsc_hpc) = mean(allrs);
    end

    %plot
    subplot(2,1, rsc_hpc); hold on
    histogram(shuf_out(:,rsc_hpc), -1:.001:1, 'normalization', 'probability')
    plot([1 1].*obs_out(rsc_hpc), ylim, 'r-')
    plot([0 0], ylim, 'k--')
    title(num2str(rsc_hpc))
end

