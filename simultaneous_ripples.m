function [rip_nums, rip_times] = simultaneous_ripples(datamtx, varargin)
%identify ripple ids for simultaneous ripples
%output vects: rsc, hpc
%times are mins
%
%second input should be a 1 for a neat plot

%stage index
stage_idx = ismember(datamtx(:,7), [1 5]);
stillness_idx = datamtx(:,11)==1;

%ripple flag numbers
rsc_rips = datamtx(stage_idx & stillness_idx,14);
hpc_rips = datamtx(stage_idx & stillness_idx,12);

%rip number output
rip_nums = unique([rsc_rips(rsc_rips>0 & hpc_rips>0) hpc_rips(rsc_rips>0 & hpc_rips>0)], 'rows');

%rip times output
region_cols = [14, 12];
rip_times = nan(size(rip_nums));
for iregion = 1:size(rip_times,2)
    for irip = 1:size(rip_times,1)
        current_rip_num = rip_nums(irip, iregion);     
        rip_times(irip, iregion) = min(datamtx(datamtx(:,region_cols(iregion))==current_rip_num  & datamtx(:,8)==0, 1));
    end
end

%plot simultaneous ripple events
if ~isempty(varargin)
    if varargin{1} == 1
        rsc_rips = datamtx(stillness_idx,14);
        hpc_rips = datamtx(stillness_idx,12);
        hpc_rips_norm = hpc_rips./max(hpc_rips);
        rsc_rips_norm = rsc_rips./max(rsc_rips);

        figure; hold on
        plot(hpc_rips_norm, rsc_rips_norm,'o')
        xlabel('hpc rips')
        ylabel('rsc rips')
        plot([1 1].*max(hpc_rips_norm(datamtx(stillness_idx,7)==1)), ylim, 'k-')
        plot([1 1].*min(hpc_rips_norm(datamtx(stillness_idx,7)==5 & hpc_rips_norm~=0)), ylim, 'k-')
        plot(xlim, [1 1].*max(rsc_rips_norm(datamtx(stillness_idx,7)==1)), 'k-')
        plot(xlim, [1 1].*min(rsc_rips_norm(datamtx(stillness_idx,7)==5 & rsc_rips_norm~=0)), 'k-')

        set(gca,'TickLength',[0, 0]); box off; axis square;
    end
end