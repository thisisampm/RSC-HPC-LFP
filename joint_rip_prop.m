function [joint_count, joint_prop, hpc_venn, rsc_venn, joint_venn] = joint_rip_prop(datamtx)
%compute proportion of hpc ripples that are accompanied by rsc ripples and
%vice versa
%

%props is proportion of ripples that occur within a joint-ripple

%minimum stillness time
joint_prop = nan;
joint_count = nan;
still_time = length(datamtx(datamtx(:,11)==1 & datamtx(:,8)==0,1))/100;
if still_time < 60 %s
    hpc_venn = nan;
    rsc_venn = nan;
    joint_venn = nan;
    return
end

%ripple numbers
[hpc_rips, rsc_rips] = true_ripples(datamtx);

%venn diagram
hpc_venn = [];
joint_venn = [];

num_rips = length(hpc_rips);
for hpc_rip = hpc_rips'
    
   time_buffer = 0.05;%s
   lo_time = min(datamtx(datamtx(:,12)==hpc_rip, 1)) - time_buffer;
   hi_time = max(datamtx(datamtx(:,12)==hpc_rip, 1)) + time_buffer;
   time_idx = datamtx(:,1)>=lo_time & datamtx(:,1)<=hi_time;
    
   if any(ismember(datamtx(time_idx, 14), rsc_rips))
       rsc_rip = mode(datamtx(time_idx & datamtx(:,14)>0, 14));
       joint_venn = [joint_venn; [hpc_rip rsc_rip]];
   else
       hpc_venn = [hpc_venn; hpc_rip];
   end
end
if ~isempty(joint_venn)
    rsc_venn = length(setdiff(rsc_rips, joint_venn(:,2)));
else
    rsc_venn = length(rsc_rips);
end
joint_venn = length(joint_venn);
hpc_venn = length(hpc_venn);

joint_count(1) = size(joint_venn,1);
joint_prop(1) = numel(joint_venn)/(length(hpc_rips)+length(rsc_rips)); 


%{
%hpc rips with rsc
num_rips = length(hpc_rips);
joint_rip_ct = 0;
for rip = hpc_rips'
   if sum(datamtx(datamtx(:,12)==rip, 14))>0
       joint_rip_ct = joint_rip_ct+1;
   end
end
joint_counts(1) = joint_rip_ct;
joint_props(1) = joint_rip_ct/num_rips;


%rsc rips with hpc
num_rips = length(rsc_rips);
joint_rip_ct = 0;
for rip = rsc_rips'
   if sum(datamtx(datamtx(:,14)==rip, 12))>0
       joint_rip_ct = joint_rip_ct+1;
   end
end
joint_counts(2) = joint_rip_ct;
joint_props(2) = joint_rip_ct/num_rips;
%}