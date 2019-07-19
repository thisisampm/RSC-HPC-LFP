function [hpc_rips, rsc_rips] = true_ripples(datamtx)
%finds the number of each ripple occuring exclusively within stillness

hpc_rips = unique(datamtx(datamtx(:,12)>0 & datamtx(:,11)==1, 12));
hpc_rips_hold = hpc_rips;
rip_ct = 0;
for rip = hpc_rips'
    rip_ct = rip_ct+1;
    if mean(datamtx(datamtx(:,12)==rip, 11))<1
       hpc_rips_hold(rip_ct) = nan;
   end
end
hpc_rips = hpc_rips_hold(~isnan(hpc_rips_hold));

rsc_rips = unique(datamtx(datamtx(:,14)>0 & datamtx(:,11)==1, 14));
rsc_rips_hold = rsc_rips;
rip_ct = 0;
for rip = rsc_rips'
    rip_ct = rip_ct+1;
    if mean(datamtx(datamtx(:,14)==rip, 11))<1
       rsc_rips_hold(rip_ct) = nan;
   end
end
rsc_rips = rsc_rips_hold(~isnan(rsc_rips_hold));