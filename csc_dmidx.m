function cscidx = csc_dmidx(csc_mtx, datamtx, dmidx)
%finds index for csc_mtx corresponding to the input datamtx idx 'dmidx'
%idx must be vector of same length as datamtx(:,1)

%set to index
dmidx_hold = dmidx>0;


%interpolate idx to match csc times
[unq_dmtime, unq_dmtime_idx] = unique(datamtx(:,1));
unq_dmidx = dmidx_hold(unq_dmtime_idx);

nnan_idx = ~isnan(unq_dmtime) & ~isnan((unq_dmidx));

cscidx = interp1(unq_dmtime(nnan_idx), double(unq_dmidx(nnan_idx)), csc_mtx(:,1));

%logic cutoff
cscidx = cscidx > 0.5;

%plot check
%figure; plot(unq_dmtime, unq_dmidx)
%hold on; plot(csc_mtx(:,1), cscidx)


if any(dmidx>1)
   
    cscidx = double(cscidx);
    
    orig_values = [];
    cur_val = nan;
    for i = 1:length(dmidx)
        if dmidx(i) == cur_val
            continue
        else
            cur_val = dmidx(i);
            orig_values = [orig_values; cur_val];
        end
    end
    
    count = 0;
    cur_val = nan;
    for i = 1:length(cscidx)
        if cscidx(i) == cur_val
            cscidx(i) = orig_values(count);
        else
            cur_val = cscidx(i);
            count = count+1;
            cscidx(i) = orig_values(count);
        end
    end  

end
    
    
    
end


