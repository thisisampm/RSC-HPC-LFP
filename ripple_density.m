function [ripples_per_sec] = ripple_density(datamtx)

%hpc before, hpc after, rsc before, rsc after
ripples_per_sec = nan(1, 4);

%stillness index
still_idx = datamtx(:,11)==1;


count = 0;
for cols = [12 14]
    for stages = [1 5]
        count = count+1;
        
        %true ripples
        [hpc_rips, rsc_rips] = true_ripples(datamtx(datamtx(:,7)==stages,:));
        if cols == 12
            num_ripples = length(hpc_rips);
        elseif cols == 14
            num_ripples = length(rsc_rips);
        end
        still_time = length(datamtx(still_idx & datamtx(:,8)==0 & datamtx(:,7)==stages,1))/100;
        
        %minimum stillness time
        if still_time < 60 %s
            continue
        end
        
        ripples_per_sec(count) = num_ripples/still_time;
        
        
    
    
    end
end



end