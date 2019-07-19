function errorbar_plot( cell_e )
%plots errorbar

%preallocate
all_mtx = cell2mat(cell_e(:));
all_mtx_grps = [];

means = nan(length(cell_e),1);
stds = nan(length(cell_e),1);
sqrt_l = nan(length(cell_e),1);


figure; 
hold on;

for connections = 1:length(cell_e{1})
    
    line_points = nan(length(cell_e));
    for ilp = 1:length(cell_e)
       line_points(ilp) = cell_e{ilp}(connections); 
    end
    
    plot(1:length(line_points), line_points, '-', 'color', [.7 .7 .7])
    
end


for imtx = 1: length(cell_e)
  
    plot(imtx, cell_e{imtx}, 'ko')
    
    all_mtx_grps = [all_mtx_grps; repmat(imtx, size(cell_e{imtx}(:)))];
    means(imtx) = nanmean(cell_e{imtx});
    stds(imtx) = nanstd(cell_e{imtx});
    sqrt_l(imtx) = sqrt(length(~isnan(cell_e{imtx})));

    
end

std_es = stds./sqrt_l;


set(gca,'TickLength',[0, 0]); box off
xlim([.5 length(cell_e)+.5])
xticks(1:4)


    
    errorbar(means, std_es, 'k')
    
    %{
    anov_mtx = [];
    for i = 1:length(cell_e)
        if size(cell_e{i},1) < size(cell_e{i},2)
            cell_e{i} = cell_e{i}';
        end
        if size(cell_e{i},2)>1
            error('cell_e items must be 1d vectors')
        end
        anov_mtx = [anov_mtx cell_e{i}];
    end
    %}

%[P,ANOVATAB] = anova1(anov_mtx)
%anovan(all_mtx, all_mtx_grps)


end

