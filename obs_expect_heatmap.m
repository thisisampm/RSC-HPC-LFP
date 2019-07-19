function [div_mtx, observed, expected] = obs_expect_heatmap(vect1, vect2, rng1, rng2, bins, plot_on)
%vect1 = zscore_mtx(datamtx(datamtx(:,11)==1 & ismember(datamtx(:,7), [1 5]),9));
%vect2 = zscore_mtx(datamtx(datamtx(:,11)==1 & ismember(datamtx(:,7), [1 5]),10));
%rng1 = [-2 20]; rng2 = [-2 20]; bins = 50;
%[div_mtx, observed, expected] = obs_expect_heatmap(vect1, vect2, rng1, rng2, bins, 1)


%evenly spaced bins of x and y coordinate ranges
edges1 = linspace(rng1(1), rng1(2)+realmin, bins+1);
edges2 = linspace(rng2(1), rng2(2)+realmin, bins+1);

%2d histogram of event counts
observed = histcounts2(vect1, vect2, edges1, edges2);
expected = histcounts(vect2, edges2).*histcounts(vect1, edges1)';


%remove single events
%observed(observed<2) = 0;

%flip
observed = flipud(observed);
expected = flipud(expected);

%prep for unexpected
expected(expected==0) = realmin;

%divide spikes by time for rate
div_mtx = observed./expected;
div_mtx(isinf(div_mtx)) = realmax;
div_mtx(isnan(div_mtx)) = 0;


%plot
if plot_on == 1
    figure; 
    subplot(1,3,1); 
    imagesc(observed); colorbar; caxis([0 0.0001]); axis square; title observed
    subplot(1,3,2);
    imagesc(expected); colorbar; caxis([0 0.0001]); axis square; title expected
    subplot(1,3,3);
    imagesc(div_mtx); colorbar; caxis([0 100]); axis square; title obs/exp
end