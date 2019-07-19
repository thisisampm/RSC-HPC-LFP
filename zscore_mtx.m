function mtx_z = zscore_mtx(mtx)


%zscore
mtx = mtx - repmat(nanmean(mtx), size(mtx,1), 1);
mtx_z = mtx./nanstd(mtx);

end