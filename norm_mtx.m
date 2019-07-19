function mtx = norm_mtx(mtx)
%subtract by min divide by (new) max
%operates on columns

mtx = (mtx-min(mtx))./max(mtx-min(mtx));
end



