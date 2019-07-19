function samprate = cscmtx_samprate(csc_mtx)
%outputs most likely sample rate from a neuralynx session

nlyx_default = 32000;
plausible_sample_rates = nlyx_default./(1:1000);
samprate = round(size(csc_mtx,1)/csc_mtx(end,1));
samprate = plausible_sample_rates(abs(plausible_sample_rates-samprate)==min(abs(plausible_sample_rates-samprate)));

end