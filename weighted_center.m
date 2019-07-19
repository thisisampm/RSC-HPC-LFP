function wc = weighted_center(vect)
%computes the center of vect according to the magnitude of the values at
%each pixle

vect = (norm_mtx(vect)).*10^5;

vect_prime = []; 
for i = 1:length(vect)
    if ~isnan(round(vect(i)))
        vect_prime = [vect_prime; repmat(i, round(vect(i)),1)]; 
    else
        vect_prime = [vect_prime; NaN]; 
    end
end

wc = mean(vect_prime);

end