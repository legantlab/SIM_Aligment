function [r_id_fit,c_id_fit] = fit_2D_phasemap(phase_map, r_id,c_id,winsize)

r_min = max(r_id - winsize,1);
r_max = min(r_id + winsize,size(phase_map,1));
if mod(r_max - r_min,2) == 1
    if r_min == 1
        r_max = r_max+1;
    elseif r_max == size(phase_map,1)
        r_min = r_min-1;
    end
end
c_min = max(c_id - winsize,1);
c_max = min(c_id + winsize,size(phase_map,2));
if mod(c_max - c_min,2) == 1
    if c_min == 1
        c_max = c_max+1;
    elseif c_max == size(phase_map,2)
        c_min = c_min-1;
    end
end
CC_subregion=phase_map(r_min:r_max,c_min:c_max);
CC_subregion = abs(CC_subregion - max(CC_subregion(:)));
CC_subregion=CC_subregion/max(CC_subregion(:));

[~,min_coord] = max(CC_subregion(:));
[r_guess,c_guess]=ind2sub(size(CC_subregion),min_coord);


%Update our shift vector guess based on the results of the Gaussian fit
FitData = Lsq_Fit_2D_2020_09_21(double(CC_subregion),r_guess,c_guess);
r_id_fit = ceil((r_min+r_max)/2) + FitData(3);
c_id_fit = ceil((c_min+c_max)/2) + FitData(2);
end
