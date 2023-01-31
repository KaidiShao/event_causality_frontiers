function [shrinked_locs] = shrink_locs_resample_uniform(loc, L)
%%
Ngen = 1;
maxNgen = length(loc);
shrinked_locs = nan(1,maxNgen);
loc_range = loc;

while Ngen < maxNgen
    if isempty(loc_range)
        break;
    end
    shrinked_locs(Ngen) = loc_range(randi(length(loc_range)));
    loc_range(abs(shrinked_locs(Ngen)-loc_range)<L)=[];

    Ngen = Ngen+1;
end

shrinked_locs(isnan(shrinked_locs))=[];
end
