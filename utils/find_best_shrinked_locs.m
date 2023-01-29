function [best_locs, distance] = find_best_shrinked_locs(D, shrinked_locs, all_locs)
[Nfull] = histcounts(D(all_locs),100,'Normalization','pdf'); % total distribution of all points D>=d0

distance = nan(1,length(shrinked_locs));
for n = 100:length(shrinked_locs)
    [Ntemp] = histcounts(D(shrinked_locs(1:n)),100,'Normalization','pdf'); % total distribution of all points D>=d0
    distance(n) = pdist2(Nfull,Ntemp);
end

[~, best_n] = min(distance);
best_locs = shrinked_locs(1:best_n);
end
    