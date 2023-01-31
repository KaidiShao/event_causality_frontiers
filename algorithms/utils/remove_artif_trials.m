% temporally only to remove artifect for opendata ca3
function [Yt_events, locs] = remove_artif_trials(Yt_events, locs, lower_thres)

idxall = find(Yt_events(1:2,:,:)<lower_thres);
[~, ~, itrial] = ind2sub(size(Yt_events(1:2,:,:)),idxall);

itrial_remove = unique(itrial);

Yt_events(:,:,itrial_remove) = [];
locs(itrial_remove) = [];

disp(strcat('removed:',int2str(length(itrial_remove)),' artifect trials'));
end