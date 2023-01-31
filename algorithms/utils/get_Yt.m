function Yt = get_Yt(y, loc, mo, tau, L_start, L_extract)
nvar = size(y,1);
Yt = nan(nvar*(mo+1), L_extract, length(loc));
idx1 = 1:nvar*(mo+1);
idx2 = repmat([1:nvar],[1,mo+1]);
delay = repmat([0:mo]*tau, [nvar,1]);
delay = delay(:);
for n=1:length(idx1)
    Yt(idx1(n),:,:) = extract_events(y(idx2(n),:),loc-delay(n),L_start, L_extract);
end
end