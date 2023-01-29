function A_event = extract_events(A, cumP, L_start, L)

A_event = nan(L, length(cumP));

for i = 1:length(cumP)
    idx = cumP(i)-L_start+1:cumP(i)+L-L_start;
    A_event(:,i)=A(idx);
end
end