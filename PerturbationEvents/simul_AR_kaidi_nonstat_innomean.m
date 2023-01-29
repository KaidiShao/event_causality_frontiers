function [y]=simul_AR_kaidi_nonstat_innomean(A, SIG, innomean, morder)

[nvar,L] = size(innomean);
y = mvnrnd(zeros(1,2), SIG,L)';
y = y + innomean;
for t = 1:L
    y(:,t) = mvnrnd(innomean(:,t), SIG);
end

for t = morder+1:L
    temp = reshape(fliplr(y(:,t-morder:t-1)),nvar*morder,1);
    y(:,t) = y(:,t)+A*temp;
end

