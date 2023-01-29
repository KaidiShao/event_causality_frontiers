
%%
function [X, Imp] = gen_ensemble_nonstat_innomean(A, SIG, ntrials, L_event, center, amp, dim, L_perturb)

[nvar, temp] = size(A);
morder = temp/nvar;

X = nan(nvar, L_event, ntrials);
Imp = nan(nvar, L_event, ntrials);

for n = 1:ntrials % for each realization
    [X(:,:,n), Imp(:,:,n)] = genvar_nonstat(A,SIG,morder,nvar,L_event,amp,dim,L_perturb,center);
end

end
function [X, Imp] = genvar_nonstat(A,SIG,morder,nvar,L_event,amp,dim,L_perturb,center)

%% initialise to Gaussian white noise
X = SIG*randn(nvar, L_event+morder); % "SIG" is actually Cholesky matrix
if L_perturb==1
    Imp_shape = amp*1;
else
    Imp_shape = amp*morlet(-4,4,L_perturb); % 101 point long morlet wave
end
    Imp = zeros(2,L_event);
    Imp(dim, center-floor(L_perturb/2):center+floor(L_perturb/2)) = Imp_shape;

for t = morder+1:L_event   % for each time step
    X_lag = reshape(fliplr(X(:,t-morder:t-1)),nvar*morder,1);
    X(:,t) = X(:,t) + A*X_lag;
    X(:,t) = X(:,t) + Imp(:,t-morder);
end
X(:,1:morder)=[];
end

