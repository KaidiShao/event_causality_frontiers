function [bt, Sigma_Et, sigma_Et] = estimate_residuals(Yt_stats)

[L, nvar, temp] = size(Yt_stats.OLS.At);

bt = nan(nvar, L);
Sigma_Et = nan(L, nvar, nvar); 
sigma_Et = nan(L, 1); 
for t = 1:L
    
    Sigma_Xt = squeeze(Yt_stats.Sigma(t,1:nvar,1:nvar));
    Sigma_Xp = squeeze(Yt_stats.Sigma(t,nvar+1:end,nvar+1:end));
    Sigma_XtXp = reshape(squeeze(Yt_stats.Sigma(t,1:nvar,nvar+1:end)), nvar, temp);
    coeff = reshape(squeeze(Yt_stats.OLS.At(t,:,:)), nvar, temp);
    
    bt(:,t) = Yt_stats.mean(1:nvar,t)-coeff*Yt_stats.mean(nvar+1:end,t);
    Sigma_Et(t,:,:) = Sigma_Xt - Sigma_XtXp*coeff' - coeff*Sigma_XtXp' +coeff*Sigma_Xp*coeff';
    sigma_Et(t,1) = trace(squeeze(Sigma_Et(t,:,:)));
end

% Sigma_Et(Sigma_Et<0)=0;
% sigma_Et(sigma_Et<0)=0;

end
