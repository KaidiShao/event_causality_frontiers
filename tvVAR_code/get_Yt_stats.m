function Yt_stats_cond = get_Yt_stats(Yt_event,mo)

nvar = size(Yt_event,1)/(mo+1);
Yt_stats_cond.mean = mean(Yt_event,3);
Yt_stats_cond.Ntrials = size(Yt_event,3);
for t = 1:size(Yt_event, 2)
    temp = squeeze(Yt_event(:,t,:))-Yt_stats_cond.mean(:,t);
    Yt_stats_cond.Sigma(t,:,:) = temp*temp'/size(Yt_event,3);
    Yt_stats_cond.OLS.At(t,:,:) = reshape(squeeze(Yt_stats_cond.Sigma(t,1:nvar,nvar+1:end)),nvar,nvar*mo)...
                                 /squeeze(Yt_stats_cond.Sigma(t,nvar+1:end,nvar+1:end));
end

[Yt_stats_cond.OLS.bt, Yt_stats_cond.OLS.Sigma_Et, Yt_stats_cond.OLS.sigma_Et] = estimate_residuals(Yt_stats_cond);

end