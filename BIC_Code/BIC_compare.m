% by Kaidi Shao, 06.03.2018, MPI Biological Cybernetics
% revised by Kaidi Shao, 12.01.2022, Shanghai

function [logL, sum_detHess] = BIC_compare(Yt_events, morder, BICParser)

[temp, nobs, ntrials] = size(Yt_events);
nvar = temp/(morder+1);

switch BICParser.Params.BIC.mode
    case 'biased'
        Yt_stats = get_Yt_stats(Yt_events, morder);
    case 'debiased'
        Params = BICParser.Params;
        Params.Options.BIC = 0;
        Params.Options.save_flag = 0;
        Params.BIC.morder = morder;
        [SnapAnalyOutput] = snapshot_detect_analysis_pipeline(BICParser.OriSignal, BICParser.DetSignal, Params);
        Yt_stats = SnapAnalyOutput.Yt_stats_debiased;
end
% Yt_stats_cond = get_RLS_estimates_Hesse(Yt_events, Yt_stats_cond, morder);


for t = 1:nobs

    C_0 = squeeze(Yt_stats.Sigma(t,nvar+1:end,nvar+1:end)); 
    switch BICParser.EstimMode
    case 'OLS'
        DSIG(t) = prod(diag(squeeze(Yt_stats.OLS.Sigma_Et(t,:,:))));
    case 'RLS'
        DSIG(t) = prod(diag(squeeze(Yt_stats.RLS.Sigma_Et(t,:,:))));
    end
    log_detHess(t) = morder*nvar^2*log(ntrials)+ nvar*log(det(C_0)) - nvar*morder*log(DSIG(t));
    
end

logL = -0.5*nobs*nvar*log(2*pi) - 0.5*sum(log(DSIG))- 0.5*nobs*nvar;
sum_detHess = sum(log_detHess); 

end




%         C = C_0 * ntrials;
%         detHess(t) = det(C)^nvar * (1/DSIG(t))^(nvar*morder);

%%
% if strcmp(time_mode,'homo')
%     C_0 = squeeze(mean(Ct_0,1));
%     C_jcr = squeeze(mean(Ct_j,1));
%     C_1rcr = squeeze(mean(Ct_1r,1));
%     
%     A = reshape(C_jcr,nvar,nvar*r)/C_1rcr;
%     
%     SIG = C_0-A*C_1rcr*A';
%     
%     sumsign = sum(sign(SIG(:)));
%     DSIG = prod(diag(SIG));
%     
%     C = C_1rcr * T * ntrials;
%     detHess = det(C)^nvar * (1/DSIG)^(nvar*morder);
%     
%     logL = - 0.5*T*nvar*log(2*pi) - 0.5*T*(log(DSIG)) - 0.5*T*nvar;
% %     sum_detHess = log(detHess); 
%     sum_detHess = nvar*log(det(C)) + (nvar*morder)*log((1/DSIG)); % do log first to avoid detHess goes to infinity
% end





    