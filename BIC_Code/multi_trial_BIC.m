%[bic, mobic, logL, pt_bic] = model_order_selection_compare(X, momax, time_mode)
%
% X:         variables (channels) * time points * trials numbers
% momax:     the maximum number of scanned model orders (scan from 1 to momax)
% time_mode: 'inhomo','homo'
%
% bic:       (format) momax*2
%            first column: BIC values calculated by our revised BIC
%            second column: BIC values calculated by BIC with hessians
% mobic:     (format) 1*2
%            first column: BIC values calculated by our revised BIC
%            second column: BIC values calculated by BIC with hessians
% logL:      (format) momax*2
%            log likelihood obtained by models of different orders
% pt_bic:    (format) momax*2 
%            first column: penalty term of our revised BIC
%            second column: penalty term of BIC with hessians

% by Kaidi Shao, 06.03.2018, MPI Biological Cybernetics

% estim_mode = 'OLS', 'RLS'
% by Kaidi Shao, 12.01.2022, MPI Biological Cybernetics

%%
function [BICoutputs] = multi_trial_BIC(Yt_events_momax, BICParser)

momax = BICParser.Params.BIC.momax;
[temp, nobs, ntrials] = size(Yt_events_momax);
nvar = temp/(momax+1);

BICoutputs.bic = nan(momax,2);
BICoutputs.pt_bic = nan(momax,2);
BICoutputs.logL = nan(momax,1);
BICoutputs.sum_detHess = nan(momax,1);

for mo = 1:momax

    display(strcat('Start calculation for model order:', int2str(mo)));
    X = squeeze(Yt_events_momax(1:nvar*(mo+1),:,:));
    [BICoutputs.logL(mo), BICoutputs.sum_detHess(mo)] = BIC_compare(X, mo, BICParser);
    
    BICoutputs.pt_bic(mo,1) = .5*nobs*mo*nvar*nvar*log(ntrials);
    BICoutputs.pt_bic(mo,2) = .5*BICoutputs.sum_detHess(mo);
    BICoutputs.pt_bic(mo,3) = .5* nobs*mo*nvar*nvar*log(ntrials*nobs);
    BICoutputs.pt_bic(mo,4) = .5* mo*nvar*nvar*log(ntrials*nobs);

    BICoutputs.bic(mo,1) = - BICoutputs.logL(mo)*ntrials + BICoutputs.pt_bic(mo,1);
    BICoutputs.bic(mo,2) = - BICoutputs.logL(mo)*ntrials + BICoutputs.pt_bic(mo,2); 
    BICoutputs.bic(mo,3) = - BICoutputs.logL(mo)*ntrials + BICoutputs.pt_bic(mo,3);
    BICoutputs.bic(mo,4) = - BICoutputs.logL(mo)*ntrials + BICoutputs.pt_bic(mo,4); 
end


[~, BICoutputs.mobic] = nanmin(BICoutputs.bic);
BICoutputs.mobic = BICoutputs.mobic';


end
  
    