% [Cs, TE] = TE_CS_nonzero_mean(Yt, morder, coeff_t, mean_Yt, cov_Yt)
% 
% Yt in the format of 2(mo+1)*time*ntrial
% by Kaidi Shao, 06.03.2020, MPI Biological Cybernetics

function [CausalOutput] = time_varying_causality(Yt_event, Yt_stats, CausalParams)

%% at the beginning Xt means the bi-variate siganl, Yt means [Xt; Xp]
[~, nobs, ntrials]= size(Yt_event);
nvar = size(Yt_stats.OLS.At,2);
ref_time = CausalParams.ref_time;

% Yt = Yt - mean_Yt;
% diag_flag=0;
for t = 1:nobs    
        
    Xt = squeeze(Yt_event(1:2,t,:));
    Xp = squeeze(Yt_event(3:end,t,:));
    
    switch CausalParams.estim_mode
        case 'OLS'
            coeff = squeeze(Yt_stats.OLS.At(t,:,:));
            SIG = squeeze(Yt_stats.OLS.Sigma_Et(t,:,:));
        case 'RLS'
            coeff = squeeze(Yt_stats.RLS.At(t,:,:));
            SIG = squeeze(Yt_stats.RLS.Sigma_Et(t,:,:));
    end

    A_square = reshape(coeff, nvar, nvar, CausalParams.morder);

    %%
    a = squeeze(A_square(1,1,:));
    b = squeeze(A_square(1,2,:));
    sigy = SIG(1,1);
    
    c = squeeze(A_square(2,1,:));
    d = squeeze(A_square(2,2,:));
    sigx = SIG(2,2);

    %% from this line on X/Y means two variables (Y is the first)
  
    cov_Xp = squeeze(Yt_stats.Sigma(t,4:2:end,4:2:end));
    cov_Yp = squeeze(Yt_stats.Sigma(t,3:2:end-1,3:2:end-1));
    C_XYp = squeeze(Yt_stats.Sigma(t,4:2:end,3:2:end-1));
    C_YXp= squeeze(Yt_stats.Sigma(t,3:2:end-1,4:2:end));   
   
    mean_Xp = Yt_stats.mean(4:2:end,t);
    mean_Yp = Yt_stats.mean(3:2:end-1,t);
    
    CausalOutput.TE(t,2) = .5*log((sigy + b'*cov_Xp*b - b'*C_XYp*inv(cov_Yp)*C_XYp'*b)/sigy );
    CausalOutput.TE(t,1) = .5*log((sigx + c'*cov_Yp*c - c'*C_YXp*inv(cov_Xp)*C_YXp'*c)/sigx );
    
    mean_X_ref =  mean(Yt_stats.mean(4:2:end,ref_time),2);
    mean_Y_ref =  mean(Yt_stats.mean(3:2:end-1,ref_time),2);
    cov_Xp_ref = cov_Xp+mean_Xp*mean_Xp'-mean_Xp*mean_X_ref'-mean_X_ref*mean_Xp'+mean_X_ref*mean_X_ref';
    cov_Yp_ref = cov_Yp+mean_Yp*mean_Yp'-mean_Yp*mean_Y_ref'-mean_Y_ref*mean_Yp'+mean_Y_ref*mean_Y_ref';
    
%      (Xp-mean(mean_Yt(3:end,1:ref_time),2))*(Xp-mean(mean_Yt(3:end,1:ref_time),2))'/ntrials;
%      cov_Xp_lag = (Xp-mean(mean_Yt(3:end,1:ref_time),2))*(Xp-mean(mean_Yt(3:end,1:ref_time),2))'/ntrials;
     cov_Xp_lag = (Xp-mean(Yt_stats.mean(3:end,ref_time),2))*(Xp-mean(Yt_stats.mean(3:end,ref_time),2))'/ntrials;
       
    if ~CausalParams.diag_flag
        CausalOutput.DCS(t,2) = .5*log( (sigy + b'*cov_Xp*b)/sigy );
        CausalOutput.DCS(t,1) = .5*log( (sigx + c'*cov_Yp*c)/sigx );
        if CausalParams.old_version 
        CausalOutput.rDCS(t,2) = .5*log( (sigy + b'*squeeze(mean(Yt_stats.Sigma(ref_time,4:2:end,4:2:end)))*b)/sigy ) - 0.5 + 0.5*(sigy + b'*cov_Xp_lag(2:2:end,2:2:end)*b)/(sigy + b'*squeeze(mean(Yt_stats.Sigma(ref_time,4:2:end,4:2:end)))*b);
        CausalOutput.rDCS(t,1) = .5*log( (sigx + c'*squeeze(mean(Yt_stats.Sigma(ref_time,3:2:end-1,3:2:end-1)))*c)/sigx ) - 0.5 + 0.5*(sigx + c'*cov_Xp_lag(1:2:end-1,1:2:end-1)*c)/(sigx + c'*squeeze(mean(Yt_stats.Sigma(ref_time,3:2:end-1,3:2:end-1)))*c);
        else
        CausalOutput.rDCS(t,2) = .5*log( (sigy + b'*squeeze(mean(Yt_stats.Sigma(ref_time,4:2:end,4:2:end),1))*b)/sigy ) - 0.5 + 0.5*(sigy + b'*cov_Xp_ref*b)/(sigy + b'*squeeze(mean(Yt_stats.Sigma(ref_time,4:2:end,4:2:end)))*b);
        CausalOutput.rDCS(t,1) = .5*log( (sigx + c'*squeeze(mean(Yt_stats.Sigma(ref_time,3:2:end-1,3:2:end-1),1))*c)/sigx ) - 0.5 + 0.5*(sigx + c'*cov_Yp_ref*c)/(sigx + c'*squeeze(mean(Yt_stats.Sigma(ref_time,3:2:end-1,3:2:end-1)))*c);
        end
    end

    if CausalParams.diag_flag
        CausalOutput.DCS(t,2) = .5*log( (sigy + b'*diag(diag(cov_Xp))*b)/sigy );
        CausalOutput.DCS(t,1) = .5*log( (sigx + c'*diag(diag(cov_Yp))*c)/sigx );
        if CausalParams.old_version 
        CausalOutput.rDCS(t,2) = .5*log( (sigy + b'*diag(diag(squeeze(mean(Yt_stats.Sigma(ref_time,4:2:end,4:2:end)))))*b)/sigy ) - 0.5 + 0.5*(sigy + b'*diag(diag(squeeze(cov_Xp_lag(2:2:end,2:2:end))))*b)/(sigy + b'*diag(diag(squeeze(mean(Yt_stats.Sigma(ref_time,4:2:end,4:2:end)))))*b);
        CausalOutput.rDCS(t,1) = .5*log( (sigx + c'*diag(diag(squeeze(mean(Yt_stats.Sigma(ref_time,3:2:end-1,3:2:end-1)))))*c)/sigx ) - 0.5 + 0.5*(sigx + c'*diag(diag(squeeze(cov_Xp_lag(1:2:end-1,1:2:end-1))))*c)/(sigx + c'*diag(diag(squeeze(mean(Yt_stats.Sigma(ref_time,3:2:end-1,3:2:end-1)))))*c);
        else
        CausalOutput.rDCS(t,2) = .5*log( (sigy + b'*diag(diag(squeeze(mean(Yt_stats.Sigma(ref_time,4:2:end,4:2:end)))))*b)/sigy ) - 0.5 + 0.5*(sigy + b'*diag(diag(squeeze(cov_Xp_ref)))*b)/(sigy + b'*diag(diag(squeeze(mean(Yt_stats.Sigma(ref_time,4:2:end,4:2:end)))))*b);
        CausalOutput.rDCS(t,1) = .5*log( (sigx + c'*diag(diag(squeeze(mean(Yt_stats.Sigma(ref_time,3:2:end-1,3:2:end-1)))))*c)/sigx ) - 0.5 + 0.5*(sigx + c'*diag(diag(squeeze(cov_Yp_ref)))*c)/(sigx + c'*diag(diag(squeeze(mean(Yt_stats.Sigma(ref_time,3:2:end-1,3:2:end-1)))))*c);
        end
    end

 end 
 


end
    