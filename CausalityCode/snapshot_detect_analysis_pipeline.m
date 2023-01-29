function [SnapAnalyOutput, Params, Yt_events] = snapshot_detect_analysis_pipeline(OriSignal, DetSignal, Params)

if Params.Options.Detection == 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SET THRESHOLD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D = DetSignal{1};
if isfield(Params.Detection, 'd0')
    d0 = Params.Detection.d0;
else
    d0 = nanmean(D) + Params.Detection.ThresRatio * nanstd(D);
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FIND REFERENCE POINTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
temp_loc = find(D>=d0);
switch Params.Detection.AlignType
    case "peak"
        locs = find_peak_loc(DetSignal{2}, temp_loc, Params.Detection.L_extract); % align either on detction signal or orginal 
    case "pooled"
        if Params.Detection.ShrinkFlag
            locs = shrink_locs_resample_uniform(temp_loc, ceil(Params.Detection.L_extract/2));
            [locs, distance] = find_best_shrinked_locs(D, locs, temp_loc);
        else
            locs = temp_loc;
        end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PRE-DEFINED REFERENCE POINTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    locs = Params.Detection.locs;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% remove border points %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
locs(locs<Params.Detection.L_extract) = [];
locs(locs>size(OriSignal,2)-Params.Detection.L_extract)= [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BIC model estimation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Params.Options.BIC
disp('performing BIC model selection');
BICParser.OriSignal = OriSignal;
BICParser.DetSignal = DetSignal;
BICParser.Params = Params;
BICParser.EstimMode ='OLS';
switch Params.BIC.mode
    case 'biased'
        Yt_events_momax = get_Yt(OriSignal, locs, Params.BIC.momax, Params.BIC.tau, Params.Detection.L_start, Params.Detection.L_extract); % Xt in snapshots
        BICoutputs = multi_trial_BIC(Yt_events_momax, BICParser);  % for empirical data
        morder = BICoutputs.mobic(2);
end
save(strcat(Params.Output.FileKeyword,'_BIC.mat'), 'Params', 'BICoutputs','-v7.3');
else
    morder = Params.BIC.morder; % preassigned model order
    BICoutputs = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% extract event snapshots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Yt_events = get_Yt(OriSignal, locs, morder, Params.BIC.tau, Params.Detection.L_start, Params.Detection.L_extract); 
if Params.Detection.remove_artif
    [Yt_events, locs] = remove_artif_trials(Yt_events, locs, -15000);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% calculate d-dependent snapshot statistics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Yt_stats = get_Yt_stats(Yt_events, morder);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% calculate d-dependent snapshot statistics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Params.Options.CausalAnlysis
CausalParams = Params.CausalParams;
CausalParams.morder = morder;
[CausalOutput.OLS] = time_varying_causality(Yt_events, Yt_stats, CausalParams);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% perform bootstrapping %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Params.Options.Bootstrap
    disp('Start Bootstrapping!')
    Params.MonteC_Params.morder = morder;
%     Params.MonteC_Params.Ntrials = length(locs);
    [Et] = get_residuals(Yt_events, Yt_stats); % calculate residuals to resample
for n_btsp = 1:Params.MonteC_Params.Nbtsp
    display(strcat('calculating bootstrap trial: ',int2str(n_btsp)));
    Yt_events_btsp = simul_AR_event_btsp(Params.MonteC_Params, Yt_events, Yt_stats, Et); % simulate resampled snapshots
    Yt_stats_btsp = get_Yt_stats(Yt_events_btsp, morder); % gets resampled snapshot statistics
    [CausalOutput_btsp.OLS] = time_varying_causality(Yt_events_btsp, Yt_stats_btsp, CausalParams);

    save(strcat(Params.Output.FileKeyword,'_btsp_',int2str(n_btsp),'_model_causality.mat'), 'Params', 'CausalOutput_btsp', 'Yt_stats_btsp','-v7.3')
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% calculate power spectral density %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
if Params.Options.PSD
    if ~Params.PSD.MonteC_flag
        Yt_stats = get_simul_timefreq(Yt_events, Yt_stats, Params.PSD);
    else
        Params.PSD.simobj.morder = morder;
        [Yt_events_mc] = simul_AR_event(Params.PSD.simobj, Yt_stats);
        Yt_stats = get_simul_timefreq(Yt_events_mc, Yt_stats, Params.PSD);
    end 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% save %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Params.Options.Detection
    SnapAnalyOutput.d0 = d0;
end
SnapAnalyOutput.locs = locs;
SnapAnalyOutput.morder = morder;
SnapAnalyOutput.Yt_stats = Yt_stats;
if Params.Options.CausalAnlysis
    SnapAnalyOutput.CausalOutput = CausalOutput;
end
if Params.Options.BIC
    SnapAnalyOutput.BICoutputs =BICoutputs;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Params.Options.save_flag
    Yt_stats.mean(3:end,:)=[];
    Yt_stats.Sigma(:,3:end,:)=[];
    Yt_stats.Sigma(:,:,3:end)=[];
    Params.DeSnap_inputs.x = [];
    Params.DeSnap_inputs.D = [];
    Params.DeSnap_inputs.Yt_stats_cond = [];
    
    if Params.Options.PSD
        PSD = Yt_stats.spectr;
        save(strcat(Params.Output.FileKeyword,'_psd.mat'), 'Params', 'PSD', '-v7.3');
    else
        if Params.Options.CausalAnlysis
            save(strcat(Params.Output.FileKeyword,'_model_causality.mat'), 'Params', 'Yt_stats', 'CausalOutput', '-v7.3')
        else
            save(strcat(Params.Output.FileKeyword,'_model.mat'), 'Params', 'Yt_stats', '-v7.3');
        end
    end
end

