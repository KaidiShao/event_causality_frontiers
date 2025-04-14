% adapted from MVGC toolbox by Kaidi Shao, 10.06.2018, MPI Biological Cybernetics
% long signal instead of ensemble simulation
clear;

signals = randn(2,10000);
winsize = 50;

momax = 10;
Yt_events_momax = get_Yt_window(signals, winsize, momax);
BICParser.Params.BIC.momax = momax;
BICParser.Params.BIC.mode = 'biased';
BICParser.EstimMode = 'OLS';
BICoutputs = multi_trial_BIC(Yt_events_momax, BICParser);  % for empirical data
morder = BICoutputs.mobic(2);
% morder = 5; % or define manually

Yt_window = get_Yt_window(signals, winsize, morder);
Yt_stats = get_Yt_stats(Yt_window, morder);

CausalParams.old_version = 0; % old version or new version
CausalParams.diag_flag = 0;
CausalParams.ref_time = 1:30;
CausalParams.estim_mode = 'OLS';
CausalParams.morder = morder;
[CausalOutput.OLS] = time_varying_causality(Yt_window, Yt_stats, CausalParams);

