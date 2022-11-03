% adapted from MVGC toolbox by Kaidi Shao, 10.06.2018, MPI Biological Cybernetics
% long signal instead of ensemble simulation

clear

% specify folder path for timeVaryingCausality/Code
restoredefaultpath;
% onedrive_path = 'D:\Kaidi\Onedrive\';
onedrive_path = 'C:\Users\skd\OneDrive\';
proj_path = strcat(onedrive_path,'\updated-desnap-with-causality\');
addpath(genpath(proj_path));
addpath(strcat(onedrive_path,'\util_functions\'));
addpath(genpath(strcat(onedrive_path,'\Toolbox\eeglab10_2_5_8b\')));
%% Parameters

nvar       = 2;       % number of variables
ntrials     = 50000;    % number of trials
L_perturb   = 1;

morder      = 1;       % number of lags in VAR model
fs          = 1000;    % sampling frequency

align_flag_range       = [0 1 1 ]; % 1 for aligning; 0 for not aligning (aligning on innovations)
dim_align_range        = [1 1 2 ]; % 1 for alignment on effect; 2 for alignment on cause

momax    = 10;
BIC_flag = 0;


%% Construct a long stationary VAR time series with event number as ntrials
amp = 10; % amplitude of the input events through innovations
dim = 2; % dimension of variables that the input goes into

a = [0.5];
b = [0.5];
c = [2];

SIG_origin = [.1 0;0 .1]; % constant "minimal VAR" residuals covariance

for d0_scale = [-4:1:4]
%%
A = [a(1),c(1);
           0, b(1)];
SIG = SIG_origin;

disp('Start generating VAR process!');
L_event_gen = 500;
center_perturb = 350;
[X, Imp] = gen_ensemble_nonstat_innomean(A, SIG, 100, L_event_gen, center_perturb, amp, dim, L_perturb);

X_template = mean(X(:,center_perturb-50:center_perturb+50,:),3);
Imp_grdt = mean(Imp(:,center_perturb-50:center_perturb+50,:),3);

%%

for nRep = 1:2 % no. of repetitions
    clear x xf
temp_t_diff = 100 + randi(200,1,ntrials);
loc = 1000+cumsum(temp_t_diff);
L = loc(end)+1000; 

Imp_long = zeros(nvar,L);
% [x_stat] = simul_AR_kaidi_nonstat_innomean(A, SIG, Imp_long, morder);

for n = 1:length(loc)
    Imp_long(:,loc(n)-50:loc(n)+50) = Imp_grdt;
end

[x] = simul_AR_kaidi_nonstat_innomean(A, SIG, Imp_long, morder);
Imp_long(:,1:morder)=[];

%%
% xf(1,:) = conv(x(1,:),fliplr(X_template(1,:)));
% xf(2,:) = conv(x(2,:),fliplr(X_template(2,:)));
% xf(:,1:52)=[];
xf = x;

for i = 1:length(align_flag_range) % no. of cases


[i nRep]
align_flag = align_flag_range(i); % 1 for alignment; 0 for non-alignment
dim_align = dim_align_range(i);        % 1 is the effect; 2 is the cause
L_start = 100;
L_extract = 200;
momax = 10;
tau = 1;

%% preprocessing --> alignment and jittering
disp('Start preprocessing --> alignment & jittering!');

% alignment and jitter choices
if ~align_flag  % if we don't align the trials
    Yt_events_momax = get_Yt(x, loc, momax, tau, L_start, L_extract);
    Imp_event = get_Yt(Imp_long, loc, 0, tau, L_start, L_extract);
else
     D = x(dim_align,:); 
%      loc=loc-10;
     d0 = mean(D(loc))+1*std(D(loc));
     idx = find(D(loc)>d0); loc_peak = loc(idx);
%      loc_peak = temp_loc;
%      loc_peak = find_peak_loc(D, temp_loc, 50);
%     idx = find(D(loc)>d0); loc_peak = loc(idx);
%      d0 = mean(D)+3*std(D);
%     temp_loc = find(D>=d0);
%     loc_peak = temp_loc(1:end-100);
%     loc_peak = find_peak_loc(D, temp_loc, 50);
%     idx = find(D(loc)>d0); loc_peak = loc(idx);
    Yt_events_momax = get_Yt(x, loc_peak, momax, tau, L_start, L_extract);
    Imp_event = get_Yt(Imp_long, loc_peak, 0, tau, L_start, L_extract);
end


%% BIC model order selection
if BIC_flag
    disp('Start model order selection procedure!');
    BICoutputs = multi_trial_BIC(Yt_events_momax, momax,'OLS');  % for empirical data
    mobic = BICoutputs.mobic(2);
end
mobic = 1; % ground-truth model order

%% 
Yt_events = Yt_events_momax(1:nvar*(mobic+1),:,:);
Yt_stats_cond = get_Yt_stats_cond(Yt_events, mobic);



%% calculation of causality measures
CausalParams.old_version = 0; % old version or new version
CausalParams.diag_flag = 0;
CausalParams.ref_time = 1:30;
CausalParams.estim_mode = 'OLS';
CausalParams.morder = mobic;
[CausalOutput.OLS] = time_varying_causality(Yt_events, Yt_stats_cond, CausalParams);

% mkdir(strcat('C:\Users\skd\OneDrive\updated-desnap-with-causality\Simulations\perturbation_events\test_inequality\d0_',num2str(d0_scale),'\'));
% cd(strcat('C:\Users\skd\OneDrive\updated-desnap-with-causality\Simulations\perturbation_events\test_inequality\d0_',num2str(d0_scale),'\'));
% % save(strcat('perturbation_events_align_',int2str(align_flag_range(i)),'_dim_',int2str(dim_align),'_trial_',int2str(n),'_BIC.mat'),'BICoutputs');
% save(strcat('perturbation_events_align_',int2str(align_flag_range(i)),'_dim_',int2str(dim_align),'_trial_',int2str(nRep),'_model.mat'),'Yt_stats_cond','Imp_event','A','SIG');
% save(strcat('perturbation_events_align_',int2str(align_flag_range(i)),'_dim_',int2str(dim_align),'_trial_',int2str(nRep),'_causality.mat'),'CausalParams','CausalOutput');

if ~align_flag  % if we don't align the trials
    Yt_events_grdt = Yt_events;
    Yt_stats_cond_grdt = Yt_stats_cond;
    CausalOutput_grdt = CausalOutput.OLS;
else
    if dim_align == 1
        Yt_events_peak1 = Yt_events;
        Yt_stats_cond_peak1 = Yt_stats_cond;
            CausalOutput_peak1 = CausalOutput.OLS;
    end
    if dim_align == 2
        Yt_events_peak2 = Yt_events;
        Yt_stats_cond_peak2 = Yt_stats_cond;
            CausalOutput_peak2 = CausalOutput.OLS;
    end
end

end
end
end

