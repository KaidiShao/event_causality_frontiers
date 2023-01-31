% adapted from MVGC toolbox by Kaidi Shao, 10.06.2018, MPI Biological Cybernetics
% long signal instead of ensemble simulation
clear;

path = '...\event_causality_frontiers\'; % change according to individual needs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parameter settings for simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Params.VAR.nvar = 2;         % number of variables
Params.VAR.ntrials = 5000;   % number of trials
Params.VAR.L_perturb = 101;  % length of deterministic perturbation
Params.VAR.morder = 4;       % number of lags in VAR model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Params.VAR.amp = 5; % amplitude of the input events through innovations
Params.VAR.dim = 2;  % dimension of variables that the input goes into
align_flag_range = [0 1 1 2 2]; % 2 for smoothed alignment, 1 for single-time alignment; 0 for not alignment (aligning on innovations)
dim_align_range = [1 1 2 1 2];  % 1 for alignment on effect; 2 for alignment on cause
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a = [-0.55 -0.45 -0.55 -0.85];
b = [0.9 -0.25 0.01 0.25];
c = [1.4 -0.3 1.5 1.7];
Params.VAR.A = [a(1),c(1),a(2),c(2),a(3),c(3),a(4),c(4);
     0, b(1), 0, b(2), 0, b(3), 0, b(4)]; % VAR coefficients
Params.VAR.SIG = [1 0;0 1]; % constant "minimal VAR" residuals covariance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NRep = 20; % no. of repetitions
rng_range = 1:NRep; % random seed to control simulation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parameter settings for analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DetectSigType = 'FIR'; % 'FIR', 'template', 'original'
Fs = 1000;
f_max = 60; % peak of spectrum
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DEFINE TASKS TO BE ANALYZED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Params.Options.generate_flag = 1;
Params.Options.Detection = 1;
Params.Options.BIC = 0; % 0 for predefined orders, 1 for BIC for given snapshot
Params.Options.PSD = 0;
Params.Options.CausalAnlysis = 1;
Params.Options.Bootstrap = 0;
Params.Options.save_flag = 1;
Params.Output.save_path = path; % change according to individual needs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Snapshot Analysis Pipeline %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Params.Detection.L_start = 100;
Params.Detection.L_extract = 200;
Params.Detection.ThresRatio = 3;
Params.Detection.AlignType = 'peak';
Params.Detection.locs = [];
Params.Detection.ShrinkFlag = 1;
Params.Detection.remove_artif = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Params.BIC.momax = 15;
Params.BIC.tau = 1;
Params.BIC.mode = 'biased'; 
Params.BIC.morder = 4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Params.CausalParams.old_version = 0; % old version or new version
Params.CausalParams.diag_flag = 0;
Params.CausalParams.ref_time = 1:30;
Params.CausalParams.estim_mode = 'OLS';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Params.MonteC_Params.L = Params.Detection.L_extract;
Params.MonteC_Params.nvar = 2;
Params.MonteC_Params.Ntrials = 1000;
Params.MonteC_Params.Nbtsp = 100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Simulation of templates %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Start generating event ensembles for templates!');
L_event_gen = 500;
center_perturb = 350;
[X, Imp] = gen_ensemble_nonstat_innomean(Params.VAR.A, Params.VAR.SIG, Params.VAR.ntrials, L_event_gen, center_perturb, Params.VAR.amp, Params.VAR.dim, Params.VAR.L_perturb);

X_template = mean(X(:,center_perturb-50:center_perturb+50,:),3);
Imp_grdt = mean(Imp(:,center_perturb-50:center_perturb+50,:),3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Start Repetitive Simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for nRep = 1:NRep
    rng(rng_range(nRep))

    clear x xf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Define where events occur %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
temp_t_diff = 400 + randi(400,1,Params.VAR.ntrials);
locs = cumsum(temp_t_diff); % ground truth peak in hidden states
L = locs(end)+1000; 

Imp_long = zeros(Params.VAR.nvar,L);
for n = 1:length(locs)
    Imp_long(:,locs(n)-50:locs(n)+50) = Imp_grdt;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% long signal with events %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[x] = simul_AR_kaidi_nonstat_innomean(Params.VAR.A, Params.VAR.SIG, Imp_long, Params.VAR.morder);
Imp_long(:,1:Params.VAR.morder)=[];

plot(x(1,:))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OBTAIN DETECTION SIGNAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(DetectSigType, 'FIR') 
B = fir1(49, [f_max-10,f_max+10]/(.5*Fs), 'bandpass'); % define FIR dilter
xf = filter(B,1,x')'; % forward filtering
xf(:,1:24)=[];
end
if strcmp(DetectSigType, 'template') 
xf(1,:) = conv(x(1,:),fliplr(X_template(1,:)));
xf(2,:) = conv(x(2,:),fliplr(X_template(2,:)));
xf(:,1:100)=[];
end
if strcmp(DetectSigType, 'original') 
xf = x;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% loop for different alignments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:length(align_flag_range) % no. of cases
disp(strcat('Start analysis for alignment method:', int2str(i), ' in Repetition No.', int2str(nRep)));
align_flag = align_flag_range(i); % 1 for alignment; 0 for non-alignment
dim_align = dim_align_range(i);        % 1 is the effect; 2 is the cause
Params.Output.FileKeyword = strcat(Params.Output.save_path, 'perturbation_events_align_',int2str(align_flag_range(i)),'_dim_',int2str(dim_align),'_trial_',int2str(nRep));

save('locs.mat','locs')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% process for alignments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch align_flag
    case 0
        Params.Options.Detection = 0;
        Params.Detection.locs = locs;
        [SnapAnalyOutput] = snapshot_detect_analysis_pipeline(x, {xf(1,:); xf(1,:)}, Params);
    case 1
        Params.Options.Detection = 0;
        D = x(dim_align,locs);
        Params.Detection.ThresRatio = 2.1; %
        d0 = mean(D)+Params.Detection.ThresRatio*std(D);
        loc_idx = find(x(dim_align,locs)>d0);
        Params.Detection.locs = locs(loc_idx);   
        [SnapAnalyOutput] = snapshot_detect_analysis_pipeline(x, {xf(1,:); xf(1,:)}, Params);
    case 2
        Params.Options.Detection = 1; 
        Params.Detection.ThresRatio = 4.2;
        [SnapAnalyOutput] = snapshot_detect_analysis_pipeline(x, {xf(dim_align,:); xf(dim_align,:)}, Params);
end

end
end
