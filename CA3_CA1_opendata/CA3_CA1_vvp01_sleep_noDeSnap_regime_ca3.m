clear;
close all;

restoredefaultpath;
% onedrive_path = 'D:\Kaidi\Onedrive\';%%%%%%%% MPI KYB DESKTOP %%%%%%%%
onedrive_path = 'C:\Users\skd\OneDrive\'; %%%%%%%% DELL LAPTOP %%%%%%%%
proj_path = strcat(onedrive_path,'\updated-desnap-with-causality\');
addpath(genpath(proj_path));
addpath(strcat(onedrive_path,'\util_functions\'));
addpath(genpath(strcat(onedrive_path,'\Toolbox\eeglab10_2_5_8b\')));%%%%%%%%% FOR PSD ANALYSIS %%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Data & channel information %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_path = strcat(proj_path,'opendata_hc3\matlab_data\');
sess_name = 'vvp01_2006-4-9_18-43-47';
load(strcat(sess_name, '_CA3.mat')); 
load(strcat(sess_name, '_CA1.mat')); 
ca3_ch = [1:8]; % HERE FOR PGO IS PBN
ca1_ch = [1:32]; % HERE FOR PGO IS LGN

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DEFINE TASKS TO BE ANALYZED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Params.Options.Detection = 1;
Params.Options.DebiasedStats = 0;
Params.Options.BIC = 1; % 0 for predefined orders, 1 for BIC for given snapshot
Params.Options.PSD = 0;
Params.Options.CausalAnlysis = 1;
Params.Options.Bootstrap = 0;
Params.Options.save_flag = 1;
Params.Output.save_path = strcat(proj_path,'\opendata_hc3\saved_data\bic_test\');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LOAD DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_chpair = (ca3_ch(1)-1)*32+1; % channel pair counter 
L_signal = size(CA1_lfp,1);

for iii = 1:length(ca3_ch)
for jjj = 1:length(ca1_ch)

    display(strcat('calculating: ',int2str(n_chpair)));
    y(1,:) = CA3_lfp(1:L_signal,ca3_ch(iii))';
    y(2,:) = CA1_lfp(1:L_signal,ca1_ch(jjj))';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OBTAIN DETECTION SIGNAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Params.Fs = 1252;
Params.passband = [140,230]; % according to Mizuseki 2009

B = fir1(49, Params.passband/(.5*Params.Fs), 'bandpass');
yf = filter(B,1,y')'; % forward filtering
yf(:, 1:24) = []; % remove FIR-induced delays 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Snapshot Analysis Pipeline %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Params.Detection.L_start = 200;
Params.Detection.L_extract = 401;
Params.Detection.ThresRatio = 5;
Params.Detection.AlignType = 'peak';
Params.Detection.locs = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Params.BIC.momax = 10;
Params.BIC.tau = 1;
Params.BIC.mode = 'biased';
Params.BIC.morder = 7;
Params.Detection.ShrinkFlag = 0;
Params.Detection.remove_artif = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Params.CausalParams.old_version = 0; % old version or new version
Params.CausalParams.diag_flag = 0;
Params.CausalParams.ref_time = 1:200;
Params.CausalParams.estim_mode = 'OLS';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Params.MonteC_Params.L = Params.Detection.L_extract;
Params.MonteC_Params.nvar = 2;
Params.MonteC_Params.Ntrials = 1000;
Params.MonteC_Params.Nbtsp = 100;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Params.DeSnap_inputs.maxStdRatio = 7;
Params.DeSnap_inputs.N_d = 50;
Params.DeSnap_inputs.L_start = Params.Detection.L_start;
Params.DeSnap_inputs.L_extract = Params.Detection.L_extract;
Params.DeSnap_inputs.diff_flag = 1;
Params.DeSnap_inputs.tau = Params.BIC.tau;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% align on sr %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Params.Output.FileKeyword = strcat(Params.Output.save_path, 'ca3_ca1_',sess_name,'_chpair_no_', int2str(n_chpair), '_ca3_', Params.Detection.AlignType);
Params.Detection.ThresRatio = 5;
temp_idx = find(yf(2,:)>=mean(yf(2,:))+Params.Detection.ThresRatio*std(yf(2,:)));
temp_idx = rm_bdry_locs(yf(2,:), temp_idx, Params.Detection.L_start, Params.Detection.L_extract)
yf_labels = zeros(1,length(yf(1,:)));
for i = 1:length(temp_idx)
    yf_labels(temp_idx(i)-Params.Detection.L_start+1:temp_idx(i)-Params.Detection.L_start+Params.Detection.L_extract) = ones(1, Params.Detection.L_extract);
end
yf2 = yf(1,:).*yf_labels;
[SnapAnalyOutput_ca1] = snapshot_detect_analysis_pipeline(y, {yf2; yf2}, Params);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% align on pl %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Params.Output.FileKeyword = strcat(Params.Output.save_path, 'ca3_ca1_',sess_name,'_chpair_no_', int2str(n_chpair), '_ca1_', Params.Detection.AlignType);
Params.Detection.ThresRatio = 5;
temp_idx = find(yf(2,:)>=mean(yf(2,:))+Params.Detection.ThresRatio*std(yf(2,:)));
temp_idx = rm_bdry_locs(yf(2,:), temp_idx, Params.Detection.L_start, Params.Detection.L_extract)
yf_labels = zeros(1,length(yf(2,:)));
for i = 1:length(temp_idx)
    yf_labels(temp_idx(i)-Params.Detection.L_start+1:temp_idx(i)-Params.Detection.L_start+Params.Detection.L_extract) = ones(1, Params.Detection.L_extract);
end
yf2 = yf(2,:).*yf_labels;
[SnapAnalyOutput_ca1] = snapshot_detect_analysis_pipeline(y, {yf2; yf2}, Params);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% align on sr %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Params.Output.FileKeyword = strcat(Params.Output.save_path, 'ca3_ca1_',sess_name,'_chpair_no_', int2str(n_chpair), '_ca3_', Params.Detection.AlignType);
% Params.Detection.ThresRatio = 5;
% [SnapAnalyOutput_ca3] = snapshot_detect_analysis_pipeline(y, {yf(1,:); yf(1,:)}, Params);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% align on pl %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Params.Output.FileKeyword = strcat(Params.Output.save_path, 'ca3_ca1_',sess_name,'_chpair_no_', int2str(n_chpair), '_ca1_', Params.Detection.AlignType);
% Params.Detection.ThresRatio = 5;
% [SnapAnalyOutput_ca1] = snapshot_detect_analysis_pipeline(y, {yf(2,:); yf(2,:)}, Params);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% align on both %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Params.Output.FileKeyword = strcat(Params.Output.save_path, 'ca3_ca1_',sess_name,'_chpair_no_', int2str(n_chpair), '_ca1_', Params.Detection.AlignType);
% Params.Detection.ThresRatio = 5;
% [SnapAnalyOutput_ca3] = snapshot_detect_analysis_pipeline(y, {yf(2,:); yf(2,:)}, Params);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end current loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_chpair=n_chpair+1;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end reminder %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mailTome('finished computation',Params.Output.FileKeyword)
load handel
sound(y,Fs)



