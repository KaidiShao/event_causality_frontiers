% adapted from MVGC toolbox by Kaidi Shao, 10.06.2018, MPI Biological Cybernetics
clear

% specify folder path for timeVaryingCausality/Code
% restoredefaultpath;
% onedrive_path = 'D:\Kaidi\Onedrive\';
onedrive_path = 'C:\Users\skd\OneDrive\';
proj_path = strcat(onedrive_path,'\updated-desnap-with-causality\');
addpath(strcat(onedrive_path,'\util_functions\'));

%% load data
data_path = strcat(proj_path,'\Simulations\perturbation_events\saved_data\data_final\pooled_amp5\center\');
cd(data_path);

% load(strcat('perturbation_events_align_1_dim_2_trial_1_BIC.mat'));
load('perturbation_events_align_0_dim_1_trial_1_model.mat'); 
Yt_stats{1} = Yt_stats_cond; d0all(1)=nan;Imp{1} =Imp_event;At{1}=Yt_stats_cond.OLS.At;bt{1}=Yt_stats_cond.OLS.bt;eta{1}=Yt_stats_cond.OLS.Sigma_Et;
load('perturbation_events_align_1_dim_1_trial_1_model.mat'); 
Yt_stats{3} = Yt_stats_cond; d0all(3)=d0;Imp{3} =Imp_event;At{3}=Yt_stats_cond.OLS.At;bt{3}=Yt_stats_cond.OLS.bt;eta{3}=Yt_stats_cond.OLS.Sigma_Et;
load('perturbation_events_align_1_dim_2_trial_1_model.mat'); 
Yt_stats{2} = Yt_stats_cond; d0all(2)=d0;Imp{2} =Imp_event;At{1}=Yt_stats_cond.OLS.At;bt{1}=Yt_stats_cond.OLS.bt;eta{2}=Yt_stats_cond.OLS.Sigma_Et;

data_path = strcat(proj_path,'\Simulations\perturbation_events\saved_data\data_final\pooled_amp5\smoothed\');
cd(data_path);
load('perturbation_events_align_1_dim_1_trial_1_model.mat'); 
Yt_stats{5} = Yt_stats_cond; d0all(5)=d0;Imp{5} =Imp_event;At{5}=Yt_stats_cond.OLS.At;bt{5}=Yt_stats_cond.OLS.bt;eta{5}=Yt_stats_cond.OLS.Sigma_Et;
load('perturbation_events_align_1_dim_2_trial_1_model.mat'); 
Yt_stats{4} = Yt_stats_cond; d0all(4)=d0;Imp{4} =Imp_event;At{4}=Yt_stats_cond.OLS.At;bt{4}=Yt_stats_cond.OLS.bt;eta{4}=Yt_stats_cond.OLS.Sigma_Et;

%%% plot
figure;
set(gcf,'outerposition',[1, 1, 2400, 1200])
t=-99:100;

subplot(3,3,1);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mean_temp = mean(squeeze(Imp{1}(2,:,:)),2);
std_temp = std(squeeze(Imp{1}(2,:,:)),[],2);
fill_std_known(mean_temp, std_temp, nan, t, 'k', 'std', 2); hold on; 
plot(zeros(1,2),gca().YLim,'k--','LineWidth',2);legend('k^2_t'); legend boxoff;title('Ground truth');ylabel({'Hidden';'state'});ylim([-5,5.5])

subplot(3,3,2);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mean_temp = mean(squeeze(Imp{3}(2,:,:)),2);
std_temp = std(squeeze(Imp{3}(2,:,:)),[],2);
fill_std_known(mean_temp, std_temp, nan, t, 'b', 'std', 2); hold on;

mean_temp = mean(squeeze(Imp{2}(2,:,:)),2);
std_temp = std(squeeze(Imp{2}(2,:,:)),[],2);
fill_std_known(mean_temp, std_temp, nan, t, 'r', 'std', 2); hold on; 
plot(zeros(1,2),gca().YLim,'k--','LineWidth',2);%legend('k^2_t|X^1','k^2_t|X^2'); legend boxoff;
title('Single-time selection');ylim([-5,5.5])

subplot(3,3,3);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mean_temp = mean(squeeze(Imp{5}(2,:,:)),2);
std_temp = std(squeeze(Imp{5}(2,:,:)),[],2);
fill_std_known(mean_temp, std_temp, nan, t, 'b', 'std', 2); hold on;

mean_temp = mean(squeeze(Imp{4}(2,:,:)),2);
std_temp = std(squeeze(Imp{4}(2,:,:)),[],2);
fill_std_known(mean_temp, std_temp, nan, t, 'r', 'std', 2); hold on;
plot(zeros(1,2),gca().YLim,'k--','LineWidth',2);legend('k^2_t|X^1','k^2_t|X^2'); legend boxoff;title('Smoothed selection');ylim([-5,5.5])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(3,3,4);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mean_temp = Yt_stats{1}.mean(1,:)';
std_temp = sqrt(squeeze(Yt_stats{1}.Sigma(:,1,1)));
fill_std_known(mean_temp, std_temp, nan, t, 'b', 'std', 2);hold on;ylim([-20,21]);
plot(zeros(1,2),gca().YLim,'k--','LineWidth',2);legend('X^1_t'); legend boxoff;ylabel({'Event'; 'waveforms'})

subplot(3,3,5);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(t(1,[1,end]),d0all(3)*ones(1,2),'b','LineWidth',1);hold on;

mean_temp = Yt_stats{3}.mean(1,:)';
std_temp = sqrt(squeeze(Yt_stats{3}.Sigma(:,1,1)));
fill_std_known(mean_temp, std_temp, nan, t, 'b', 'std', 2);hold on

mean_temp = Yt_stats{3}.mean(2,:)';
std_temp = sqrt(squeeze(Yt_stats{3}.Sigma(:,2,2)));
fill_std_known(mean_temp, std_temp, nan, t, 'r', 'std', 2);hold on;ylim([-20,21]);
plot(zeros(1,2),gca().YLim,'k--','LineWidth',2); legend('threshold d_0'); legend boxoff;


subplot(3,3,6);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mean_temp = Yt_stats{5}.mean(1,:)';
std_temp = sqrt(squeeze(Yt_stats{5}.Sigma(:,1,1)));
fill_std_known(mean_temp, std_temp, nan, t, 'b', 'std', 2);hold on

mean_temp = Yt_stats{5}.mean(2,:)';
std_temp = sqrt(squeeze(Yt_stats{5}.Sigma(:,2,2)));
fill_std_known(mean_temp, std_temp, nan, t, 'r', 'std', 2);hold on;ylim([-20,21])
plot(t(1,[1,end]),d0all(5)*ones(1,2),'b','LineWidth',1);
plot(zeros(1,2),gca().YLim,'k--','LineWidth',2);legend('X^1_t|X^1','X^2_t|X^1'); legend boxoff;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(3,3,7);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mean_temp = Yt_stats{1}.mean(2,:)';
std_temp = sqrt(squeeze(Yt_stats{1}.Sigma(:,2,2)));
fill_std_known(mean_temp, std_temp, nan, t, 'r', 'std', 2);hold on;ylabel('X^2_t');ylim([-16,22]);
plot(zeros(1,2),gca().YLim,'k--','LineWidth',2);legend('X^2_t'); legend boxoff;ylabel({'Event'; 'waveforms'})

subplot(3,3,8);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(t(1,[1,end]),d0all(2)*ones(1,2),'r','LineWidth',1);hold on
mean_temp = Yt_stats{2}.mean(1,:)';
std_temp = sqrt(squeeze(Yt_stats{2}.Sigma(:,1,1)));
fill_std_known(mean_temp, std_temp, nan, t, 'b', 'std', 2);hold on;

mean_temp = Yt_stats{2}.mean(2,:)';
std_temp = sqrt(squeeze(Yt_stats{2}.Sigma(:,2,2)));
fill_std_known(mean_temp, std_temp, nan, t, 'r', 'std', 2);hold on;ylim([-16,22]);
plot(zeros(1,2),gca().YLim,'k--','LineWidth',2); legend('threshold d_0'); legend boxoff;


subplot(3,3,9);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mean_temp = Yt_stats{4}.mean(1,:)';
std_temp = sqrt(squeeze(Yt_stats{4}.Sigma(:,1,1)));
fill_std_known(mean_temp, std_temp, nan, t, 'b', 'std', 2);hold on

mean_temp = Yt_stats{4}.mean(2,:)';
std_temp = sqrt(squeeze(Yt_stats{4}.Sigma(:,2,2)));
fill_std_known(mean_temp, std_temp, nan, t, 'r', 'std', 2);hold on;ylim([-16,22]);
plot(t(1,[1,end]),d0all(4)*ones(1,2),'r','LineWidth',1);
plot(zeros(1,2),gca().YLim,'k--','LineWidth',2);legend('X^1_t|X^2','X^2_t|X^2'); legend boxoff;

for n=1:9
    subplot(3,3,n);
    ax = gca; ax.FontSize = 20;box off;
    if n>6 xlabel('Peri-event time t^{\prime} (ms)'); end
end



