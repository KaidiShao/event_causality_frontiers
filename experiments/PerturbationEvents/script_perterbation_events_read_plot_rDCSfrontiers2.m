% adapted from MVGC toolbox by Kaidi Shao, 10.06.2018, MPI Biological Cybernetics
clear

% restoredefaultpath;
% onedrive_path = 'D:\Kaidi\Onedrive\';
onedrive_path = 'C:\Users\skd\OneDrive\';
proj_path = strcat(onedrive_path,'\updated-desnap-with-causality\');
addpath(strcat(onedrive_path,'\util_functions\'));

%% Parameters


cd('C:\Users\skd\OneDrive\updated-desnap-with-causality\Simulations\perturbation_events\saved_data\data_final\pooled_amp5\center\')

N = 100;
At=nan(3,N,200,2);

for n=[1:N]
load(strcat('perturbation_events_align_0_dim_1_trial_',int2str(n),'_model.mat'));
At(1,n,:,1)=Yt_stats_cond.OLS.At(:,1,2);
At(1,n,:,2)=Yt_stats_cond.OLS.At(:,2,7);

load(strcat('perturbation_events_align_1_dim_1_trial_',int2str(n),'_model.mat'));
At(3,n,:,1)=Yt_stats_cond.OLS.At(:,1,2);
At(3,n,:,2)=Yt_stats_cond.OLS.At(:,2,7);

load(strcat('perturbation_events_align_1_dim_2_trial_',int2str(n),'_model.mat'));
At(2,n,:,1)=Yt_stats_cond.OLS.At(:,1,2);
At(2,n,:,2)=Yt_stats_cond.OLS.At(:,2,7);
end


% cd('C:\Users\skd\OneDrive\updated-desnap-with-causality\Simulations\perturbation_events\data_final\pooled_amp5\smoothed\')
cd('C:\Users\skd\OneDrive\updated-desnap-with-causality\Simulations\perturbation_events\saved_data\data_final\pooled_amp5\smoothed\')

N = 100;
At2=nan(3,N,200,2);

for n=[1:N]
load(strcat('perturbation_events_align_0_dim_1_trial_',int2str(n),'_model.mat'));
At2(1,n,:,1)=Yt_stats_cond.OLS.At(:,1,2);
At2(1,n,:,2)=Yt_stats_cond.OLS.At(:,2,7);

load(strcat('perturbation_events_align_1_dim_1_trial_',int2str(n),'_model.mat'));
At2(3,n,:,1)=Yt_stats_cond.OLS.At(:,1,2);
At2(3,n,:,2)=Yt_stats_cond.OLS.At(:,2,7);

load(strcat('perturbation_events_align_1_dim_2_trial_',int2str(n),'_model.mat'));
At2(2,n,:,1)=Yt_stats_cond.OLS.At(:,1,2);
At2(2,n,:,2)=Yt_stats_cond.OLS.At(:,2,7);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
set(gcf,'outerposition',[1, 1, 2400, 400])
subplot(1,3,1);
mean_temp = squeeze(nanmean(At(1,:,:,2),2));
std_temp = squeeze(nanstd(At(1,:,:,2),[],2));
fill_std_known(mean_temp,std_temp,100,-99:100,'b','std',2);
mean_temp = squeeze(nanmean(At(1,:,:,1),2));
std_temp = squeeze(nanstd(At(1,:,:,1),[],2));
fill_std_known(mean_temp,std_temp,100,-99:100,'r','std',2);
title('Coupling strength (Ground truth)')
    plot([0,0],[-.2,1.5],'k--','LineWidth',2);
legend('X^1\rightarrowX^2','X^2\rightarrowX^1');legend boxoff 

subplot(1,3,2);
mean_temp = squeeze(nanmean(At(3,:,:,2),2));
std_temp = squeeze(nanstd(At(3,:,:,2),[],2));
fill_std_known(mean_temp,std_temp,100,-99:100,'b','std',2);
mean_temp = squeeze(nanmean(At(2,:,:,1),2));
std_temp = squeeze(nanstd(At(2,:,:,1),[],2));
fill_std_known(mean_temp,std_temp,100,-99:100,'r','std',2);
xlim([-100,100]);ylim([0,1.5]);title('Coupling strength (single-time)')
plot([0,0],[-.2,1.5],'k--','LineWidth',2);
legend('X^1\rightarrowX^2|X^1','X^2\rightarrowX^1|X^2');legend boxoff 


subplot(1,3,3);
mean_temp = squeeze(nanmean(At2(3,:,:,2),2));
std_temp = squeeze(nanstd(At2(3,:,:,2),[],2));
fill_std_known(mean_temp,std_temp,100,-99:100,'b','std',2);
mean_temp = squeeze(nanmean(At2(2,:,:,1),2));
std_temp = squeeze(nanstd(At2(2,:,:,1),[],2));
fill_std_known(mean_temp,std_temp,100,-99:100,'r','std',2);
xlim([-100,100]);ylim([0,1.5]);title('Coupling strength (local peak)');
plot([0,0],[-.2,1.5],'k--','LineWidth',2);
legend('X^1\rightarrowX^2|X^1','X^2\rightarrowX^1|X^2');legend boxoff 




for i=1:3
    subplot(1,3,i);
    xlim([-100,100]);ylim([-.2,1.5]);
    if i==1 ylabel('Coefficients'); end
    ax = gca; ax.FontSize = 20;box off;
    xlabel('Peri-event time t^{\prime} (ms)');
end