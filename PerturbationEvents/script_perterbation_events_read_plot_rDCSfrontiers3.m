% adapted from MVGC toolbox by Kaidi Shao, 10.06.2018, MPI Biological Cybernetics
clear

restoredefaultpath;
% onedrive_path = 'D:\Kaidi\Onedrive\';
onedrive_path = 'C:\Users\skd\OneDrive\';
proj_path = strcat(onedrive_path,'\updated-desnap-with-causality\');
addpath(strcat(onedrive_path,'\util_functions\'));

%% Parameters
set(gcf,'outerposition',[1, 1, 2400, 1200])

cd('C:\Users\skd\OneDrive\updated-desnap-with-causality\Simulations\perturbation_events\saved_data\data_final\pooled_amp5\center\')

N = 100;
CausalStructure.case1.TE = nan(N,200,2);
CausalStructure.case2.TE = nan(N,200,2);
CausalStructure.case3.TE = nan(N,200,2);
CausalStructure.case1.DCS = nan(N,200,2);
CausalStructure.case2.DCS = nan(N,200,2);
CausalStructure.case3.DCS = nan(N,200,2);
CausalStructure.case1.rDCS = nan(N,200,2);
CausalStructure.case2.rDCS = nan(N,200,2);
CausalStructure.case3.rDCS = nan(N,200,2);

for n=[1:N]
load(strcat('perturbation_events_align_0_dim_1_trial_',int2str(n),'_causality.mat'));
CausalStructure.case1.TE(n,:,:) = CausalOutput.OLS.TE;
CausalStructure.case1.DCS(n,:,:) = CausalOutput.OLS.DCS;
CausalStructure.case1.rDCS(n,:,:) = CausalOutput.OLS.rDCS;

load(strcat('perturbation_events_align_1_dim_1_trial_',int2str(n),'_causality.mat'));
CausalStructure.case3.TE(n,:,:) = CausalOutput.OLS.TE;
CausalStructure.case3.DCS(n,:,:) = CausalOutput.OLS.DCS;
CausalStructure.case3.rDCS(n,:,:) = CausalOutput.OLS.rDCS;

load(strcat('perturbation_events_align_1_dim_2_trial_',int2str(n),'_causality.mat'));
CausalStructure.case2.TE(n,:,:) = CausalOutput.OLS.TE;
CausalStructure.case2.DCS (n,:,:)= CausalOutput.OLS.DCS;
CausalStructure.case2.rDCS(n,:,:) = CausalOutput.OLS.rDCS;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

subplot(3,3,1);
mean_temp = squeeze(nanmean(CausalStructure.case1.TE,1));
std_temp = squeeze(nanstd(CausalStructure.case1.TE,[],1));
fill_std_known(mean_temp(:,1),std_temp(:,1),100,-99:100,'b','std',2);
fill_std_known(mean_temp(:,2),std_temp(:,2),100,-99:100,'r','std',2);xlim([-100,100]);ylim([0,20]);ylabel('TE');title('Causality (Ground truth)')

subplot(3,3,2);
mean_temp = squeeze(nanmean(CausalStructure.case3.TE,1));
std_temp = squeeze(nanstd(CausalStructure.case3.TE,[],1));
fill_std_known(mean_temp(:,1),std_temp(:,1),100,-99:100,'b','std',2);
mean_temp = squeeze(nanmean(CausalStructure.case2.TE,1));
std_temp = squeeze(nanstd(CausalStructure.case2.TE,[],1));
fill_std_known(mean_temp(:,2),std_temp(:,2),100,-99:100,'r','std',2);ylim([0,20]);ylabel('TE');title('Causality (single-time)')
subplot(3,3,4);
mean_temp = squeeze(nanmean(CausalStructure.case1.DCS,1));
std_temp = squeeze(nanstd(CausalStructure.case1.DCS,[],1));
fill_std_known(mean_temp(:,1),std_temp(:,1),100,-99:100,'b','std',2);
fill_std_known(mean_temp(:,2),std_temp(:,2),100,-99:100,'r','std',2);xlim([-100,100]);ylim([0,20]);ylabel('DCS');%title('Ground truth)');

subplot(3,3,5);
mean_temp = squeeze(nanmean(CausalStructure.case3.DCS,1));
std_temp = squeeze(nanstd(CausalStructure.case3.DCS,[],1));
fill_std_known(mean_temp(:,1),std_temp(:,1),100,-99:100,'b','std',2);
mean_temp = squeeze(nanmean(CausalStructure.case2.DCS,1));
std_temp = squeeze(nanstd(CausalStructure.case2.DCS,[],1));
fill_std_known(mean_temp(:,2),std_temp(:,2),100,-99:100,'r','std',2);ylim([0,20]);ylabel('DCS');%title('TE (aligned on the putative cause)')

subplot(3,3,7);
mean_temp = squeeze(nanmean(CausalStructure.case1.rDCS,1));
std_temp = squeeze(nanstd(CausalStructure.case1.rDCS,[],1));
fill_std_known(mean_temp(:,1),std_temp(:,1),100,-99:100,'b','std',2);
fill_std_known(mean_temp(:,2),std_temp(:,2),100,-99:100,'r','std',2);xlim([-100,100]);ylim([0,20]);ylabel('rDCS');%title('TE (aligned on the putative cause)')

subplot(3,3,8);
mean_temp = squeeze(nanmean(CausalStructure.case3.rDCS,1));
std_temp = squeeze(nanstd(CausalStructure.case3.rDCS,[],1));
fill_std_known(mean_temp(:,1),std_temp(:,1),100,-99:100,'b','std',2);
mean_temp = squeeze(nanmean(CausalStructure.case2.rDCS,1));
std_temp = squeeze(nanstd(CausalStructure.case2.rDCS,[],1));
fill_std_known(mean_temp(:,2),std_temp(:,2),100,-99:100,'r','std',2);ylim([0,20]);ylabel('rDCS');%title('TE (aligned on the putative cause)')


cd('C:\Users\skd\OneDrive\updated-desnap-with-causality\Simulations\perturbation_events\saved_data\data_final\pooled_amp5\smoothed\')


for n=[1:N]
load(strcat('perturbation_events_align_0_dim_1_trial_',int2str(n),'_causality.mat'));
CausalStructure.case1.TE(n,:,:) = CausalOutput.OLS.TE;
CausalStructure.case1.DCS(n,:,:) = CausalOutput.OLS.DCS;
CausalStructure.case1.rDCS(n,:,:) = CausalOutput.OLS.rDCS;

load(strcat('perturbation_events_align_1_dim_1_trial_',int2str(n),'_causality.mat'));
CausalStructure.case3.TE(n,:,:) = CausalOutput.OLS.TE;
CausalStructure.case3.DCS(n,:,:) = CausalOutput.OLS.DCS;
CausalStructure.case3.rDCS(n,:,:) = CausalOutput.OLS.rDCS;

load(strcat('perturbation_events_align_1_dim_2_trial_',int2str(n),'_causality.mat'));
CausalStructure.case2.TE(n,:,:) = CausalOutput.OLS.TE;
CausalStructure.case2.DCS (n,:,:)= CausalOutput.OLS.DCS;
CausalStructure.case2.rDCS(n,:,:) = CausalOutput.OLS.rDCS;
end


subplot(3,3,3);
mean_temp = squeeze(nanmean(CausalStructure.case3.TE,1));
std_temp = squeeze(nanstd(CausalStructure.case3.TE,[],1));
fill_std_known(mean_temp(:,1),std_temp(:,1),100,-99:100,'b','std',2);
mean_temp = squeeze(nanmean(CausalStructure.case2.TE,1));
std_temp = squeeze(nanstd(CausalStructure.case2.TE,[],1));
fill_std_known(mean_temp(:,2),std_temp(:,2),100,-99:100,'r','std',2);ylim([0,20]);ylabel('TE');title('Causality (smoothed)')


subplot(3,3,6);
mean_temp = squeeze(nanmean(CausalStructure.case3.DCS,1));
std_temp = squeeze(nanstd(CausalStructure.case3.DCS,[],1));
fill_std_known(mean_temp(:,1),std_temp(:,1),100,-99:100,'b','std',2);
mean_temp = squeeze(nanmean(CausalStructure.case2.DCS,1));
std_temp = squeeze(nanstd(CausalStructure.case2.DCS,[],1));
fill_std_known(mean_temp(:,2),std_temp(:,2),100,-99:100,'r','std',2);ylim([0,20]);ylabel('DCS')%title('TE (aligned on the putative cause)')


subplot(3,3,9);
mean_temp = squeeze(nanmean(CausalStructure.case3.rDCS,1));
std_temp = squeeze(nanstd(CausalStructure.case3.rDCS,[],1));
fill_std_known(mean_temp(:,1),std_temp(:,1),100,-99:100,'b','std',2);
mean_temp = squeeze(nanmean(CausalStructure.case2.rDCS,1));
std_temp = squeeze(nanstd(CausalStructure.case2.rDCS,[],1));
fill_std_known(mean_temp(:,2),std_temp(:,2),100,-99:100,'r','std',2);ylim([0,20]);ylabel('rDCS')%title('TE (aligned on the putative cause)')



for i=1:9
    subplot(3,3,i);ylim([0,16]);
    plot([0,0],[0,16],'k--','LineWidth',2);
    ax = gca; ax.FontSize = 20;box off;
    if i>6    xlabel('Peri-event time t^{\prime} (ms)'); end
    if i==7 legend('X^1\rightarrowX^2','X^2\rightarrowX^1');legend boxoff 
    end
    if i==9 legend('X^1\rightarrowX^2|X^1','X^2\rightarrowX^1|X^2');legend boxoff
    end
end