% adapted from MVGC toolbox by Kaidi Shao, 10.06.2018, MPI Biological Cybernetics
clear

restoredefaultpath;
% onedrive_path = 'D:\Kaidi\Onedrive\';
onedrive_path = 'C:\Users\skd\OneDrive\';
proj_path = strcat(onedrive_path,'\updated-desnap-with-causality\');
addpath(strcat(onedrive_path,'\util_functions\'));

%% Parameters

cd('C:\Users\skd\OneDrive\updated-desnap-with-causality\Simulations\perturbation_events\data_final\pooled_amp5\center\')

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
Yt_stats.mean_sr = nan(N,200,2);
Yt_stats.mean_pl = nan(N,200,2);
Yt_stats.Sigma_sr= nan(N,200,2);
Yt_stats.Sigma_sr =nan(N,200,2);
Yt_stats.Sigma_pl = nan(N,200,2);
Yt_stats.Sigma_pl = nan(N,200,2);


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

CausalStructure.case1.TE(isinf(CausalStructure.case1.TE)) = nan;
CausalStructure.case1.DCS(isinf(CausalStructure.case1.DCS)) = nan;
CausalStructure.case1.rDCS(isinf(CausalStructure.case1.rDCS)) = nan;


CausalStructure.case2.TE(isinf(CausalStructure.case2.TE)) = nan;
CausalStructure.case2.DCS(isinf(CausalStructure.case2.DCS)) = nan;
CausalStructure.case2.rDCS(isinf(CausalStructure.case2.rDCS)) = nan;


CausalStructure.case3.TE(isinf(CausalStructure.case3.TE)) = nan;
CausalStructure.case3.DCS(isinf(CausalStructure.case3.DCS)) = nan;
CausalStructure.case3.rDCS(isinf(CausalStructure.case3.rDCS)) = nan;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
figure;

subplot(3,5,1);
mean_temp = squeeze(nanmean(CausalStructure.case1.TE,1));
std_temp = squeeze(nanstd(CausalStructure.case1.TE,[],1));
fill_std_known(mean_temp(:,1),std_temp(:,1),100,-99:100,'b','std',2);
fill_std_known(mean_temp(:,2),std_temp(:,2),100,-99:100,'r','std',2);xlim([-100,100]);ylim([0,20]);title('TE (ground truth)')

subplot(3,5,2);
mean_temp = squeeze(nanmean(CausalStructure.case3.TE,1));
std_temp = squeeze(nanstd(CausalStructure.case3.TE,[],1));
fill_std_known(mean_temp(:,1),std_temp(:,1),100,-99:100,'b','std',2);
mean_temp = squeeze(nanmean(CausalStructure.case2.TE,1));
std_temp = squeeze(nanstd(CausalStructure.case2.TE,[],1));
fill_std_known(mean_temp(:,2),std_temp(:,2),100,-99:100,'r','std',2);ylim([0,20]);title('TE (aligned on the putative cause)')

subplot(3,5,3);
mean_temp = squeeze(nanmean(CausalStructure.case2.TE,1));
std_temp = squeeze(nanstd(CausalStructure.case2.TE,[],1));
fill_std_known(mean_temp(:,1),std_temp(:,1),100,-99:100,'b','std',2);
mean_temp = squeeze(nanmean(CausalStructure.case3.TE,1));
std_temp = squeeze(nanstd(CausalStructure.case3.TE,[],1));
fill_std_known(mean_temp(:,2),std_temp(:,2),100,-99:100,'r','std',2);ylim([0,20]);title('TE (aligned on the putative effect)')

subplot(3,5,4);
mean_temp = squeeze(nanmean(CausalStructure.case3.TE,1));
std_temp = squeeze(nanstd(CausalStructure.case3.TE,[],1));
fill_std_known(mean_temp(:,2),std_temp(:,2),100,-99:100,'b','std',2);
mean_temp = squeeze(nanmean(CausalStructure.case2.TE,1));
std_temp = squeeze(nanstd(CausalStructure.case3.TE,[],1));
fill_std_known(mean_temp(:,2),std_temp(:,2),100,-99:100,'r','std',2);ylim([0,20]);title('TE (ground truth direction)')

subplot(3,5,5);
mean_temp = squeeze(nanmean(CausalStructure.case3.TE,1));
std_temp = squeeze(nanstd(CausalStructure.case3.TE,[],1));
fill_std_known(mean_temp(:,1),std_temp(:,1),100,-99:100,'b','std',2);
mean_temp = squeeze(nanmean(CausalStructure.case2.TE,1));
std_temp = squeeze(nanstd(CausalStructure.case2.TE,[],1));
fill_std_known(mean_temp(:,1),std_temp(:,1),100,-99:100,'r','std',2);ylim([0,20]);title('TE (non ground truth direction)')

subplot(3,5,6);
mean_temp = squeeze(nanmean(CausalStructure.case1.DCS,1));
std_temp = squeeze(nanstd(CausalStructure.case1.DCS,[],1));
fill_std_known(mean_temp(:,1),std_temp(:,1),100,-99:100,'b','std',2);
fill_std_known(mean_temp(:,2),std_temp(:,2),100,-99:100,'r','std',2);xlim([-100,100]);ylim([0,20]);title('TE (ground truth)')

subplot(3,5,7);
mean_temp = squeeze(nanmean(CausalStructure.case3.DCS,1));
std_temp = squeeze(nanstd(CausalStructure.case3.DCS,[],1));
fill_std_known(mean_temp(:,1),std_temp(:,1),100,-99:100,'b','std',2);
mean_temp = squeeze(nanmean(CausalStructure.case2.DCS,1));
std_temp = squeeze(nanstd(CausalStructure.case2.DCS,[],1));
fill_std_known(mean_temp(:,2),std_temp(:,2),100,-99:100,'r','std',2);ylim([0,20]);title('TE (aligned on the putative cause)')

subplot(3,5,8);
mean_temp = squeeze(nanmean(CausalStructure.case2.DCS,1));
std_temp = squeeze(nanstd(CausalStructure.case2.DCS,[],1));
fill_std_known(mean_temp(:,1),std_temp(:,1),100,-99:100,'b','std',2);
mean_temp = squeeze(nanmean(CausalStructure.case3.DCS,1));
std_temp = squeeze(nanstd(CausalStructure.case3.DCS,[],1));
fill_std_known(mean_temp(:,2),std_temp(:,2),100,-99:100,'r','std',2);ylim([0,20]);title('TE (aligned on the putative effect)')

subplot(3,5,9);
mean_temp = squeeze(nanmean(CausalStructure.case3.DCS,1));
std_temp = squeeze(nanstd(CausalStructure.case3.DCS,[],1));
fill_std_known(mean_temp(:,2),std_temp(:,2),100,-99:100,'b','std',2);
mean_temp = squeeze(nanmean(CausalStructure.case2.DCS,1));
std_temp = squeeze(nanstd(CausalStructure.case3.DCS,[],1));
fill_std_known(mean_temp(:,2),std_temp(:,2),100,-99:100,'r','std',2);ylim([0,20]);title('TE (ground truth direction)')

subplot(3,5,10);
mean_temp = squeeze(nanmean(CausalStructure.case3.DCS,1));
std_temp = squeeze(nanstd(CausalStructure.case3.DCS,[],1));
fill_std_known(mean_temp(:,1),std_temp(:,1),100,-99:100,'b','std',2);
mean_temp = squeeze(nanmean(CausalStructure.case2.DCS,1));
std_temp = squeeze(nanstd(CausalStructure.case2.DCS,[],1));
fill_std_known(mean_temp(:,1),std_temp(:,1),100,-99:100,'r','std',2);ylim([0,20]);title('TE (non ground truth direction)')

subplot(3,5,11);
mean_temp = squeeze(nanmean(CausalStructure.case1.rDCS,1));
std_temp = squeeze(nanstd(CausalStructure.case1.rDCS,[],1));
fill_std_known(mean_temp(:,1),std_temp(:,1),100,-99:100,'b','std',2);
fill_std_known(mean_temp(:,2),std_temp(:,2),100,-99:100,'r','std',2);xlim([-100,100]);ylim([0,20]);title('TE (ground truth)')

subplot(3,5,12);
mean_temp = squeeze(nanmean(CausalStructure.case3.rDCS,1));
std_temp = squeeze(nanstd(CausalStructure.case3.rDCS,[],1));
fill_std_known(mean_temp(:,1),std_temp(:,1),100,-99:100,'b','std',2);
mean_temp = squeeze(nanmean(CausalStructure.case2.rDCS,1));
std_temp = squeeze(nanstd(CausalStructure.case2.rDCS,[],1));
fill_std_known(mean_temp(:,2),std_temp(:,2),100,-99:100,'r','std',2);ylim([0,20]);title('TE (aligned on the putative cause)')

subplot(3,5,13);
mean_temp = squeeze(nanmean(CausalStructure.case2.rDCS,1));
std_temp = squeeze(nanstd(CausalStructure.case2.rDCS,[],1));
fill_std_known(mean_temp(:,1),std_temp(:,1),100,-99:100,'b','std',2);
mean_temp = squeeze(nanmean(CausalStructure.case3.rDCS,1));
std_temp = squeeze(nanstd(CausalStructure.case3.rDCS,[],1));
fill_std_known(mean_temp(:,2),std_temp(:,2),100,-99:100,'r','std',2);ylim([0,20]);title('TE (aligned on the putative effect)')

subplot(3,5,14);
mean_temp = squeeze(nanmean(CausalStructure.case3.rDCS,1));
std_temp = squeeze(nanstd(CausalStructure.case3.rDCS,[],1));
fill_std_known(mean_temp(:,2),std_temp(:,2),100,-99:100,'b','std',2);
mean_temp = squeeze(nanmean(CausalStructure.case2.rDCS,1));
std_temp = squeeze(nanstd(CausalStructure.case3.rDCS,[],1));
fill_std_known(mean_temp(:,2),std_temp(:,2),100,-99:100,'r','std',2);ylim([0,20]);title('TE (ground truth direction)')

subplot(3,5,15);
mean_temp = squeeze(nanmean(CausalStructure.case3.rDCS,1));
std_temp = squeeze(nanstd(CausalStructure.case3.rDCS,[],1));
fill_std_known(mean_temp(:,1),std_temp(:,1),100,-99:100,'b','std',2);
mean_temp = squeeze(nanmean(CausalStructure.case2.rDCS,1));
std_temp = squeeze(nanstd(CausalStructure.case2.rDCS,[],1));
fill_std_known(mean_temp(:,1),std_temp(:,1),100,-99:100,'r','std',2);ylim([0,20]);title('TE (non ground truth direction)')


for i=1:15
    subplot(3,5,i);ylim([0,15]);
    plot([0,0],[0,15],'k--','LineWidth',2);
    xlabel('Time (ms)')
    if i==11 legend('X^1\rightarrowX^2','X^2\rightarrowY^1');legend boxoff 
    end
    if i==12 legend('X^1\rightarrowX^2|X^1','X^2\rightarrowY^1|X^2');legend boxoff
    end
    if i==13 legend('X^1\rightarrowX^2|X^2','X^2\rightarrowY^1|X^2');legend boxoff 
    end
    if i==14 legend('X^2\rightarrowY^1|X^1','X^2\rightarrowY^1|X^2');legend boxoff
    end
    if i==15 legend('X^1\rightarrowX^2|X^1','X^1\rightarrowX^2|X^2');legend boxoff 
    end
end