L = size(x,2);

D1 = xf(1,:);
d0_1 = mean(D1)+3.2*std(D1);
temp_loc1 = find(D1>=d0_1)-25;
loc_peak1 = find_peak_loc(x(1,:), temp_loc1, 100);

D2 = xf(2,:);
d0_2 =14;
temp_loc2 = find(D2>=d0_2)-25;
loc_peak2 = find_peak_loc(x(2,:), temp_loc2, 100);

xf1=xf(1,26:end);
xf2=xf(2,26:end);
x1=x(1,1:end-25);
x2=x(2,1:end-25);

momax = 4;
tau=1;
% locs = loc;
% locs=
% figure;
% subplot(1,2,1);plot(x(1,:),x(2,:),'k');hold on;plot(x(1,temp_loc1 ),x(2,temp_loc1),'b*','LineWidth',2);plot(x(1,loc),x(2,loc),'y*','LineWidth',2)
% subplot(1,2,2);plot(x(1,:),x(2,:),'k');hold on;plot(x(1,temp_loc2 ),x(2,temp_loc2),'r*','LineWidth',2);plot(x(1,loc),x(2,loc),'y*','LineWidth',2)
i=12;
% figure;
% plot3(x2(1:end-2*i),x2(i+1:end-i),x2(2*i+1:end));hold on;
% plot3(x2(loc+3),x2(loc+i+3),x2(loc+2*i+3),'r*','LineWidth',2);
% plot3(x2(temp_loc1),x2(temp_loc1+i),x2(temp_loc1+2*i),'y*-','LineWidth',2);
% plot3(x2(temp_loc2),x2(temp_loc2+i),x2(temp_loc2+2*i),'g*-','LineWidth',2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

d0 = 14;
or = [0.8500, 0.3250, 0.0980];
orl = [244,185,159]/256;
figure; set(gcf,'outerposition',[1, 1, 1200*2, 900])
subplot(2,4,[1,2,5,6]);
Yt_events= get_Yt(xf, locs+24+2, momax, tau, 100, 201);
temp_event = squeeze(Yt_events(2,101-50:101+50,1:1000));t1=-50:50;
Yt_events_xf = get_Yt(xf, locs+1, momax, tau, 100, 201);
temp_event_xf = squeeze(Yt_events_xf(2,101-50:101+50,1:1000));t1=-50:50;
plot(t1, temp_event,'Color',[.8,.8,.8]);xlim([t1(1),t1(end)]); hold on;
% for i=1:size(temp_event,2)
%     [temp_max,temp_idx] = max(squeeze(temp_event(:,i)));
%     if temp_max>d0
%         plot(t1(temp_idx),temp_event(temp_idx,i),'k*');
%     end
% end
plot([t1(1),t1(end)],ones(1,2)*d0,'Color',or,'LineWidth',2);
Ntrials = size(temp_event,2);
plot([0,0],[-20,20],'k--','LineWidth',2)
idx=find(temp_event(51,:)>=d0);idx2=find(temp_event(51,:)<d0);
plot(zeros(1,length(idx2)), temp_event(51,idx2),'+','Color',orl,'LineWidth',2);
plot(zeros(1,length(idx)), temp_event(51,idx),'+','Color',or);xlabel('Peri-event Time t^{\prime}(ms)');ylabel('X^2_t');ax = gca; ax.FontSize = 15;
title('Event ensemble for cause variable X^2_{t^{\prime}}'); box off;
plot([-2,2,2,-2,-2],[11,11,16.5, 16.5,11],'k')

subplot(2,4,3);
plot([-2,2],ones(1,2)*d0,'Color',or,'LineWidth',2);hold on
temp_event = squeeze(Yt_events(2,101-2:101+2,1:2000));t1=-2:2;
idx=find(temp_event(3,:)>=d0);idx2=find(temp_event(3,:)<d0);
plot(zeros(1,length(idx)), temp_event(3,idx),'+','Color',or);hold on;
plot(zeros(1,length(idx2)), temp_event(3,idx2),'+','Color',orl,'LineWidth',2);

plot(t1, temp_event,'Color',[.8,.8,.8]);xlim([t1(1),t1(end)]); hold on;ylim([11,16.5]);
plot([0,0],[11,16.5],'k--','LineWidth',2)
plot(zeros(1,length(idx)), temp_event(3,idx),'+','Color',or);hold on;
plot(zeros(1,length(idx2)), temp_event(3,idx2),'+','Color',orl,'LineWidth',2);
plot([-2,2],ones(1,2)*d0,'Color',or,'LineWidth',2);xlabel('Peri-event Time t^{\prime}(ms)');ylabel('X^2_t');ax = gca; ax.FontSize = 15;box off;
title('Biased selection at t^{\prime}=0');legend('threshold d_0');legend box off; 

subplot(2,4,4); histogram(temp_event(3,idx),'FaceColor',or,'LineStyle','none','BinWidth',.1); hold on;
histogram(temp_event(3,idx2),'FaceColor',orl,'LineStyle','none','BinWidth',.1);view([90 -90]);ylabel('Counts');xlabel('X^2_t');
xlim([11,16.5]);ax = gca; ax.FontSize = 15; box off; title('Sample Distribution');legend('Selected samples','Unselected samples');legend box off;

bll = [134,206,255]/255;
bl = [0, 0.4470, 0.7410];
subplot(2,4,7);
plot([-2,2],ones(1,2)*d0,'Color',or,'LineWidth',2);hold on;
temp_event = squeeze(Yt_events(2,101-2:101+2,1:2000));t1=-2:2;
idx=find(temp_event(3,:)>=d0);idx2=find(temp_event(3,:)<d0);
plot(zeros(1,length(idx)), temp_event(3,idx),'+','Color',bl);hold on;
plot(zeros(1,length(idx2)), temp_event(3,idx2),'+','Color',bll,'LineWidth',2);
plot(t1, temp_event,'Color',[.8,.8,.8]);xlim([t1(1),t1(end)]); hold on;ylim([11,16.5]);
plot([0,0],[11,16.5],'k--','LineWidth',2)
plot(zeros(1,length(idx)), temp_event(3,idx),'+','Color',bl);hold on;
plot(zeros(1,length(idx2)), temp_event(3,idx2),'+','Color',bll,'LineWidth',2);

idx_pos =find(temp_event(4,:)>=d0);idx_pos2=find(temp_event(4,:)<d0);
plot(ones(1,length(idx_pos)), temp_event(4,idx_pos),'+','Color',bl);hold on;
plot(ones(1,length(idx_pos2)), temp_event(4,idx_pos2),'+','Color',bll,'LineWidth',2);

idx_neg =find(temp_event(2,:)>=d0);idx_neg2=find(temp_event(2,:)<d0);
plot(-ones(1,length(idx_neg)), temp_event(2,idx_neg),'+','Color',bl);hold on;
plot(-ones(1,length(idx_neg2)), temp_event(2,idx_neg2),'+','Color',bll,'LineWidth',2);

% idx_pos21 =find(temp_event(6,:)>=d0);idx_pos22=find(temp_event(6,:)<d0);
% plot(2*ones(1,length(idx_pos21)), temp_event(6,idx_pos21),'+','Color',bl);hold on;
% plot(2*ones(1,length(idx_pos22)), temp_event(6,idx_pos22),'+','Color',bll,'LineWidth',2);
% 
% idx_neg21 =find(temp_event(2,:)>=d0);idx_neg22=find(temp_event(2,:)<d0);
% plot(-2*ones(1,length(idx_neg21)), temp_event(2,idx_neg21),'+','Color',bl);hold on;
% plot(-2*ones(1,length(idx_neg22)), temp_event(2,idx_neg22),'+','Color',bll,'LineWidth',2);

% for i=1:size(temp_event,2)
%     [temp_max(i),temp_idx] = max(squeeze(temp_event(:,i)));
%     if temp_max>d0
%         plot(t1(temp_idx),temp_event(temp_idx,i),'k*');
%     end
% end
plot([-2,2],ones(1,2)*d0,'Color',or,'LineWidth',2);xlabel('Time (ms)');ylabel('X^2_t');ax = gca; ax.FontSize = 15;box off;
title('Biased selection at t^{\prime}=0 (smoothed)');

subplot(2,4,8); 
histogram([temp_event(3,idx),temp_event(4,idx_pos),temp_event(2,idx_neg)],'FaceColor',bl,'LineStyle','none','BinWidth',.1); hold on;
histogram([temp_event(3,idx2),temp_event(4,idx_pos2),temp_event(2,idx_neg2)],'FaceColor',bll,'LineStyle','none','BinWidth',.1);view([90 -90]);ylabel('Counts');xlabel('X^2_t');
xlim([11,16.5]);ax = gca; ax.FontSize = 15; box off; title('Sample Distribution (smoothed)');legend('Selected samples','Unselected samples');legend box off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;plot3(x1(1:end-2*i),x1(i+1:end-i),x1(2*i+1:end));hold on;
plot3(x1(locs+3),x1(locs+i+3),x1(locs+2*i+3),'r*','LineWidth',2);
plot3(x1(temp_loc1),x1(temp_loc1+i),x1(temp_loc1+2*i),'y*-','LineWidth',2);
plot3(x1(temp_loc2+1),x1(temp_loc2+i+1),x1(temp_loc2+2*i+1),'g*-','LineWidth',2);

figure;
bl=[0, 0.4470, 0.7410];
plot3(xf2(1:end-2*i),xf2(i+1:end-i),xf2(2*i+1:end),'Color',[.8,.8,.8]);hold on;
plot3(xf2(locs+2),xf2(locs+i+2),xf2(locs+2*i+2),'*','Color',or,'LineWidth',2);
% plot3(xf2(temp_loc1),xf2(temp_loc1+i),xf2(temp_loc1+2*i),'y*-','LineWidth',2);
plot3(xf2(temp_loc2),xf2(temp_loc2+i),xf2(temp_loc2+2*i),'*','Color',bl,'LineWidth',2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; set(gcf,'outerposition',[1, 1, 1200, 800])
bl=[0, 0.4470, 0.7410];
subplot(2,3,[1,2,4,5]);
plot(xf2(1:end-2*i),xf2(i+1:end-i),'Color',[.8,.8,.8]);hold on;
plot(xf2(locs+2),xf2(locs+i+2),'*','Color',or,'LineWidth',2);
plot(xf2(temp_loc2),xf2(temp_loc2+i),'*','Color',bl,'LineWidth',2);ylim([-17.5,20]);xlim([-17.5,20]);
plot([10,17,17,10,10],[-10,-10,5,5,-10],'k')
title('Reconstructed state-space trajectory');legend('trajectory for all time t','t^{\prime}=0 (ground truth)','t^{\prime}=0 (threshold detection)');legend boxoff
xlabel('X^2_t');ylabel('X^2_{t-\tau}'); ax = gca; ax.FontSize = 15;box off

subplot(2,3,3);
plot(xf2(1:end-2*i),xf2(i+1:end-i),'Color',[.8,.8,.8]);hold on;
plot(xf2(locs+2),xf2(locs+i+2),'*','Color',or,'LineWidth',2);xlim([10,17]);ylim([-10,5]);ax = gca; ax.FontSize = 15;box off;title('Zoomed trajectory');xlabel('X^2_t');ylabel('X^2_{t-\tau}');
subplot(2,3,6);
plot(xf2(1:end-2*i),xf2(i+1:end-i),'Color',[.8,.8,.8]);hold on;
plot(xf2(locs+2),xf2(locs+i+2),'*','Color',or,'LineWidth',2);
plot(xf2(temp_loc2),xf2(temp_loc2+i),'*','Color',bl,'LineWidth',2);xlim([10,17]);ylim([-10,5]);xlabel('X^2_t');ylabel('X^2_{t-\tau}');ax = gca; ax.FontSize = 15;box off;title('Zoomed trajectory')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; set(gcf,'outerposition',[1, 1, 1200, 800])
bl=[0, 0.4470, 0.7410];
subplot(2,3,[1,2,4,5]);
plot(x2(1:end-2*i),x2(i+1:end-i),'Color',[.8,.8,.8]);hold on;
plot(x2(locs+2),x2(locs+i+2),'*','Color',or,'LineWidth',2);
plot(x2(temp_loc2),xf2(temp_loc2+i),'*','Color',bl,'LineWidth',2);ylim([-17.5,20]);xlim([-17.5,20]);
plot([10,17,17,10,10],[-10,-10,5,5,-10],'k')
title('Reconstructed state-space trajectory');legend('trajectory for all time t','t^{\prime}=0 (ground truth)','t^{\prime}=0 (threshold detection)');legend boxoff
xlabel('X^2_t');ylabel('X^2_{t-\tau}'); ax = gca; ax.FontSize = 15;box off

subplot(2,3,3);
plot(x2(1:end-2*i),xf2(i+1:end-i),'Color',[.8,.8,.8]);hold on;
plot(xf2(locs+2),xf2(locs+i+2),'*','Color',or,'LineWidth',2);xlim([10,17]);ylim([-10,5]);ax = gca; ax.FontSize = 15;box off;title('Zoomed trajectory');xlabel('X^2_t');ylabel('X^2_{t-\tau}');
subplot(2,3,6);
plot(xf2(1:end-2*i),xf2(i+1:end-i),'Color',[.8,.8,.8]);hold on;
plot(xf2(locs+2),xf2(locs+i+2),'*','Color',or,'LineWidth',2);
plot(xf2(temp_loc2),xf2(temp_loc2+i),'*','Color',bl,'LineWidth',2);xlim([10,17]);ylim([-10,5]);xlabel('X^2_t');ylabel('X^2_{t-\tau}');ax = gca; ax.FontSize = 15;box off;title('Zoomed trajectory')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;plot3(xf1(1:end-2*i),xf1(i+1:end-i),xf1(2*i+1:end));hold on;
plot3(xf1(locs+3),xf1(locs+i+3),xf1(locs+2*i+3),'r*','LineWidth',2);
plot3(xf1(temp_loc1),xf1(temp_loc1+i),xf1(temp_loc1+2*i),'y*-','LineWidth',2);
plot3(xf1(temp_loc2+2),xf1(temp_loc2+i+2),xf1(temp_loc2+2*i+2),'g*-','LineWidth',2);

xf1_smoth = smooth(xf1,2);
xf2_smoth = smooth(xf2,2);
D1 = xf1_smoth;
d0_1 = mean(D1)+3.2*std(D1);
temp_loc1 = find(D1>=d0_1)-25;

D2 = xf2_smoth;
d0_2 = mean(D2)+3.2*std(D2);
temp_loc2 = find(D2>=d0_2)-25;


figure;
plot3(xf2_smoth(1:end-2*i),xf2_smoth(i+1:end-i),xf2_smoth(2*i+1:end));hold on;
plot3(xf2_smoth(locs+3),xf2(locs+i+3),xf2_smoth(locs+2*i+3),'r*','LineWidth',2);
plot3(xf2_smoth(temp_loc1),xf2(temp_loc1+i),xf2_smoth(temp_loc1+2*i),'y*-','LineWidth',2);
plot3(xf2_smoth(temp_loc2),xf2(temp_loc2+i),xf2_smoth(temp_loc2+2*i),'g*-','LineWidth',2);

figure;plot3(xf1_smoth(1:end-2*i),xf1_smoth(i+1:end-i),xf1_smoth(2*i+1:end));hold on;
plot3(xf1_smoth(locs+3),xf1_smoth(locs+i+3),xf1_smoth(locs+2*i+3),'r*','LineWidth',2);
plot3(xf1_smoth(temp_loc1),xf1_smoth(temp_loc1+i),xf1_smoth(temp_loc1+2*i),'y*-','LineWidth',2);
plot3(xf1_smoth(temp_loc2+2),xf1_smoth(temp_loc2+i+2),xf1_smoth(temp_loc2+2*i+2),'g*-','LineWidth',2);











% 
% D1 = xf(1,:);
% d0_1 = mean(D1)+4*std(D1);
% temp_loc1 = find(D1>=d0_1);
% loc_peak1 = find_peak_loc(x(1,:), temp_loc1, 100);
% 
% D2 = xf(2,:);
% d0_2 = mean(D2)+4*std(D2);
% temp_loc2 = find(D2>=d0_2);
% loc_peak2 = find_peak_loc(x(2,:), temp_loc2, 100);
% 
% mobic = 1; tau = 1;
% Yt_events_grdt = get_Yt(x, locs, mobic, tau, L_start, L_extract);
% Yt_stats_cond_grdt = get_Yt_stats_cond(Yt_events_grdt,mobic);
% 
% Yt_events_peak1 = get_Yt(x, loc_peak1, mobic, tau, L_start, L_extract);
% Yt_stats_cond1_full = get_Yt_stats_cond(Yt_events_peak1,mobic);
% Yt_events_peak1(1,:,:)=[];
% Yt_stats_cond1 = get_Yt_stats_cond_uni(Yt_events_peak1,mobic);
% % Yt_events_peak1(2,:,:)=[];
% % Yt_stats_cond1_part = get_Yt_stats_cond_uni(Yt_events_peak1,mobic);
% 
% Yt_events_peak2 = get_Yt(x, loc_peak2, mobic, tau, L_start, L_extract);
% Yt_stats_cond2_full = get_Yt_stats_cond(Yt_events_peak2,mobic);
% Yt_events_peak2(2,:,:)=[];
% Yt_stats_cond2 = get_Yt_stats_cond_uni(Yt_events_peak2,mobic);
% 
% figure;
% subplot(1,2,1); plot(squeeze(Yt_stats_cond1_full.OLS.At(:,2,:)));hold on;plot(squeeze(Yt_stats_cond1.OLS.At));
% subplot(1,2,2); plot(squeeze(Yt_stats_cond2_full.OLS.At(:,1,:)));hold on;plot(squeeze(Yt_stats_cond2.OLS.At));


function Yt_stats_cond = get_Yt_stats_cond_uni(Yt_event,mo)
nvar = 1;
Yt_stats_cond.mean = mean(Yt_event,3);
for t = 1:size(Yt_event, 2)
    temp = squeeze(Yt_event(:,t,:))-Yt_stats_cond.mean(:,t);
    Yt_stats_cond.Sigma(t,:,:) = temp*temp'/size(Yt_event,3);
    Yt_stats_cond.OLS.At(t,:,:) = reshape(squeeze(Yt_stats_cond.Sigma(t,1:nvar,nvar+1:end)),1,2*mo)...
                                 /squeeze(Yt_stats_cond.Sigma(t,nvar+1:end,nvar+1:end));
end

[Yt_stats_cond.OLS.bt, Yt_stats_cond.OLS.Sigma_Et, Yt_stats_cond.OLS.sigma_Et] = estimate_residuals_uni(Yt_stats_cond);
end

function [bt, Sigma_Et, sigma_Et] = estimate_residuals_uni(Yt_stats)
[L, nvar, ~] = size(Yt_stats.OLS.At);
nvar = 1;

bt = nan(nvar, L);
Sigma_Et = nan(L, nvar, nvar); 
sigma_Et = nan(L, 1); 
for t = 1:L
    
    Sigma_Xt = squeeze(Yt_stats.Sigma(t,1:nvar,1:nvar));
    Sigma_Xp = squeeze(Yt_stats.Sigma(t,nvar+1:end,nvar+1:end));
    Sigma_XtXp = reshape(squeeze(Yt_stats.Sigma(t,1:nvar,nvar+1:end)),1,2);
    coeff = reshape(squeeze(Yt_stats.OLS.At(t,:,:)),1,2);
    
    bt(:,t) = Yt_stats.mean(1:nvar,t)-coeff*Yt_stats.mean(nvar+1:end,t);
%     Sigma_Et(t,:,:) = Sigma_Xt -coeff*Sigma_Xp*coeff';
    Sigma_Et(t,:,:) = Sigma_Xt - Sigma_XtXp*coeff' - coeff*Sigma_XtXp' +coeff*Sigma_Xp*coeff';
%     Sigma_Et(t,:,:) = .5*(squeeze(Sigma_Et(t,:,:))+squeeze(Sigma_Et(t,:,:))');
%     if det(squeeze(Sigma_Et(t,:,:)))<0 && disp_flag
%         disp(strcat('Non Semi-definite time point: t = ',int2str(t)))
%     end
    sigma_Et(t,1) = trace(squeeze(Sigma_Et(t,:,:)));
end

Sigma_Et(Sigma_Et<0)=0;
sigma_Et(sigma_Et<0)=0;

end

% subplot(1,2,1);
% patch(x(1,1:end-1),x(2,1:end-1),x(1,2:end),x(2,2:end),'EdgeColor','k','FaceColor','none','EdgeAlpha',.1);hold on
% patch(x_stat(1,1:end-1),x_stat(2,1:end-1),x_stat(1,2:end),x_stat(2,2:end),'EdgeColor','c','FaceColor','none','EdgeAlpha',.1);view(27,-3)
% plot3(x(1,loc),x(2,loc),x(1,loc+1),'b*');
% plot3(x(1,loc_peak1),x(2,loc_peak1),x(1,loc_peak1+1),'r*');
% 
% 
% subplot(1,2,2);
% patch(x(1,1:end-1),x(2,1:end-1),x(2,2:end),x(2,2:end),'EdgeColor','k','FaceColor','none','EdgeAlpha',.1);hold on;
% patch(x_stat(1,1:end-1),x_stat(2,1:end-1),x_stat(2,2:end),x_stat(2,2:end),'EdgeColor','c','FaceColor','none','EdgeAlpha',.1);view(27,-3)
% plot3(x(1,loc),x(2,loc),x(2,loc+1),'b*');
% plot3(x(1,loc_peak2),x(2,loc_peak2),x(2,loc_peak2+1),'r*');
% 
% figure;
% subplot(1,2,1);
% patch(x(1,1:end-1),x(1,2:end),x(1,2:end),'EdgeColor','k','FaceColor','none','EdgeAlpha',.1);hold on
% patch(x_stat(1,1:end-1),x_stat(1,2:end),x_stat(2,2:end),'EdgeColor','c','FaceColor','none','EdgeAlpha',.1);
% plot(x(1,loc),x(1,loc+1),'b*');
% 
% subplot(1,2,2);
% patch(x(1,1:end-1),x(1,2:end),x(1,2:end),'EdgeColor','k','FaceColor','none','EdgeAlpha',.1);hold on
% patch(x_stat(1,1:end-1),x_stat(1,2:end),x_stat(2,2:end),'EdgeColor','c','FaceColor','none','EdgeAlpha',.1);
% plot(x(1,loc),x(1,loc+1),'b*');
% plot(x(1,loc_peak1),x(1,loc_peak1+1),'r*');
% 
% figure;
% subplot(1,2,1);
% patch(x(1,1:end-1),x(2,2:end),x(2,2:end),'EdgeColor','k','FaceColor','none','EdgeAlpha',.1);hold on
% patch(x_stat(1,1:end-1),x_stat(2,2:end),x_stat(2,2:end),'EdgeColor','c','FaceColor','none','EdgeAlpha',.1);
% plot(x(1,loc),x(2,loc+1),'b*');
% 
% subplot(1,2,2);
% patch(x(1,1:end-1),x(2,2:end),x(1,2:end),'EdgeColor','k','FaceColor','none','EdgeAlpha',.1);hold on
% patch(x_stat(1,1:end-1),x_stat(2,2:end),x_stat(2,2:end),'EdgeColor','c','FaceColor','none','EdgeAlpha',.1);
% plot(x(1,loc),x(2,loc+1),'b*');
% plot(x(1,loc_peak2),x(2,loc_peak2+1),'r*');
% 
% figure;
% subplot(1,2,2);
% patch(x(1,1:end-1),x(2,1:end-1),x(2,2:end),x(2,2:end),'EdgeColor','k','FaceColor','none','EdgeAlpha',.1);hold on;
% patch(x_stat(1,1:end-1),x_stat(2,1:end-1),x_stat(2,2:end),x_stat(2,2:end),'EdgeColor','c','FaceColor','none','EdgeAlpha',.1);view(27,-3)
% plot3(x(1,loc),x(2,loc),x(2,loc+1),'b*');
% plot3(x(1,loc_peak2),x(2,loc_peak2),x(2,loc_peak2+1),'r*');
