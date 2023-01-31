clear

% specify folder path for timeVaryingCausality/Code
folder_path = 'C:\Users\skd\OneDrive\Lab rotation\timeVaryingCausality\Code\';

restoredefaultpath;
addpath(genpath(folder_path)); % other functions related to parameter estimation

Color(1,:) = [0, 0.4470, 0.7410];
Color(2,:) = [0.8500, 0.3250, 0.0980];
Color(3,:) = [0.9290, 0.6940, 0.1250];
Color(4,:) = [0.4940, 0.1840, 0.5560];
Color(5,:) = [0.5, 0.5, 0.5];

y_up = 0.025;
y_lo = 0;

T1 = 200;% 200
T2 = 20;
zeta1 = 0.015722;
zeta2 = .2;

Omega1 = 2*pi/T1;
Omega2 = 2*pi/T2;
gamma1 = 2*zeta1*Omega1;
gamma2 = 2*zeta2*Omega2;

T = 1000;
Ntrial = 2000;

for i = 1:20

for N = 1:Ntrial

x(1) = rand;
x(2) = rand;
y(1) = rand;
y(2) = rand;

h = 1;

% ns_x = 0.02*[ones(1,650) ones(1,201)-1*morlet(-0.29,0.29,201) ones(1,150)];
ns_x = 0.02*ones(1,T+1);
ns_y = 0.005*ones(1,T+1);
c2 = 0; %0.0007 ;
c1 = .098;

for t = 2:T-1
   x(t+1) = (2-gamma1*h) * x(t) + (-1+gamma1*h-h^2*Omega1^2)*x(t-1) + h^2*ns_x(t)*randn + h^2*c2*y(t-1);
   y(t+1) = (2-gamma2*h) * y(t) + (-1+gamma2*h-h^2*Omega2^2)*y(t-1) + h^2*ns_y(t)*randn + h^2*c1*x(t-1); 
end
u = [x;y];
u = u(:,501:end);
X(:,:,N) = u;
end

[bic, mobic] = model_order_selection_compare(X, 10, 'inhomo'); 
[Cs(i,:,:),TE(i,:,:),GC(i,:,:)] = CS_nonzero_mean(X, 2,'inhomo', 0);
end

figure
for n = 1:Ntrial
temp = squeeze(X(:,:,n));
se(n) = sum((temp(1,:)-temp(2,:)).^2);
max_n(n) = max(temp(:));
end
idx1 = se>mean(se);
idx2 = max_n<15;
idx3 = idx1.*idx2;
idx_set = find(idx2==1);

subplot(4,2,1);hold off
plot(squeeze(X(:,:,idx_set(randi(length(idx_set)))))','LineWidth',2);
plot(squeeze(X(:,:,idx_set(randi(length(idx_set)))))','LineWidth',2);
hold on
% legend('X','Y')
title('Example Time Series (control)')
set(gca,'FontSize',20);box off

subplot(4,2,3)
plot(ns_x(501:1000),'LineWidth',2)
hold on
plot(ns_y(501:1000),'LineWidth',2)
% legend('residual variance of X','residual variance of Y')
ylim([y_lo,y_up])
set(gca,'FontSize',20);box off

subplot(4,2,[5 7]);hold off
TE_mean = squeeze(nanmean(TE,1));
CS_mean = squeeze(nanmean(Cs,1));
TE_std = squeeze(nanstd(TE,[],1));
CS_std = squeeze(nanstd(Cs,[],1));

[h(1), h(2)] = plot_shade(1:500, TE_mean(:,1)', TE_std(:,1)', 2, Color(1,:));
hold on
[h(3), h(4)] = plot_shade(1:500, TE_mean(:,2)', TE_std(:,2)', 2, Color(2,:));
[h(5), h(6)] = plot_shade(1:500, CS_mean(:,1)', CS_std(:,1)', 2, Color(3,:));
[h(7), h(8)] = plot_shade(1:500, CS_mean(:,2)', CS_std(:,2)', 2, Color(4,:));
ylim([-0.5 6.5])
plot(150*ones(1,2),get(gca,'Ylim'),'--k')
plot(350*ones(1,2),get(gca,'Ylim'),'--k')

% legend(h(2:2:8),'TE (X \rightarrow Y)','TE (Y \rightarrow X)','CS (X \rightarrow Y)','CS (Y \rightarrow X)')
set(gca,'FontSize',20);box off
%%
for i = 1:20

for N = 1:Ntrial

x(1) = rand;
x(2) = rand;
y(1) = rand;
y(2) = rand;

h = 1;

ns_x = 0.02*[ones(1,650) ones(1,201)-1*morlet(-0.29,0.29,201) ones(1,150)];
% ns_x = 0.02*ones(1,T+1);
ns_y = 0.005*ones(1,T+1);
c2 = 0; %0.0007 ;
c1 = .098;

for t = 2:T-1
   x(t+1) = (2-gamma1*h) * x(t) + (-1+gamma1*h-h^2*Omega1^2)*x(t-1) + h^2*ns_x(t)*randn + h^2*c2*y(t-1);
   y(t+1) = (2-gamma2*h) * y(t) + (-1+gamma2*h-h^2*Omega2^2)*y(t-1) + h^2*ns_y(t)*randn + h^2*c1*x(t-1); 
end
u = [x;y];
u = u(:,501:end);
X(:,:,N) = u;
end

% [bic, mobic] = model_order_selection_compare(X, 10, 'inhomo'); 
%[Cs(i,:,:),TE(i,:,:),GC(i,:,:)] = CS_nonzero_mean(X, 2,'inhomo', 0);
end

for n = 1:Ntrial
temp = squeeze(X(:,:,n));
se(n) = sum((temp(1,:)-temp(2,:)).^2);
max_n(n) = max(temp(:));
end
idx1 = se>mean(se);
idx2 = max_n<15;
idx3 = idx1.*idx2;
idx_set = find(idx2==1);

subplot(4,2,2);hold off
plot(squeeze(X(:,:,idx_set(randi(length(idx_set)))))','LineWidth',2);
% fill_std(squeeze(X(1,:,:)),1:500,'b','std',2)
plot(squeeze(X(:,:,idx_set(randi(length(idx_set)))))','LineWidth',2);
hold on
plot(150*ones(1,2),get(gca,'Ylim'),'--k')
plot(350*ones(1,2),get(gca,'Ylim'),'--k')
legend('x(t)','y(t)');legend boxoff;
title('Example Time Series (Transiently Synchronized)')
set(gca,'FontSize',20);box off

subplot(4,2,4)
plot(ns_x(501:1000),'LineWidth',2)
title('Innovations variance' )
hold on
plot(ns_y(501:1000),'LineWidth',2)
title('Innovations variance' )
ylim([y_lo,y_up])
plot(150*ones(1,2),get(gca,'Ylim'),'--k')
plot(350*ones(1,2),get(gca,'Ylim'),'--k')
legend('x(t)','y(t)');legend boxoff;
set(gca,'FontSize',20);box off

subplot(4,2,[6 8]);hold off
TE_mean = squeeze(nanmean(TE,1));
CS_mean = squeeze(nanmean(Cs,1));
TE_std = squeeze(nanstd(TE,[],1));
CS_std = squeeze(nanstd(Cs,[],1));

[h(1), h(2)] = plot_shade(1:500, TE_mean(:,1)', TE_std(:,1)', 2, Color(1,:));
hold on
[h(3), h(4)] = plot_shade(1:500, TE_mean(:,2)', TE_std(:,2)', 2, Color(2,:));
[h(5), h(6)] = plot_shade(1:500, CS_mean(:,1)', CS_std(:,1)', 2, Color(3,:));
[h(7), h(8)] = plot_shade(1:500, CS_mean(:,2)', CS_std(:,2)', 2, Color(4,:));
ylim([-0.5 6.5])
plot(150*ones(1,2),get(gca,'Ylim'),'--k')
plot(350*ones(1,2),get(gca,'Ylim'),'--k')

legend(h(2:2:8),'TE (x(t) \rightarrow y(t))','TE (y(t) \rightarrow x(t))','DCS (x(t) \rightarrow y(t))','DCS (y(t) \rightarrow x(t))')
set(gca,'FontSize',20);box off; legend boxoff;
%%

cd('D:\Kaidi\Onedrive\Lab rotation\Updated\Results\Figures\Effect of jitter\Data\Figure 1')
set(gcf,'outerposition',get(0,'screensize'));
filename = 'figure_TE_decrease_during_synchrony';
fig_path = strcat(folder_path,'plot\saved_figs\');
addpath('D:\Kaidi\Onedrive\Toolbox\export_fig\')
set(gcf, 'Color', 'w');
saveas(gcf, strcat(fig_path,filename,'.fig'))
saveas(gcf, strcat(fig_path,filename),'epsc')
export_fig(gcf, strcat(fig_path,filename), '-pdf','-png');