clear;
% close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DEFINE PATHS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
onedrive_path = 'C:\Users\skd\OneDrive\';
% onedrive_path = 'D:\Kaidi\OneDrive\';
proj_path = strcat(onedrive_path,'\updated-desnap-with-causality\');
save_path = strcat(proj_path,'\opendata_hc3\saved_data\newtest9\');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DATASET INFO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_chpair = 1024;
L_event = 201;
align_mode = {'ca3_peak','ca1_peak'};
causal_methods = {'TE', 'DCS', 'rDCS'};
sess_name = 'vvp01_2006-4-9_18-43-47';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ALLOCATE SPACE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n = 1:size(align_mode, 2)
for k = 1:size(causal_methods,2)
    VarAlign = matlab.lang.makeValidName(align_mode{n});
    VarCausalM = matlab.lang.makeValidName(causal_methods{k});
    CausalStructure.OLS.(VarAlign).(VarCausalM) = nan(n_chpair, L_event, 2); 
    CausalStructure_btsp.OLS.(VarAlign).(VarCausalM) = nan(n_chpair, L_event, 2, 100); 
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% load data for different conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n_ch = [1:n_chpair]
    save_path_chpair = strcat(save_path,'chpair_no_',int2str(n_ch),'\');
for n = 1:size(align_mode,2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% load saved data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename = strcat(save_path_chpair, 'ca3_ca1_',sess_name,'_chpair_no_', int2str(n_ch), '_', align_mode{n}, '_model_causality.mat');
load(filename);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% extract Yt_stats %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
VarAlign = matlab.lang.makeValidName(align_mode{n});
Yt_stats_all.mean.(VarAlign)(n_ch,:,:) = Yt_stats.mean(1:2,:)';
Yt_stats_all.Sigma.(VarAlign)(n_ch,:,1) = Yt_stats.Sigma(:,1,1)'/Yt_stats.Ntrials;
Yt_stats_all.Sigma.(VarAlign)(n_ch,:,2) = Yt_stats.Sigma(:,2,2)'/Yt_stats.Ntrials;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% extract causality results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k = 1:size(causal_methods,2)
VarCausalM = matlab.lang.makeValidName(causal_methods{k});
CausalStructure.OLS.(VarAlign).(VarCausalM)(n_ch,:,:) = CausalOutput.OLS.(VarCausalM);   
if Params.Options.Bootstrap
    for n_btsp = 1:Params.MonteC_Params.Nbtsp
        load(strcat(save_path_chpair, 'ca3_ca1_',sess_name,'_chpair_no_', int2str(n_ch),'_', align_mode{n}, '_btsp_',int2str(n_btsp), '_model_causality'));  
        CausalStructure_btsp.OLS.(VarAlign).(VarCausalM)(n_ch,:,:,n_btsp) = CausalOutput_btsp.OLS.(VarCausalM);
    end
end
end

end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plot for Yt stats (different alignments) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; set(gcf,'outerposition',[1, 1, 1200, 900])
mean_var1_temp = squeeze(nanmean(Yt_stats_all.mean.(mkname(align_mode{1}))));    
std_var1_temp = sqrt(squeeze(nanmean(Yt_stats_all.Sigma.(mkname(align_mode{1})))));
mean_var2_temp = squeeze(nanmean(Yt_stats_all.mean.(mkname(align_mode{2}))));    
std_var2_temp = sqrt(squeeze(nanmean(Yt_stats_all.Sigma.(mkname(align_mode{2})))));

t = -100:100;
t_idx = 1:201;
Fs = 1252;
subplot(2,2,1);
fill_std_known(mean_var1_temp(t_idx,1), std_var2_temp(t_idx,1), 2000, t/Fs, 'b', 'std', 2);
subplot(2,2,2);
fill_std_known(mean_var1_temp(t_idx,2), std_var2_temp(t_idx,2), 2000, t/Fs, 'r', 'std', 2);
subplot(2,2,3);
fill_std_known(mean_var2_temp(t_idx,1), std_var2_temp(t_idx,1), 2000, t/Fs, 'b', 'std', 2);
subplot(2,2,4);
fill_std_known(mean_var2_temp(t_idx,2), std_var2_temp(t_idx,2), 2000, t/Fs, 'r', 'std', 2);

event_titles = {'(sr|sr)', '(pl|sr)', '(sr|pl)', '(pl|pl)'};
for i = 1:4
    subplot(2,2,i);xlim(t(1,[1,end])/Fs); 
    ax = gca; ax.FontSize = 15; 
    temp = ax.YLim;
    plot([0,0]/Fs,ax.YLim,'k--','LineWidth',2);
    title(strcat('Event waveforms', event_titles{i}));
    %xlim([-50,50]/667);
    if i>2 xlabel('Peri-event time t^{\prime} (ms)'); end
    ylabel('LFP (mV)'); box off;
    ylim(temp);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plot for Yt stats (different alignments) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; set(gcf,'outerposition',[1, 1, 1200, 900])
y_max = 0.4;
condition_label = {'Aligned by putative cause', 'Aligned by putative effect'};
legend_Col = {'ca3\rightarrowca1|ca3','ca1\rightarrowca3|ca1';
              'ca3\rightarrowca1|ca1','ca1\rightarrowca3|ca3'};

t = -100:100;
for k = 1:size(causal_methods,2)
    for i_Col = 1:2
    subplot(3,2,(k-1)*2+i_Col)
    Varname3 = matlab.lang.makeValidName(causal_methods{k});

    switch i_Col

        case 1
            mean_temp1 = squeeze(nanmean(CausalStructure.OLS.(mkname(align_mode{1})).(Varname3)(:,t_idx,1),1))';
%             std_temp1 = squeeze(nanstd(CausalStructure.OLS.sr_peak.(Varname3)(:,t_idx,1),1))';
            std_temp1 = squeeze(nanmean(nanstd(CausalStructure_btsp.OLS.(mkname(align_mode{1})).(Varname3)(:,t_idx,1,:),[],4),1));
            fill_std_known(mean_temp1, std_temp1, 200, t/Fs, 'b', 'std', 2);hold on;
            
            mean_temp2 = squeeze(nanmean(CausalStructure.OLS.(mkname(align_mode{2})).(Varname3)(:,t_idx,2),1))';
            std_temp2 = squeeze(nanmean(nanstd(CausalStructure_btsp.OLS.(mkname(align_mode{2})).(Varname3)(:,t_idx,2,:),[],4),1))';
%             std_temp2 = squeeze(nanstd(CausalStructure.OLS.pl_peak.(Varname3)(:,t_idx,2),1))';
            fill_std_known(mean_temp2, std_temp2, 1000, t/Fs, 'r', 'std', 2);hold on;
            
    
        case 2 
            mean_temp1 = squeeze(nanmean(CausalStructure.OLS.(mkname(align_mode{2})).(Varname3)(:,t_idx,1),1))';
%             std_temp1 = squeeze(nanstd(CausalStructure.OLS.pl_peak.(Varname3)(:,t_idx,1),1))';
            std_temp1 = squeeze(nanmean(nanstd(CausalStructure_btsp.OLS.(mkname(align_mode{2})).(Varname3)(:,t_idx,1,:),[],4),1))';
            fill_std_known(mean_temp1, std_temp1, 1000, t/Fs, 'b', 'std', 2);hold on;
            
            mean_temp2 = squeeze(nanmean(CausalStructure.OLS.(mkname(align_mode{1})).(Varname3)(:,t_idx,2),1))';
            std_temp2 = squeeze(nanmean(nanstd(CausalStructure_btsp.OLS.(mkname(align_mode{1})).(Varname3)(:,t_idx,2,:),[],4),1))';
%             std_temp2 = squeeze(nanstd(CausalStructure.OLS.sr_peak.(Varname3)(:,t_idx,2),1))';
            fill_std_known(mean_temp2, std_temp2, 1000, t/Fs, 'r', 'std', 2);hold on;
    end
  
    ylim([0,y_max]); xlim(t(1,[1,end])/Fs); 
    ylabel(causal_methods{k});
    if k==1 title(condition_label{i_Col}); end
    plot(zeros(1,2),[0,y_max],'k--','LineWidth',2)
    if k==1 l = legend(legend_Col{i_Col,:}); set(l,'color','none'); legend boxoff; xlabel('Time (ms)'); end 
      ax = gca; ax.FontSize = 15; box off;
    if k==3 xlabel('Peri-event time t^{\prime} (ms)');end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% gaussian confidence interval
figure; set(gcf,'outerposition',[1, 1, 1200, 900])
y_max = 0.75;
condition_label = {'Aligned by putative cause', 'Aligned by putative effect'};
legend_Col = {'ca3\rightarrowca1|ca3','ca1\rightarrowca3|ca1';
              'ca3\rightarrowca1|ca1','ca1\rightarrowca3|ca3'};

t = -100:100;
for k = 1:size(causal_methods,2)
    for i_Col = 1:2
    subplot(3,2,(k-1)*2+i_Col)
    Varname3 = matlab.lang.makeValidName(causal_methods{k});

    switch i_Col

        case 1
            mean_temp1 = squeeze(nanmean(CausalStructure.OLS.(mkname(align_mode{1})).(Varname3)(:,t_idx,1),1))';
%             std_temp1 = squeeze(nanstd(CausalStructure.OLS.sr_peak.(Varname3)(:,t_idx,1),1))';
            mean_btsp_temp1 = squeeze(nanmean(nanmean(CausalStructure_btsp.OLS.(mkname(align_mode{1})).(Varname3)(:,t_idx,1,:),4),1));
            std_btsp_temp1 = squeeze(nanmean(nanstd(CausalStructure_btsp.OLS.(mkname(align_mode{1})).(Varname3)(:,t_idx,1,:),[],4),1));
            fill_ci_btsp(mean_temp1, mean_btsp_temp1, std_btsp_temp1, 100, t/Fs, 'b', 2);hold on;
            
            mean_temp2 = squeeze(nanmean(CausalStructure.OLS.(mkname(align_mode{2})).(Varname3)(:,t_idx,2),1))';
            mean_btsp_temp2 = squeeze(nanmean(nanmean(CausalStructure_btsp.OLS.(mkname(align_mode{2})).(Varname3)(:,t_idx,2,:),4),1))';
            std_btsp_temp2 = squeeze(nanmean(nanstd(CausalStructure_btsp.OLS.(mkname(align_mode{2})).(Varname3)(:,t_idx,2,:),[],4),1))';
%             std_temp2 = squeeze(nanstd(CausalStructure.OLS.pl_peak.(Varname3)(:,t_idx,2),1))';
            fill_ci_btsp(mean_temp2, mean_btsp_temp2, std_btsp_temp2,100, t/Fs, 'r', 2);hold on;
            patch(t/Fs,.65*ones(1,201),100*sum(squeeze(h.(mkname(strcat('pos',int2str(k))))(1,:,:)),1)/1024,'FaceColor','none','EdgeColor','interp','LineStyle','none','Marker','*');
    
        case 2 
            mean_temp1 = squeeze(nanmean(CausalStructure.OLS.(mkname(align_mode{2})).(Varname3)(:,t_idx,1),1))';
%             std_temp1 = squeeze(nanstd(CausalStructure.OLS.sr_peak.(Varname3)(:,t_idx,1),1))';
            mean_btsp_temp1 = squeeze(nanmean(nanmean(CausalStructure_btsp.OLS.(mkname(align_mode{2})).(Varname3)(:,t_idx,1,:),4),1));
            std_btsp_temp1 = squeeze(nanmean(nanstd(CausalStructure_btsp.OLS.(mkname(align_mode{2})).(Varname3)(:,t_idx,1,:),[],4),1));
            fill_ci_btsp(mean_temp1, mean_btsp_temp1, std_btsp_temp1, 100, t/Fs, 'b', 2);hold on;
            
            mean_temp2 = squeeze(nanmean(CausalStructure.OLS.(mkname(align_mode{1})).(Varname3)(:,t_idx,2),1))';
            mean_btsp_temp2 = squeeze(nanmean(nanmean(CausalStructure_btsp.OLS.(mkname(align_mode{1})).(Varname3)(:,t_idx,2,:),4),1))';
            std_btsp_temp2 = squeeze(nanmean(nanstd(CausalStructure_btsp.OLS.(mkname(align_mode{1})).(Varname3)(:,t_idx,2,:),[],4),1))';
%             std_temp2 = squeeze(nanstd(CausalStructure.OLS.pl_peak.(Varname3)(:,t_idx,2),1))';
            fill_ci_btsp(mean_temp2, mean_btsp_temp2, std_btsp_temp2, 100, t/Fs, 'r', 2);hold on;

            patch(t/Fs,.65*ones(1,201),100*sum(squeeze(h.(mkname(strcat('pos',int2str(k))))(2,:,:)),1)/1024,'FaceColor','none','EdgeColor','interp','LineStyle','none','Marker','*');
            
    end
    
    caxis([0,100])
    ylim([0,y_max]); xlim(t(1,[1,end])/Fs); 
    ylabel(causal_methods{k});
    if k==1 title(condition_label{i_Col}); end
    plot(zeros(1,2),[0,y_max],'k--','LineWidth',2)
    if k==1 l = legend(legend_Col{i_Col,:}); set(l,'color','none'); legend boxoff; xlabel('Time (ms)'); end 
      ax = gca; ax.FontSize = 15; box off;
    if k==3 xlabel('Peri-event time t^{\prime} (ms)');end
    end
end
colormap(othercolor('Set34'));a=colorbar;a.Label.String = '% of significance';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% empirical confidence interval
figure; set(gcf,'outerposition',[1, 1, 1200, 900])
y_max = 0.75;
condition_label = {'Aligned by putative cause', 'Aligned by putative effect'};
legend_Col = {'ca3\rightarrowca1|ca3','ca1\rightarrowca3|ca1';
              'ca3\rightarrowca1|ca1','ca1\rightarrowca3|ca3'};

t = -100:100;
for k = 1:size(causal_methods,2)
    for i_Col = 1:2
    subplot(3,2,(k-1)*2+i_Col)
    Varname3 = matlab.lang.makeValidName(causal_methods{k});

    switch i_Col

        case 1
            mean_temp1 = squeeze(nanmean(CausalStructure.OLS.(mkname(align_mode{1})).(Varname3)(:,t_idx,1),1))';
            for tt=t_idx
                for i_chpair = 1:1024
                    pop1 = squeeze(CausalStructure_btsp.OLS.(mkname(align_mode{1})).(Varname3)(i_chpair,tt,1,:));
                    [CI_temp(i_chpair,1),CI_temp(i_chpair,2)] = ci_ecdf(pop1);
                end
                CI_t(tt,:) = mean(CI_temp,1);
            end
            fill_known_ci_btsp(mean_temp1, CI_t, 100, t/Fs, 'b', 2);hold on;

            mean_temp2 = squeeze(nanmean(CausalStructure.OLS.(mkname(align_mode{2})).(Varname3)(:,t_idx,2),1))';
            for tt=t_idx
                for i_chpair = 1:1024
                    pop1 = squeeze(CausalStructure_btsp.OLS.(mkname(align_mode{2})).(Varname3)(i_chpair,tt,2,:));
                    [CI_temp(i_chpair,1),CI_temp(i_chpair,2)] = ci_ecdf(pop1);
                end
                CI_t(tt,:) = mean(CI_temp,1);
            end
            fill_known_ci_btsp(mean_temp2, CI_t, 100, t/Fs, 'r', 2);hold on;
            patch(t/Fs,.65*ones(1,201),100*sum(squeeze(h.(mkname(strcat('pos',int2str(k))))(1,:,:)),1)/1024,'FaceColor','none','EdgeColor','interp','LineStyle','none','Marker','*');
    
        case 2 
            mean_temp1 = squeeze(nanmean(CausalStructure.OLS.(mkname(align_mode{2})).(Varname3)(:,t_idx,1),1))';
            for tt=t_idx
                for i_chpair = 1:1024
                    pop1 = squeeze(CausalStructure_btsp.OLS.(mkname(align_mode{2})).(Varname3)(i_chpair,tt,1,:));
                    [CI_temp(i_chpair,1),CI_temp(i_chpair,2)] = ci_ecdf(pop1);
                end
                CI_t(tt,:) = mean(CI_temp,1);
            end
            fill_known_ci_btsp(mean_temp1, CI_t, 100, t/Fs, 'b', 2);hold on;

            mean_temp2 = squeeze(nanmean(CausalStructure.OLS.(mkname(align_mode{1})).(Varname3)(:,t_idx,2),1))';
            for tt=t_idx
                for i_chpair = 1:1024
                    pop1 = squeeze(CausalStructure_btsp.OLS.(mkname(align_mode{1})).(Varname3)(i_chpair,tt,2,:));
                    [CI_temp(i_chpair,1),CI_temp(i_chpair,2)] = ci_ecdf(pop1);
                end
                CI_t(tt,:) = mean(CI_temp,1);
            end
            fill_known_ci_btsp(mean_temp2, CI_t, 100, t/Fs, 'r', 2);hold on;
            patch(t/Fs,.65*ones(1,201),100*sum(squeeze(h.(mkname(strcat('pos',int2str(k))))(2,:,:)),1)/1024,'FaceColor','none','EdgeColor','interp','LineStyle','none','Marker','*');
            
    end
  
    ylim([0,y_max]); xlim(t(1,[1,end])/Fs); 
    ylabel(causal_methods{k});
    caxis([0,100])
    if k==1 title(condition_label{i_Col}); end
    plot(zeros(1,2),[0,y_max],'k--','LineWidth',2)
    if k==1 l = legend(legend_Col{i_Col,:}); set(l,'color','none'); legend boxoff; xlabel('Time (ms)'); end 
      ax = gca; ax.FontSize = 15; box off;
    if k==3 xlabel('Peri-event time t^{\prime} (ms)');end
    end
end
colormap(othercolor('Set34'));a=colorbar;a.Label.String = '% of significance';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% significance figure (ranksum)
for k = 1:size(causal_methods,2)
for i_align = 1:2
subplot(3,2,(k-1)*2+i_Col)
Varname3 = matlab.lang.makeValidName(causal_methods{k});
for n_chpair = 1:1024
for t=t_idx
pop1 = squeeze(CausalStructure_btsp.OLS.(mkname(align_mode{1})).(Varname3)(n_chpair,t,1,:));
pop2 = squeeze(CausalStructure_btsp.OLS.(mkname(align_mode{2})).(Varname3)(n_chpair,t,2,:));
mean_pop1 = squeeze(CausalStructure.OLS.(mkname(align_mode{1})).(Varname3)(n_chpair,t,1));
mean_pop2 = squeeze(CausalStructure.OLS.(mkname(align_mode{2})).(Varname3)(n_chpair,t,2));
[h, p] = ranksum(pop1,pop2);
h = h*sign(mean_pop1-mean_pop2);
CausalStructure_btsp.OLS.sigtest.(mkname(strcat(causal_methods{k},'_sig')))(1,n_chpair,t) = h;
CausalStructure_btsp.OLS.sigtest.(mkname(strcat(causal_methods{k},'_pvalue')))(1,n_chpair,t) = p;
pop3 = squeeze(CausalStructure_btsp.OLS.(mkname(align_mode{2})).(Varname3)(n_chpair,t,1,:));
pop4 = squeeze(CausalStructure_btsp.OLS.(mkname(align_mode{1})).(Varname3)(n_chpair,t,2,:));
mean_pop3 = squeeze(CausalStructure.OLS.(mkname(align_mode{2})).(Varname3)(n_chpair,t,1));
mean_pop4 = squeeze(CausalStructure.OLS.(mkname(align_mode{1})).(Varname3)(n_chpair,t,2));
[h, p] = ranksum(pop3,pop4);
h = h*sign(mean_pop3-mean_pop4);
CausalStructure_btsp.OLS.sigtest.(mkname(strcat(causal_methods{k},'_sig')))(2,n_chpair,t) = h;
CausalStructure_btsp.OLS.sigtest.(mkname(strcat(causal_methods{k},'_pvalue')))(2,n_chpair,t) = p;
end
end
end
end

clear h
figure; set(gcf,'outerposition',[1, 1, 1800, 900])
h1 = CausalStructure_btsp.OLS.sigtest.TE_sig;h.pos1 = h1>0;
h2 = CausalStructure_btsp.OLS.sigtest.DCS_sig;h.pos2 = h2>0;
h3 = CausalStructure_btsp.OLS.sigtest.rDCS_sig;h.pos3 = h3>0;
subplot(3,3,1)
plot(-100:100,100*sum(squeeze(h.pos1(1,:,:)),1)/1024);hold on;
plot(-100:100,100*sum(squeeze(h.pos1(2,:,:)),1)/1024);
subplot(3,3,2)
imagesc(-100:100,1:1024,squeeze(h.pos1(1,:,:)));hold on;
subplot(3,3,3)
imagesc(-100:100,1:1024,squeeze(h.pos1(2,:,:)));hold on;
subplot(3,3,4)
plot(-100:100,100*sum(squeeze(h.pos2(1,:,:)),1)/1024);hold on;
plot(-100:100,100*sum(squeeze(h.pos2(2,:,:)),1)/1024);
subplot(3,3,5)
imagesc(-100:100,1:1024,squeeze(h.pos2(1,:,:)));hold on;
subplot(3,3,6)
imagesc(-100:100,1:1024,squeeze(h.pos2(2,:,:)));hold on;
subplot(3,3,7)
plot(-100:100,100*sum(squeeze(h.pos3(1,:,:)),1)/1024);hold on;
plot(-100:100,100*sum(squeeze(h.pos3(2,:,:)),1)/1024);
subplot(3,3,8)
imagesc(-100:100,1:1024,squeeze(h.pos3(1,:,:)));hold on;
subplot(3,3,9)
imagesc(-100:100,1:1024,squeeze(h.pos3(2,:,:)));hold on;

for i=1:9
    subplot(3,3,i);
    ax = gca; y_max = ax.YLim(2); plot(zeros(1,2),[0,y_max],'k--','LineWidth',1)
    ax.FontSize = 15; box off;
    if i==9 colorbar('east'); end
    if i==1 title('Significant Percentage'); end
    if mod(i,3)==1 ylabel('% of channel pairs'); end
    if i==2 || i==3 title('Timewise Significance'); end
    if mod(i,3)==2 || mod(i,3)==3 ylabel('no. of channel pair'); end
    if i>6 xlabel('Peri-event Time (ms)'); end
    if i==7 legend('Aligned by putative cause', 'Aligned by putative effects'); legend boxoff; end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% significance figure (ttest2)
for k = 1:size(causal_methods,2)
for i_align = 1:2
subplot(3,2,(k-1)*2+i_Col)
Varname3 = matlab.lang.makeValidName(causal_methods{k});
for n_chpair = 1:1024
for t=t_idx
pop1 = squeeze(CausalStructure_btsp.OLS.(mkname(align_mode{1})).(Varname3)(n_chpair,t,1,:));
pop2 = squeeze(CausalStructure_btsp.OLS.(mkname(align_mode{2})).(Varname3)(n_chpair,t,2,:));
mean_pop1 = squeeze(CausalStructure.OLS.(mkname(align_mode{1})).(Varname3)(n_chpair,t,1));
mean_pop2 = squeeze(CausalStructure.OLS.(mkname(align_mode{2})).(Varname3)(n_chpair,t,2));
[h, p] = ttest2(pop1,pop2);
h = h*sign(mean_pop1-mean_pop2);
CausalStructure_btsp.OLS.sigtest.(mkname(strcat(causal_methods{k},'_sig')))(1,n_chpair,t) = h;
CausalStructure_btsp.OLS.sigtest.(mkname(strcat(causal_methods{k},'_pvalue')))(1,n_chpair,t) = p;
pop3 = squeeze(CausalStructure_btsp.OLS.(mkname(align_mode{2})).(Varname3)(n_chpair,t,1,:));
pop4 = squeeze(CausalStructure_btsp.OLS.(mkname(align_mode{1})).(Varname3)(n_chpair,t,2,:));
mean_pop3 = squeeze(CausalStructure.OLS.(mkname(align_mode{2})).(Varname3)(n_chpair,t,1));
mean_pop4 = squeeze(CausalStructure.OLS.(mkname(align_mode{1})).(Varname3)(n_chpair,t,2));
[h, p] = ttest2(pop3,pop4);
h = h*sign(mean_pop3-mean_pop4);
CausalStructure_btsp.OLS.sigtest.(mkname(strcat(causal_methods{k},'_sig')))(2,n_chpair,t) = h;
CausalStructure_btsp.OLS.sigtest.(mkname(strcat(causal_methods{k},'_pvalue')))(2,n_chpair,t) = p;
end
end
end
end

clear h
figure; set(gcf,'outerposition',[1, 1, 1800, 900])
h1 = CausalStructure_btsp.OLS.sigtest.TE_sig;h.pos1 = h1>0;
h2 = CausalStructure_btsp.OLS.sigtest.DCS_sig;h.pos2 = h2>0;
h3 = CausalStructure_btsp.OLS.sigtest.rDCS_sig;h.pos3 = h3>0;
subplot(3,3,1)
plot(-100:100,100*sum(squeeze(h.pos1(1,:,:)),1)/1024);hold on;
plot(-100:100,100*sum(squeeze(h.pos1(2,:,:)),1)/1024);
subplot(3,3,2)
imagesc(-100:100,1:1024,squeeze(h.pos1(1,:,:)));hold on;
subplot(3,3,3)
imagesc(-100:100,1:1024,squeeze(h.pos1(2,:,:)));hold on;
subplot(3,3,4)
plot(-100:100,100*sum(squeeze(h.pos2(1,:,:)),1)/1024);hold on;
plot(-100:100,100*sum(squeeze(h.pos2(2,:,:)),1)/1024);
subplot(3,3,5)
imagesc(-100:100,1:1024,squeeze(h.pos2(1,:,:)));hold on;
subplot(3,3,6)
imagesc(-100:100,1:1024,squeeze(h.pos2(2,:,:)));hold on;
subplot(3,3,7)
plot(-100:100,100*sum(squeeze(h.pos3(1,:,:)),1)/1024);hold on;
plot(-100:100,100*sum(squeeze(h.pos3(2,:,:)),1)/1024);
subplot(3,3,8)
imagesc(-100:100,1:1024,squeeze(h.pos3(1,:,:)));hold on;
subplot(3,3,9)
imagesc(-100:100,1:1024,squeeze(h.pos3(2,:,:)));hold on;

for i=1:9
    subplot(3,3,i);
    ax = gca; y_max = ax.YLim(2); plot(zeros(1,2),[0,y_max],'k--','LineWidth',1)
    ax.FontSize = 15; box off;
    if i==9 colorbar('east'); end
    if i==1 title('Significant Percentage'); end
    if mod(i,3)==1 ylabel('% of channel pairs'); end
    if i==2 || i==3 title('Timewise Significance'); end
    if mod(i,3)==2 || mod(i,3)==3 ylabel('no. of channel pair'); end
    if i>6 xlabel('Peri-event Time (ms)'); end
    if i==7 legend('Aligned by putative cause', 'Aligned by putative effects'); legend boxoff; end
end
