figure;
% % subplot(2,4,1);
% % plot(-49:50,.5*ones(1,100),'k');hold on;plot(-49:50,Yt_stats_cond_peak1.OLS.At(51:150,1,1),'b');ylabel('Coefficient X^1_{t-1}\rightarrowX^1_{t}|X^1_t');title('Effect\rightarrowEffect|Effect');
% % subplot(2,4,2);
% % plot(-49:50,2*ones(1,100),'k');hold on;plot(-49:50,Yt_stats_cond_peak1.OLS.At(51:150,1,2),'b');ylabel('Coefficient X^2_{t-1}\rightarrowX^1_{t}|X^1_t');title('Cause\rightarrowEffect|Effect');
% % subplot(2,4,3);
% % plot(-49:50,0*ones(1,100),'k');hold on;plot(-49:50,Yt_stats_cond_peak1.OLS.At(51:150,2,1),'g');ylabel('Coefficient X^1_{t-1}\rightarrowX^2_{t}|X^1_t');title('Effect\rightarrowCause|Effect');
% % subplot(2,4,4);
% % plot(-49:50,.5*ones(1,100),'k');hold on;plot(-49:50,Yt_stats_cond_peak1.OLS.At(51:150,2,2),'g');ylabel('Coefficient X^2_{t-1}\rightarrowX^2_{t}|X^1_t');title('Cause\rightarrowCause|Effect');
% % subplot(2,4,5);
% % plot(-49:50,.5*ones(1,100),'k');hold on;plot(-49:50,Yt_stats_cond_peak2.OLS.At(51:150,1,1),'r');ylabel('Coefficient X^1_{t-1}\rightarrowX^1_{t}|X^2_t');title('Effect\rightarrowEffect|Cause');
% % subplot(2,4,6);
% % plot(-49:50,2*ones(1,100),'k');hold on;plot(-49:50,Yt_stats_cond_peak2.OLS.At(51:150,1,2),'r');ylabel('Coefficient X^2_{t-1}\rightarrowX^1_{t}|X^2_t');title('Cause\rightarrowEffect|Cause');
% % subplot(2,4,7);
% % plot(-49:50,0*ones(1,100),'k');hold on;plot(-49:50,Yt_stats_cond_peak2.OLS.At(51:150,2,1),'c');ylabel('Coefficient X^1_{t-1}\rightarrowX^2_{t}|X^2_t');title('Effect\rightarrowCause|Cause');
% % subplot(2,4,8);
% % plot(-49:50,.5*ones(1,100),'k');hold on;plot(-49:50,Yt_stats_cond_peak2.OLS.At(51:150,2,2),'c');ylabel('Coefficient X^2_{t-1}\rightarrowX^2_{t}|X^2_t');title('Cause\rightarrowCause|Cause');

for i=1:8
    subplot(2,4,i);
    xlabel('time (ms)');ylim([-0.1,2.1]);plot([0,0],[-0.1,2.1],'k--')
end


grl = [.8,.8,.8];
bl = 	[0, 0.4470, 0.7410];
or = 	[0.8500, 0.3250, 0.0980];
yel =    	[0.9290, 0.6940, 0.1250];
gr = [0.4660, 0.6740, 0.1880];

for t = 99:99

x0 = squeeze(Yt_events_grdt(1,t-1,:));
y0 = squeeze(Yt_events_grdt(2,t-1,:));
z01 = squeeze(Yt_events_grdt(1,t,:));
z02 = squeeze(Yt_events_grdt(2,t,:));

x1 = squeeze(Yt_events_peak1(1,t-1,:));
y1 = squeeze(Yt_events_peak1(2,t-1,:));
z11 = squeeze(Yt_events_peak1(1,t,:));
z12 = squeeze(Yt_events_peak1(2,t,:));

x2 = squeeze(Yt_events_peak2(1,t-1,:));
y2 = squeeze(Yt_events_peak2(2,t-1,:));
z21 = squeeze(Yt_events_peak2(1,t,:));
z22 = squeeze(Yt_events_peak2(2,t,:));


% figure;
% subplot(2,4,1);plot3(x0,y0,z01,'k.');hold on;plot3(x1,y1,z11,'b.');xlabel('X^1_{t-1}|X^1');ylabel('X^2_{t-1}|X^1');zlabel('X^1_{t}|X^1');title('P(X^1_t|X^1_{t-1},X^2_{t-1}) with S|X^1_t')
% subplot(2,4,2);plot3(x0,y0,z02,'k.');hold on;plot3(x1,y1,z12,'g.');xlabel('X^1_{t-1}|X^1');ylabel('X^2_{t-1}|X^1');zlabel('X^2_{t}|X^1');title('P(X^2_t|X^1_{t-1},X^2_{t-1}) with S|X^1_t')
% subplot(2,4,3);plot3(x0,y0,z01,'k.');hold on;plot3(x2,y2,z21,'r.');xlabel('X^1_{t-1}|X^2');ylabel('X^2_{t-1}|X^2');zlabel('X^1_{t}|X^2');title('P(X^1_t|X^1_{t-1},X^2_{t-1}) with S|X^2_t')
% subplot(2,4,4);plot3(x0,y0,z02,'k.');hold on;plot3(x2,y2,z22,'c.');xlabel('X^1_{t-1}|X^2');ylabel('X^2_{t-1}|X^2');zlabel('X^2_{t}|X^2');title('P(X^2_t|X^1_{t-1},X^2_{t-1}) with S|X^2_t')
% legend('S-independent','S-dependent');
% 
% subplot(2,4,5);plot(y0,z02,'k.');hold on;plot(y1,z12,'g.');xlabel('X^2_{t-1}|X^1');ylabel('X^2_{t}|X^1');title('P(X^2_t|X^2_{t-1}) with S|X^1_t')
% subplot(2,4,6);plot(x0,z02,'k.');hold on;plot(x1,z12,'g.');xlabel('X^1_{t-1}|X^1');ylabel('X^2_{t}|X^1');title('P(X^2_t|X^1_{t-1}) with S|X^1_t')
% subplot(2,4,7);plot(y0,z02,'k.');hold on;plot(y2,z22,'c.');xlabel('X^2_{t-1}|X^2');ylabel('X^2_{t}|X^2');title('P(X^2_t|X^2_{t-1}) with S|X^2_t')
% subplot(2,4,8);plot(x0,z02,'k.');hold on;plot(x2,z22,'c.');xlabel('X^1_{t-1}|X^2');ylabel('X^2_{t}|X^2');title('P(X^2_t|X^1_{t-1}) with S|X^2_t')

figure; set(gcf,'outerposition',[1, 1, 1200*2, 1100])
subplot(2,3,1);plot3(x0,y0,z01,'k.');hold on;view([-112,13]);ax = gca; ax.FontSize = 15;legend('S-independent');legend boxoff;
        xlabel('X^1_{t-1}');ylabel('X^2_{t-1}');zlabel('X^1_{t}');title('P(X^1_t|X^1_{t-1},X^2_{t-1}) (Ground truth)');
subplot(2,3,3);plot3(x0,y0,z01,'.','Color', grl);hold on;plot3(x1,y1,z11,'*','Color',bl);ax = gca; ax.FontSize = 15;legend('S-independent','S-dependent');legend boxoff;
        xlabel('X^1_{t-1}');ylabel('X^2_{t-1}');zlabel('X^1_{t}');title('P(X^1_t|X^1_{t-1},X^2_{t-1}) with S|X^1_t (Aligned by effect)');view([-112,13]);
subplot(2,3,4);plot(y0,z02,'k.'); box off;ax = gca; ax.FontSize = 15;legend('S-independent');legend boxoff;
        xlabel('X^2_{t-1}');ylabel('X^2_{t}');title('P(X^2_t|X^1_{t-1},X^2_{t-1}) (Ground truth)')
subplot(2,3,6);plot(y0,z02,'.','Color', grl);hold on;plot(y1,z12,'*','Color',gr);box off;ax = gca; ax.FontSize = 15;legend('S-independent','S-dependent');legend boxoff;
         xlabel('X^2_{t-1}');ylabel('X^2_{t}');title('P(X^2_t|X^1_{t-1},X^2_{t-1}) with S|X^1_t (Aligned by effect)')
subplot(2,3,2);plot3(x2,y2,z21,'*','Color',or);view([-112,13]);ax = gca; ax.FontSize = 15;legend('S-dependent');legend boxoff;
        xlabel('X^1_{t-1}');ylabel('X^2_{t-1}');zlabel('X^1_{t}');title('P(X^1_t|X^1_{t-1},X^2_{t-1}) with S|X^2_t (Aligned by cause)')
subplot(2,3,5);plot(y0,z02,'.','Color', grl);hold on;plot(y2,z22,'*','Color',yel);box off;ax = gca; ax.FontSize = 15;legend('S-independent','S-dependent');legend boxoff;
         xlabel('X^2_{t-1}');ylabel('X^2_{t}');title('P(X^2_t|X^1_{t-1},X^2_{t-1}) with S|X^2_t (Aligned by cause)')
    

% subplot(2,4,5);plot(y0,z02,'k.');hold on;plot(y1,z12,'g.');xlabel('X^2_{t-1}|X^1');ylabel('X^2_{t}|X^1');title('P(X^2_t|X^2_{t-1}) with S|X^1_t')
% subplot(2,4,6);plot(x0,z02,'k.');hold on;plot(x1,z12,'g.');xlabel('X^1_{t-1}|X^1');ylabel('X^2_{t}|X^1');title('P(X^2_t|X^1_{t-1}) with S|X^1_t')
% subplot(2,4,7);plot(y0,z02,'k.');hold on;plot(y2,z22,'c.');xlabel('X^2_{t-1}|X^2');ylabel('X^2_{t}|X^2');title('P(X^2_t|X^2_{t-1}) with S|X^2_t')
% subplot(2,4,8);plot(x0,z02,'k.');hold on;plot(x2,z22,'c.');xlabel('X^1_{t-1}|X^2');ylabel('X^2_{t}|X^2');title('P(X^2_t|X^1_{t-1}) with S|X^2_t')

%suptitle(strcat('Peri-event time t=',num2str(t-100)));
end

