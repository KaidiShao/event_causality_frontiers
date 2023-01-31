%%ca3 channel 2
figure; set(gcf,'outerposition',[1, 1, 1200*2, 900])
% L=length(y(1,:))-24;t_range=1:length(y(1,:))-24;y_short=y(:,25:end);yf_short=yf;
L = 3500;t_range =1125500+[1:L];y_short = y(:,t_range-24);yf_short = yf(:,t_range);%1125500
Fs=1;

ax1 = subplot(4,8,[1:6]);plot([1:L]/Fs/Fs,y_short(1,:),'k');hold on; ax = gca; ax.FontSize = 15;
ax2 = subplot(4,8,[1:6]+8);plot([1:L]/Fs,y_short(2,:),'k');hold on; ax = gca; ax.FontSize = 15;
ax3 = subplot(4,8,[1:6]+16);plot([1:L]/Fs,yf_short(1,:),'k');hold on;plot([1,L]/Fs,(mean(yf(1,:))+5*std(yf(1,:)))*ones(1,2),'k'); ax = gca; ax.FontSize = 15;
ax4 = subplot(4,8,[1:6]+24);plot([1:L]/Fs,yf_short(2,:),'k');hold on;plot([1,L]/Fs,(mean(yf(2,:))+5*std(yf(2,:)))*ones(1,2),'k'); ax = gca; ax.FontSize = 15;

loc = [1035, 1860];
%%%% event 1
L = 200;t_range_event1 = loc(1)-99:+loc(1)+100;

subplot(4,8,[1:6]);plot(t_range_event1/Fs/Fs,y_short(1,t_range_event1),'b');box off;ylabel('CA3 LFP (mV)');title('Original Signal');
subplot(4,8,[1:6]+8);plot(t_range_event1/Fs,y_short(2,t_range_event1),'b');box off;ylabel('CA1 LFP (mV)'); 
subplot(4,8,[1:6]+16);plot(t_range_event1/Fs,yf_short(1,t_range_event1),'b');box off;ylabel({'Filtered';'CA3 LFP (mV)'});title('Detection Signal (Filtered in [140,230]Hz)'); 
subplot(4,8,[1:6]+24);plot(t_range_event1/Fs,yf_short(2,t_range_event1),'b');box off;ylabel({'Filtered';'CA1 LFP (mV)'});xlabel('Time (ms)');
subplot(4,8,7);plot(t_range_event1/Fs,y_short(1,t_range_event1),'b','LineWidth',1);hold on;box off;ax = gca; ax.FontSize = 15; title('Detected Events')
                plot(loc(1)*ones(1,2)/Fs,ax1.YLim,'k--','LineWidth',1);xlim(t_range_event1(1,[1,end])/Fs);ylim(ax1.YLim);
subplot(4,8,7+8);plot(t_range_event1/Fs,y_short(2,t_range_event1),'b','LineWidth',1);hold on;box off; ax = gca; ax.FontSize = 15;
                plot(loc(1)*ones(1,2)/Fs,ax2.YLim,'k--','LineWidth',1);xlim(t_range_event1(1,[1,end])/Fs);ylim(ax2.YLim)
subplot(4,8,7+16);plot(t_range_event1/Fs,yf_short(1,t_range_event1),'b','LineWidth',1);hold on;plot(t_range_event1([1,end])/Fs,(mean(yf(1,:))+5*std(yf(1,:)))*ones(1,2),'k')
                plot(loc(1)*ones(1,2)/Fs,ax3.YLim,'k--','LineWidth',1);xlim(t_range_event1(1,[1,end])/Fs);ylim(ax3.YLim);
                locs=loc(1);
                plot(locs/Fs,yf_short(1,locs),'b*','LineWidth',2);box off; ax = gca; ax.FontSize = 15;
subplot(4,8,7+24);plot(t_range_event1/Fs,yf_short(2,t_range_event1),'b','LineWidth',1);hold on;plot(t_range_event1([1,end])/Fs,(mean(yf(2,:))+5*std(yf(2,:)))*ones(1,2),'k')
                plot(loc(1)*ones(1,2)/Fs,ax4.YLim,'k--','LineWidth',1);xlim(t_range_event1(1,[1,end])/Fs);ylim(ax4.YLim);box off; ax = gca; ax.FontSize = 15;xlabel('Time (ms)');


%%%%%%%% event 2
L = 200;t_range_event2 = loc(2)-99:+loc(2)+100;
subplot(4,8,[1:6]);plot(t_range_event2/Fs/Fs,y_short(1,t_range_event2),'r')
subplot(4,8,[1:6]+8);plot(t_range_event2/Fs,y_short(2,t_range_event2),'r')
subplot(4,8,[1:6]+16);plot(t_range_event2/Fs,yf_short(1,t_range_event2),'r');
subplot(4,8,[1:6]+24);plot(t_range_event2/Fs,yf_short(2,t_range_event2),'r');

subplot(4,8,8);plot(t_range_event2/Fs,y_short(1,t_range_event2),'r','LineWidth',1);hold on;box off; ax = gca; ax.FontSize = 15;
                plot(loc(2)*ones(1,2)/Fs,ax1.YLim,'k--','LineWidth',1);xlim(t_range_event2(1,[1,end])/Fs);ylim(ax1.YLim)
subplot(4,8,8+8);plot(t_range_event2/Fs,y_short(2,t_range_event2),'r','LineWidth',1);hold on;box off; ax = gca; ax.FontSize = 15;
                plot(loc(2)*ones(1,2)/Fs,ax2.YLim,'k--','LineWidth',1);xlim(t_range_event2(1,[1,end])/Fs);ylim(ax2.YLim)
subplot(4,8,8+16);plot(t_range_event2/Fs,yf_short(1,t_range_event2),'r','LineWidth',1);hold on;plot(t_range_event2([1,end])/Fs,(mean(yf(1,:))+5*std(yf(1,:)))*ones(1,2),'k')
                plot(loc(2)*ones(1,2)/Fs,ax3.YLim,'k--','LineWidth',1);xlim(t_range_event2(1,[1,end])/Fs);ylim(ax3.YLim);
%                 locs=[1830,1831,loc(2),1867,1868,1875,1876,1883];     
                locs=[loc(2)];
                plot(locs/Fs,yf_short(1,locs),'r*','LineWidth',2);box off; ax = gca; ax.FontSize = 15;
subplot(4,8,8+24);plot(t_range_event2/Fs,yf_short(2,t_range_event2),'r','LineWidth',1);hold on;plot(t_range_event2([1,end])/Fs,(mean(yf(2,:))+5*std(yf(2,:)))*ones(1,2),'k')
                plot(loc(2)*ones(1,2)/Fs,ax4.YLim,'k--','LineWidth',1);xlim(t_range_event2(1,[1,end])/Fs);ylim(ax4.YLim);box off; ax = gca; ax.FontSize = 15;xlabel('Time (ms)');