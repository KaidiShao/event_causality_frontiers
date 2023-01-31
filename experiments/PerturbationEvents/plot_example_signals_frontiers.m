figure; set(gcf,'outerposition',[1, 1, 1200*2, 650])
L=3000;t_range=12501:15500;y_short = x(:,t_range);%yf_short = yf(:,t_range);
Fs=1;

d01=mean(x(1,:))+2.7*std(x(1,:));d02=mean(x(2,:))+3*std(x(2,:));
ax1 = subplot(2,8,[1:6]);plot([1:L]/Fs/Fs,y_short(1,:),'k');hold on; plot([1,L]/Fs,d01*ones(1,2),'k'); ax = gca; ax.FontSize = 20;
ax2 = subplot(2,8,[1:6]+8);plot([1:L]/Fs,y_short(2,:),'k');hold on; plot([1,L]/Fs,d02*ones(1,2),'k'); ax = gca; ax.FontSize = 20;

loc = [485,1913];
%%%% event 1
L = 101;t_range_event1 = loc(1)-50:+loc(1)+50;

subplot(2,8,[1:6]);plot(t_range_event1/Fs/Fs,y_short(1,t_range_event1),'b');box off;ylabel('X^1_t');title('Original Signal');
subplot(2,8,[1:6]+8);plot(t_range_event1/Fs,y_short(2,t_range_event1),'b');box off;ylabel('X^2_t');xlabel('Time (ms)');
subplot(2,8,7);plot(t_range_event1/Fs,y_short(1,t_range_event1),'b','LineWidth',1);hold on;box off;ax = gca; ax.FontSize = 20; title('Detected Events')
                plot(loc(1)*ones(1,2)/Fs,ax1.YLim,'k--','LineWidth',1);xlim(t_range_event1(1,[1,end])/Fs);ylim(ax1.YLim);plot(t_range_event1([1,end])/Fs,d01*ones(1,2),'k')
subplot(2,8,7+8);plot(t_range_event1/Fs,y_short(2,t_range_event1),'b','LineWidth',1);hold on;box off; ax = gca; ax.FontSize = 20;plot(t_range_event1([1,end])/Fs,d02*ones(1,2),'k')
                plot(loc(1)*ones(1,2)/Fs,ax2.YLim,'k--','LineWidth',1);xlim(t_range_event1(1,[1,end])/Fs);ylim(ax2.YLim)
                locs=loc(1);plot(locs/Fs,y_short(2,locs),'b*','LineWidth',2);box off; ax = gca; ax.FontSize = 20;
% subplot(4,8,7+16);plot(t_range_event1/Fs,yf_short(1,t_range_event1),'b','LineWidth',1);hold on;
%                 plot(loc(1)*ones(1,2)/Fs,gca().YLim,'k--','LineWidth',1);xlim(t_range_event1(1,[1,end])/Fs);ylim([-4,4])
%                 locs=[loc(1)-4,loc(1),loc(1)+4];
%                 plot(locs/Fs,yf_short(1,locs),'b*','LineWidth',2);box off; ax = gca; ax.FontSize = 20;
% subplot(4,8,7+24);plot(t_range_event1/Fs,yf_short(2,t_range_event1),'b','LineWidth',1);hold on;
%                 plot(loc(1)*ones(1,2)/Fs,gca().YLim,'k--','LineWidth',1);xlim(t_range_event1(1,[1,end])/Fs);ylim([-4,4]);box off; ax = gca; ax.FontSize = 20;xlabel('Time (ms)');


%%%%%%%% event 2
L = 101;t_range_event2 = loc(2)-50:+loc(2)+50;
subplot(2,8,[1:6]);hold on;plot(t_range_event2/Fs/Fs,y_short(1,t_range_event2),'r')
subplot(2,8,[1:6]+8);hold on;plot(t_range_event2/Fs,y_short(2,t_range_event2),'r')


subplot(2,8,8);plot(t_range_event2/Fs,y_short(1,t_range_event2),'r','LineWidth',1);hold on;box off; ax = gca; ax.FontSize = 20;plot(t_range_event2([1,end])/Fs,d01*ones(1,2),'k')
                plot(loc(2)*ones(1,2)/Fs,ax1.YLim,'k--','LineWidth',1);xlim(t_range_event2(1,[1,end])/Fs);ylim(ax1.YLim)
subplot(2,8,8+8);plot(t_range_event2/Fs,y_short(2,t_range_event2),'r','LineWidth',1);hold on;box off; ax = gca; ax.FontSize = 20;plot(t_range_event2([1,end])/Fs,d02*ones(1,2),'k')
                plot(loc(2)*ones(1,2)/Fs,ax2.YLim,'k--','LineWidth',1);xlim(t_range_event2(1,[1,end])/Fs);ylim(ax2.YLim)
                locs=[1913];
                plot(locs/Fs,y_short(2,locs),'r*','LineWidth',2);box off; ax = gca; ax.FontSize = 20;
% subplot(4,8,8+16);plot(t_range_event2/Fs,yf_short(1,t_range_event2),'r','LineWidth',1);hold on;plot(t_range_event2([1,end])/Fs,(mean(yf(1,:))+5*std(yf(1,:)))*ones(1,2),'k')
%                 plot(loc(2)*ones(1,2)/Fs,gca().YLim,'k--','LineWidth',1);xlim(t_range_event2(1,[1,end])/Fs);ylim([-4,4])
%                 locs=[loc(2),loc(2)+9,loc(2)+14];
%                 plot(locs/Fs,yf_short(1,locs),'r*','LineWidth',2);box off; ax = gca; ax.FontSize = 20;
% subplot(4,8,8+24);plot(t_range_event2/Fs,yf_short(2,t_range_event2),'r','LineWidth',1);hold on;plot(t_range_event2([1,end])/Fs,(mean(yf(2,:))+5*std(yf(2,:)))*ones(1,2),'k')
%                 plot(loc(2)*ones(1,2)/Fs,gca().YLim,'k--','LineWidth',1);xlim(t_range_event2(1,[1,end])/Fs);ylim([-4,4]);box off; ax = gca; ax.FontSize = 20;xlabel('Time (ms)');