% Err = [0.0061 0.0033 0.0023 0.0018 0.0014]';
% t = [60.7915 121.5830 182.3745 243.1661 303.9576]';
% % logerr = log10(Err);
% % logt = log10(t);
% 
% figure(5)
% % plot(Err , t ,'-b*', 'linewidth' , 2)
% plot(t , Err ,'-b*', 'linewidth' , 2)
% xlabel('time')
% ylabel('error')
% title('||error||/||T_{exact}|| vs time')
% % ylabel('LOG (||error||/||T_{exact}||)')
% set(gca,'FontName','Times New Roman','FontSize',10,'fontWeight','bold');
% grid on


% FOR TABLE #2 IN THE REPORT
t = [10 50 100 150 200]';
% for dt = 1 sec.
Err1 = [0.0050 0.0019 0.0011 8.2241e-04 6.7226e-04]';
% for dt = 3 sec.
Err25 =[0.0049 0.0019 9.7290e-04 6.2518e-04 4.3195e-04];
% for dt = 5 sec.
Err5 = [0.0203 0.0054 0.0030 0.0021 0.0016]';

figure(5)
plot(t , Err1 ,'-g*', 'linewidth' , 2)
xlabel('time')
ylabel('error')
title('||error||/||T_{exact}|| vs time')
set(gca,'FontName','Times New Roman','FontSize',10,'fontWeight','bold');
hold on
plot(t , Err25 ,'-ro', 'linewidth' , 2)
hold on
plot(t , Err5 ,'-ks', 'linewidth' , 2)
legend('\Delta t = 1 sec.' , '\Delta t = 2.5 sec.' , '\Delta t = 5 sec.')
hold off
grid on

