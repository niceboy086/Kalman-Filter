% Kalman Filter from the Ground Up (2023) -Alex Becker.pdf
close all;
clear ; 
% the set of 30 noisy measurements of the altitude hn and acceleration an:
Zh =[6.43,1.3,39.43,45.89,41.44,48.7,78.06,80.08, ...
     61.77,75.15,110.39,127.83,158.75,156.55,213.32,229.82,262.8, ...
     297.57,335.69,367.92,377.19,411.18,460.7,468.39,553.9,583.97, ...
     655.15,723.09,736.85,787.22];

Za =[39.81,39.67,39.81,39.84,40.05,39.85,39.78,39.65, ...
     39.67,39.78,39.59,39.87,39.85,39.59,39.84,39.9,39.63, ...
     39.59,39.76,39.79,39.73,39.93,39.83,39.85,39.94,39.86, ...
     39.76,39.86,39.74,39.94];
% the measurements period
dtime = 0.25;
% the altimeter measurement error standard deviation
dm = 20;
% the accelerometer measurement error standard deviation
epsilon = 0.1;
% the state transition matrix
F = [1,dtime
     0,1];
% the control matrix
G = [dtime^2/2; dtime];
% the process noise matrix
Q = [dtime^4/4, dtime^3/2
     dtime^3/2, dtime^2] * epsilon^2;
% the measurement uncertainty
R = [dm^2];

% 
Xc = [0,0]';
Uc = 9.8;

Pc = [500,0
      0,500];
 
H = [1,0];
  
I = [1,0
     0,1];

Est = [];
N = size(Zh,2);
% init  
Xpre = F*Xc + G*(Uc - 9.8);
Ppre = F*Pc*(F')+Q;
for i=1:N
    z = Zh(i);
    Uc = Za(i);
    Kc = Ppre*(H')*(H*Ppre*(H')+R)^(-1);
    Xc = Xpre + Kc*(z - H*Xpre);
    Pc = (I - Kc*H)*Ppre*((I - Kc*H)') + Kc*R*(Kc');
    if i == 1
        Est = Xc';
    else
        Est = [Est;Xc'];
    end
    Xpre = F*Xc + G*(Uc - 9.8);
    Ppre = F*Pc*(F')+Q;
end


% True Value
Th = [];
Tv = [];
for i=1:30
    Th(i) = 30*((i-1)*dtime)^2/2 + 15;
    Tv(i) = 30*((i-1)*dtime);
end

figure;
subplot(1,1,1);
plot(Zh, '-db','MarkerFaceColor','c','MarkerSize',4);
hold on;
plot(Est(:,1)', '-dr','MarkerFaceColor','m','MarkerSize',4);
hold on;
plot(Th, '-dg','MarkerFaceColor','w','MarkerSize',4);
legend('measurements','estimates','true value','Location','NorthEastOutside');
%xlim([0 350]);
ylim([-50 850]);
%daspect([1 1 1])
xlabel('X(m)');
ylabel('Y(m)');
grid on;
set(gca,'TickDir','out')
title('vehicle position')
%{
subplot(2,2,2);
plot(Tvx, '-dg','MarkerFaceColor','w','MarkerSize',4);
hold on;
plot(Est(:,2)', '-dr','MarkerFaceColor','m','MarkerSize',4);
legend('true value','estimates','Location','NorthEast');
xlim([0 36]);
ylim([-100 250]);
daspect([0.05 1 1])
xlabel('time(s)');
ylabel('Vx(m/s)');
grid on;
set(gca,'TickDir','out');
title('vehicle X-axis velocity')
subplot(2,2,4);
plot(Tvy, '-dg','MarkerFaceColor','w','MarkerSize',4);
hold on;
plot(Est(:,5)', '-dr','MarkerFaceColor','m','MarkerSize',4);
legend('true value','estimates','Location','SouthEast');
xlim([0 36]);
ylim([-300 100]);
daspect([0.04 1 1])
xlabel('time(s)');
ylabel('Vy(m/s)');
grid on;
set(gca,'TickDir','out')
title('vehicle Y-axis velocity')
%}