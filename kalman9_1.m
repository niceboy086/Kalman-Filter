% Kalman Filter from the Ground Up - Alex Becker
close all;
clear ; 

Zx =[301.5,298.23,297.83,300.42,301.94,299.5,305.98,301.25, ...
     299.73,299.2,298.62,301.84,299.6,295.3,299.3,301.95,296.3, ...
     295.11,295.12,289.9,283.51,276.42,264.22,250.25,236.66,217.47, ...
     199.75,179.7,160,140.92,113.53,93.68,69.71,45.93,20.87];

Zy =[-401.46,-375.44,-346.15,-320.2,-300.08,-274.12,-253.45,-226.4, ...
     -200.65,-171.62,-152.11,-125.19,-93.4,-74.79,-49.12,-28.73,2.99, ...
     25.65,49.86,72.87,96.34,120.4,144.69,168.06,184.99,205.11, ...
     221.82,238.3,253.02,267.19,270.71,285.86,288.48,292.9,298.77];

Zm = [Zx;Zy];
 
F = [1,1,0.5,0,0,0
     0,1,1,  0,0,0
     0,0,1,  0,0,0
     0,0,0,  1,1,0.5
     0,0,0,  0,1,1
     0,0,0,  0,0,1];
 
Q = [0.25,0.5,0.5,0,0,0
     0.5,1,1,     0,0,0
     0.5,1,1,     0,0,0
     0,0,0,       0.25,0.5,0.5
     0,0,0,       0.5,1,1
     0,0,0,       0.5,1,1] * 0.2^2;
 
Pc = [500,0,0,0,0,0
     0,500,0,0,0,0
     0,0,500,0,0,0
     0,0,0,500,0,0
     0,0,0,0,500,0
     0,0,0,0,0,500];
 
 H = [1,0,0,0,0,0
      0,0,0,1,0,0];
  
 R = [9,0
      0,9];
 
 I = [1,0,0,0,0,0
      0,1,0,0,0,0
      0,0,1,0,0,0
      0,0,0,1,0,0
      0,0,0,0,1,0
      0,0,0,0,0,1];
  
Xc = [0,0,0,0,0,0]';

Est = [];
N = size(Zm,2);
for i=1:N
    Xpre = F*Xc;
    Ppre = F*Pc*(F')+Q;
    Kc = Ppre*(H')*(H*Ppre*(H')+R)^(-1);
    z = Zm(:,i);
    Xc = Xpre + Kc*(z - H*Xpre);
    Pc = (I - Kc*H)*Ppre*((I - Kc*H)') + Kc*R*(Kc');
    if i == 1
        Est = Xc';
    else
        Est = [Est;Xc'];
    end
end

% True Value
Tx(1:17) = 300;
Ty = linspace(-400,0, 17);

An = linspace(0, acos(20/300), 19);
An(1)=[];

Tx = [Tx, cos(An)*300];
Ty = [Ty, sin(An)*300];

Tvx = Tx(2:35) - Tx(1:34);
Tvx = [Tvx(1),Tvx];
Tvy = Ty(2:35) - Ty(1:34);
Tvy = [Tvy(1),Tvy];

figure;
subplot(2,2,[1 3]);
plot(Zx,Zy, '-db','MarkerFaceColor','c','MarkerSize',4);
hold on;
plot(Est(:,1)',Est(:,4)', '-dr','MarkerFaceColor','m','MarkerSize',4);
hold on;
plot(Tx,Ty, '-dg','MarkerFaceColor','w','MarkerSize',4);
legend('measurements','estimates','true value','Location','NorthEastOutside');
xlim([0 350]);
ylim([-450 350]);
daspect([1 1 1])
xlabel('X(m)');
ylabel('Y(m)');
grid on;
set(gca,'TickDir','out')
title('vehicle position')

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
