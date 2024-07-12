% Kalman Filter from the Ground Up - Alex Becker
close all;
clear ; 

Zr =[502.55,477.34,457.21,442.94,427.27,406.05,400.73,377.32,...
     360.27,345.93,333.34,328.07,315.48,301.41,302.87,304.25,294.46,...
     294.29,299.38,299.37,300.68,304.1,301.96,300.3,301.9,296.7,...
     297.07,295.29,296.31,300.62,292.3,298.11,298.07,298.92,298.04];

Zan =[-0.9316,-0.8977,-0.8512,-0.8114,-0.7853,-0.7392,-0.7052,-0.6478,...
      -0.59,-0.5183,-0.4698,-0.3952,-0.3026,-0.2445,-0.1626,-0.0937,0.0085,...
      0.0856,0.1675,0.2467,0.329,0.4149,0.504,0.5934,0.667,0.7537,...
      0.8354,0.9195,1.0039,1.0923,1.1546,1.2564,1.3274,1.409,1.5011];

Zm = [Zr;Zan];

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
 
 R = [5^2,0
      0,0.0087^2];

N = 6;
W = [-1 1/6 1/6 1/6 1/6 1/6 1/6 1/6 1/6 1/6 1/6 1/6 1/6];
% I = eye(2 * N + 1);
% Wdiag = diag(W);

Pc = [500,0,0,0,0,0
     0,500,0,0,0,0
     0,0,500,0,0,0
     0,0,0,500,0,0
     0,0,0,0,500,0
     0,0,0,0,0,500];
 
Xc = [400,0,0,-300,0,0]';
Est = [];
NN = size(Zm,2);
for ii = 1:NN
    L = chol(3 * Pc)';

    Xsgm(:,1) = Xc;
    for i = 1:N
        Xsgm(:,i+1) = Xc + L(:,i);
    end
    for i = 1:N
        Xsgm(:,i+N+1) = Xc - L(:,i);
    end

    Xsgmp = F*Xsgm;
    Xpre = Xsgmp * W';
    Ppre = (Xsgmp - repmat(Xpre,1,2*N+1)) * diag(W) * (Xsgmp - repmat(Xpre,1,2*N+1))' + Q;
    r = Xsgmp(1,:).^2 + Xsgmp(4,:).^2;
    Hpre = [sqrt(r); atan(Xsgmp(4,:) ./ Xsgmp(1,:))];
    Hpm =   Hpre * W';
    Php = (Hpre - repmat(Hpm,1,2*N+1)) * diag(W) * (Hpre - repmat(Hpm,1,2*N+1))' + R;
    Pxh = (Xsgmp - repmat(Xpre,1,2*N+1)) * diag(W) * (Hpre - repmat(Hpm,1,2*N+1))';

    Kc = Pxh * Php^(-1);
    Xc = Xpre + Kc * (Zm(:,ii) - Hpm);
    Pc = Ppre - Kc * Php * Kc';
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
plot(Zr.*cos(Zan),Zr.*sin(Zan), '-db','MarkerFaceColor','c','MarkerSize',4);
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
ylim([-70 30]);
%daspect([0.05 1 1])
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
ylim([-100 100]);
%daspect([0.04 1 1])
xlabel('time(s)');
ylabel('Vy(m/s)');
grid on;
set(gca,'TickDir','out');
title('vehicle Y-axis velocity');





