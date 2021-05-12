clc
clear all
close all
MONT = 500;
T = 0.5;             % sampling Time
total_time = 20;
t =0: T: total_time;        % Time
nstates = 6;         % numbers of States
nmeas = 3;           % number of measurements 
init_pos = [15; 15 ; 200]; %  inital Position
V = [7; 7; 0];
TotalScans=length(t);
R=[9 0 0;
    0 1 0;
    0 0 0.1];
x(:,1) = [init_pos;V ]+[randn(3,1);16*randn(3,1)]; % Set initial state estimate
P0 = [9*eye(3,3) zeros(3,3); zeros(3,3) 16^2*eye(3,3)]; % Set initial error covariance
noise =sqrt(R)*randn(3, TotalScans); % Generate random measurement noise 
xt(:,1) = [15; 15; 200; 7; 7; 0]; %+ sqrt(P0)*randn(6,1); % Set true initial state
ynoisy = zeros(3, TotalScans); % Initialize size of output vector for all k
for k = 2:TotalScans
xt(1,k) = xt(1,k-1)+T*V(1);
xt(2,k) = xt(2,k-1)+T*V(2);
xt(3,k) = xt(3,k-1)+T*V(3);
         
end
xt(:,1)=[15;15;200;7;7;0];
for k = 1:TotalScans
    ytrue(:,k) = [sqrt(xt(1,k)^2 + xt(2,k)^2 + xt(3,k)^2); ...
           atan2d(xt(2,k),xt(1,k)); ...
           atan2d(xt(3,k),(sqrt(xt(1,k)^2 +xt(2,k)^2)))] 
ynoisy(:,k) = ytrue(:,k) + noise(:,k);           
end
xest = zeros(6,TotalScans);
F=[1 0 0 T 0 0;          % Transition matrix
    0 1 0 0 T 0;
    0 0 1 0 0 T;
    0 0 0 1 0 0;
    0 0 0 0 1 0;
    0 0 0 0 0 1];
h(:,1)=[201,45,54]
for k = 1:TotalScans
% Prediction
if k==1
    xest(:,1) = x(:,1);
P = P0;
else
x_pred = F*xest(:,k-1);
P_pred = F*P*F';
 x=x_pred(1);
 y=x_pred(2);
 z=x_pred(3);
% Observation
y_est(:,k) = [sqrt(x_pred(1)^2 + x_pred(2)^2 + x_pred(3)^2); ...
atan2d(x_pred(2),x_pred(1)); ...
atan2d(x_pred(3),(sqrt(x_pred(1)^2 +x_pred(2)^2)))];
 %%% Jacobian Matrix  %%%%%
 H=[x/(x^2 + y^2 + z^2)^(1/2),                        y/(x^2 + y^2 + z^2)^(1/2),                   z/(x^2 + y^2 + z^2)^(1/2), 0, 0, 0;
  -y/(x^2*(y^2/x^2 + 1)),                              1/(x*(y^2/x^2 + 1)),                                           0, 0, 0, 0;
  -(x*z)/((z^2/(x^2 + y^2 + z^2) + 1)*(x^2 + y^2 + z^2)^(3/2)), -(y*z)/((z^2/(x^2 + y^2 + z^2) + 1)*(x^2 + y^2 + z^2)^(3/2)), (1/(x^2 + y^2 + z^2)^(1/2) - z^2/(x^2 + y^2 + z^2)^(3/2))/(z^2/(x^2 + y^2 + z^2) + 1), 0, 0, 0];
% Measurement Update
S=(H*P_pred*H' + R)
K = P_pred*H'*inv(S) % Calculate Kalman gain
xest(:,k) = x_pred + K*(ynoisy(:,k) - y_est(:,k)); % Update state estimate
P = (eye(nstates)-K*H)*P_pred; % Update covariance estimate
h(:,k)=[sqrt(xest(1,k)^2 + xest(2,k)^2 + xest(3,k)^2); ...
atan2d(xest(2,k),xest(1,k)); ...
atan2d(xest(3,k),(sqrt(xest(1,k)^2 +xest(2,k)^2)))];
end
end

figure, 
plot(t,xest(1,:),'r','LineWidth',1.5)
hold on
plot(t,xt(1,:),'g','LineWidth',1.5)
legend('Estimated Position on x axis', 'True Position on x axis')
xlabel('Time (s)')
ylabel('position (m)')
title('Target Position on x axis ')
grid on


figure, 
plot(t,xest(2,:),'LineWidth',1.5)
hold on
plot(t,xt(2,:),'g','LineWidth',1.5)
legend('Estimated Position on y axis', 'True Position on y axis')
xlabel('Time (s)')
ylabel('position (m)')
title('Target Position on y axis ')
grid on

figure, 
plot(t,xest(3,:),'r','LineWidth',1.5)
hold on
plot(t,xt(3,:),'g','LineWidth',1.5)
legend('Estimated Position on z axis', 'True Position on z axis')
xlabel('Time (s)')
ylabel('position (m)')
title('Target Position on z axis ')
grid on
% figure, 
% plot(xe(4,:))
% 
% figure, 
% plot(xe(5,:))
% 
% figure, 
% plot(xe(6,:))

figure, plot(ynoisy(1,2:TotalScans),'r','LineWidth',1.5)
hold on
plot(h(1,2:TotalScans),'g','LineWidth',1.5)
hold on
plot(ytrue(1,2:TotalScans),'k','LineWidth',1.5)
legend('Estimated range','Measured range','true range')
xlabel('Time (s)')
ylabel('Range (m)')
title('Target Range')
grid on

figure, plot(ynoisy(2,2:TotalScans),'r','LineWidth',1.5)
hold on
plot(h(2,2:TotalScans),'g','LineWidth',1.5)
hold on
plot(ytrue(2,2:TotalScans),'k','LineWidth',1.5)
legend('Estimated azimuth','Measured azimuth','True azimuth')
xlabel('Time (s)')
ylabel('azimuth (degree)')
title('Target azimuth')
grid on


figure, plot(ynoisy(3,2:TotalScans),'r','LineWidth',1.5)
hold on
plot(h(3,2:TotalScans),'g','LineWidth',1.5)
hold on
plot(ytrue(3,2:TotalScans),'k','LineWidth',1.5)
legend('Estimated elevation','Measured elevation','True elevation')
xlabel('Time (s)')
ylabel('elevation (degree)')
title('Target elevation')
grid on
%RMSE
figure
mse=(h(1,:)-ytrue(1,:)).^2;
plot(mse(3:TotalScans))
xlabel('Time (s)')
ylabel('MSE')
title('MSE of position Estimate')
grid on
for i=1:TotalScans
Vx=V-xest(4:6,i);
msev(i)=sqrt(Vx(1)^2+Vx(2)^2+Vx(3)^2);
end
figure
plot(msev)
grid on
xlabel('Time (s)')
ylabel('MSE')
title('MSE of Velocity Estiamte')