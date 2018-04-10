% plot the epitrochoid
% mengtang li
% Nov 21 2017

clear;clc;

% rho = 43.53, d = 4.029, e = d, m = 9, dc = 7.47;
rho = 40;
d = 9;
e = d;
m = 4; % Eqn.3: md < rho
dc = 13; % Conclusion of Eqn.6 and 7: dc < R_min

phi = 0:0.02:2*pi;
xr = rho*cos(phi);
yr = rho*sin(phi);
x = rho*cos(phi) + d*cos(m*phi);
y = rho*sin(phi) + d*sin(m*phi);
figure(1); clf;
plot(xr,yr,'r','LineWidth',2);
hold on; grid minor; grid on;
plot(x,y,'b','LineWidth',2);
axis equal;
% xlim([-15 15]);
