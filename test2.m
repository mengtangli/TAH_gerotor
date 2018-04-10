clear;
% plot the epitrochoid and parallel
% mengtang li
% Nov 21 2017

clear;clc;

% original input: rho = 43.53, d = 4.029, m = 9, dc = 7.47
rho = 40;
d = 9;
e = d;
m = 4; % Eqn.3: md < rho
dc = 10; % Conclusion of Eqn.6 and 7: dc < R_min
phi = 0:0.02:2*pi;

phi_m = 1/(m-1)*acos((rho^2+m^2*d^2)/(m*d*rho) - 3*(rho^2+m^3*d^2)/(m*d*rho*(m+1)))
R_min_num = (rho^2+m^2*d^2+2*m*d*rho*cos((m-1)*phi_m))^1.5;
R_min_den = (rho^2+m^3*d^2+m*d*rho*(m+1)*cos((m-1)*phi_m));
R_min = R_min_num/R_min_den % R_min = 3.0551 for m = 8, 4.0171 for m = 7
R_num = (rho^2+m^2*d^2+2*m*d*rho*cos((m-1)*phi)).^1.5;
R_den = (rho^2+m^3*d^2+m*d*rho*(m+1)*cos((m-1)*phi));
R = R_num./R_den;
r_min = rho+2*e-dc % Eqn.8
r_max = rho+dc % Eqn.9
r = ceil(r_min)

% figure(1);
% plot(phi, R, 'b', 'LineWidth', 2);
% grid minor; grid on;
% xlim([0 2*pi/(m-1)]);
% ylim([-2 6]);

figure(2);
xo = rho*cos(phi) + d*cos(m*phi); % Eqn.1
yo = rho*sin(phi) + d*sin(m*phi); % Eqn.2
dxo = -rho*sin(phi) - m*d*sin(m*phi);
dyo = rho*cos(phi) + m*d*cos(m*phi);
xi = xo - dc*dyo./sqrt(dxo.^2+dyo.^2); % Eqn.4
yi = yo + dc*dxo./sqrt(dxo.^2+dyo.^2); % Eqn.5
plot(xo,yo,'b','LineWidth',2);
hold on; grid minor; grid on;
plot(xi,yi,'r','LineWidth',2);
axis('equal','xy');
legend('xo and yo', 'xi and yi');
