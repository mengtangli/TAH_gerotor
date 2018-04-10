clear;
% plot the epitrochoid and parallel
% mengtang li
% Nov 21 2017


% original input: 
rho = 43.53, d = 4.029, e = d, m = 9, dc = 7.47;
% rho = 40;
% d = 9;
% e = d;
% m = 4; % Eqn.3: md < rho
% dc = 13; % 2e < dc < R_min
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
r = ceil(r_min)+1

figure(1);
xo = rho*cos(phi) + d*cos(m*phi); % Eqn.1
yo = rho*sin(phi) + d*sin(m*phi); % Eqn.2
dxo = -rho*sin(phi) - m*d*sin(m*phi);
dyo = rho*cos(phi) + m*d*cos(m*phi);
xi = xo - dc*dyo./sqrt(dxo.^2+dyo.^2); % Eqn.4
yi = yo + dc*dxo./sqrt(dxo.^2+dyo.^2); % Eqn.5
plot(xo,yo,'b:','LineWidth',2);
hold on; grid minor; grid on;
plot(xi,yi,'r','LineWidth',2);
axis('equal','xy');
legend('xo and yo', 'xi and yi');

figure(2);
plot(xo+e,yo,'b','LineWidth',2);
hold on; grid minor; grid on;
plot(xi+e,yi,'k','LineWidth',2);
xc = rho*cos(phi); 
yc = rho*sin(phi);
plot(xc,yc,'r--','LineWidth',2);
xr = r*cos(phi); 
yr = r*sin(phi);
plot(xr,yr,'r','LineWidth',2);
d_phi = pi/m;
for i = 1:2:2*m-1
   plot(rho*cos(i*d_phi),rho*sin(i*d_phi),'k+','LineWidth',2);
   x = rho*cos(i*d_phi)+dc*cos(phi);
   y = rho*sin(i*d_phi)+dc*sin(phi);
   plot(x,y,'r--','LineWidth',2);
end
plot(0,0,'k+','LineWidth',2);
text(-2,-2,'C');
plot(e,0,'k+','LineWidth',2);
text(e+1,-2,'O');
axis('equal','xy');

