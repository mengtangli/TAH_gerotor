% try to calculate the total area A_i. Dont know how to deal with Eqn.27
% mengtang li
% Nov 21 2017

clear;clc;

rho = 10;
d = 1;
m = 7; % Eqn.3: md < rho
dc = 2; % Conclusion of Eqn.6 and 7: dc < R_min
e = d; % Conclusion above Fig.3 and 4
r_min = rho+2*e-dc; % 10 for m = 7
r_max = rho+dc; % 12 for m = 7
r = 11;
theta = 0:0.05:3*pi; % 1x189
n = size(theta,2); % 189
angle = 0:0.05:3*pi;

syms xo yo dxo dyo ddxo ddyo xi yi dxi dyi;
for i = 1:1:n
    phi = angle(i);
    xo = rho*cos(phi) + d*cos(m*phi); % Eqn.1
    yo = rho*sin(phi) + d*sin(m*phi); % Eqn.2
    dxo = -rho*sin(phi) - m*d*sin(m*phi);
    dyo = rho*cos(phi) + m*d*cos(m*phi);
    ddxo = -rho*cos(phi) - m^2*d*cos(m*phi);
    ddyo = -rho*sin(phi) - m^2*d*sin(m*phi);
    xi = xo - dc*dyo/sqrt(dxo^2+dyo^2); % Eqn.4
    yi = yo + dc*dxo/sqrt(dxo^2+dyo^2); % Eqn.5
    dxi = dxo - dc/(dxo^2+dyo^2)^1.5*(ddyo*(dxo^2+dyo^2)-dyo*(dxo+dyo));
    dyi = dyo + dc/(dxo^2+dyo^2)^1.5*(ddxo*(dxo^2+dyo^2)-dxo*(dxo+dyo));
    
    f(i) = xi*dyi - dxi*yi; % Eqn.27
end
plot(angle, f, 'b', 'LineWidth', 2); grid minor; grid on;
% phi_s = acos(-m*e/rho)/(m-1) % Eqn.30
% f_theta_s =  atan(dyo/dxo) % Eqn.31
% theta_s = double(subs(f_theta_s, phi, phi_s)) % Eqn.31
% theta_i = theta(1);
% phi_L_i = -phi_s/(pi+theta_s)*(theta_i+theta_s)+phi_s; % Eqn.28
% phi_F_i = -phi_s/(pi+theta_s)*(theta_i+theta_s+2*pi)+phi_s; % Eqn.29
% A_I_i = 0.5*int(f, phi, phi_L_i, phi_F_i) % Eqn.27






