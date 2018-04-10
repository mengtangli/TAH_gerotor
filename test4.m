% Try to calculate the total area A_i. Now it works.
% But not sure whether all of those equations are correct.
% mengtang li
% Nov 21 2017

clear;clc;

% rho = 10;
% d = 1;
% m = 7; % Eqn.3: md < rho
% dc = 2; % Conclusion of Eqn.6 and 7: dc < R_min
% e = d; % Conclusion above Fig.3 and 4
% r_min = rho+2*e-dc; % 10 for m = 7
% r_max = rho+dc; % 12 for m = 7
% r = 11;
% theta = 0:0.05:3*pi; % 1x189
% n = size(theta,2); % 189
% theta_i = theta(1);

syms rho d m dc xo yo dxo dyo ddxo ddyo xi yi dxi dyi phi;

xo = rho*cos(phi) + d*cos(m*phi) % Eqn.1
yo = rho*sin(phi) + d*sin(m*phi) % Eqn.2
dxo = -rho*sin(phi) - m*d*sin(m*phi)
dyo = rho*cos(phi) + m*d*cos(m*phi)
% diff_xo = diff(xo,phi) % Correct
% diff_yo = diff(yo,phi) % Correct
ddxo = -rho*cos(phi) - m^2*d*cos(m*phi)
ddyo = -rho*sin(phi) - m^2*d*sin(m*phi)
% diff2_xo = diff(diff_xo,phi) % Correct
% diff2_yo = diff(diff_yo,phi) % Correct
xi = xo - dc*dyo/sqrt(dxo^2+dyo^2) % Eqn.4
yi = yo + dc*dxo/sqrt(dxo^2+dyo^2) % Eqn.5
dxi = dxo - dc/(dxo^2+dyo^2)^1.5*(ddyo*(dxo^2+dyo^2)-dyo*(dxo+dyo))
dyi = dyo + dc/(dxo^2+dyo^2)^1.5*(ddxo*(dxo^2+dyo^2)-dxo*(dxo+dyo))
diff_xi = diff(xi,phi) % Not the same
diff_yi = diff(yi,phi) % Not the same

% f1 = xi*dyi
% f2 = dxi*yi
f = xi*dyi - dxi*yi % Eqn.27
f2 = xi*diff_yi - diff_xi*yi % Eqn.27

% phi_s = acos(-m*e/rho)/(m-1) % Eqn.30
f_theta_s = atan(dyo/dxo) % Eqn.31
% theta_s = double(subs(f_theta_s, phi, phi_s)) % Eqn.31
% phi_L_i = -phi_s/(pi+theta_s)*(theta_i+theta_s)+phi_s % Eqn.28
% phi_F_i = -phi_s/(pi+theta_s)*(theta_i+theta_s+2*pi)+phi_s % Eqn.29
% A_I_i = double(0.5*int(f, phi, phi_L_i, phi_F_i)) % Eqn.27
% 
% A_s = pi*r^2/m % Eqn.16
% A_t = pi*dc^2 % Eqn.17
% beta_F_i = asin((m*e*sin((m-1)/m*theta_i-pi/m))/...
%     sqrt(m^2*e^2+rho^2-2*m*e*rho*cos((m-1)/m*theta_i-pi/m))) % Eqn.26
% beta_L_i = asin((m*e*sin((m-1)/m*theta_i+pi/m))/...
%     sqrt(m^2*e^2+rho^2-2*m*e*rho*cos((m-1)/m*theta_i+pi/m))) % Eqn.25
% A_SF_i = beta_F_i*dc^2/2 % Eqn.24
% A_TF_i = rho*dc*sin(beta_F_i)/2 % Eqn.23
% A_SL_i = beta_L_i*dc^2/2 % Eqn.22
% A_TL_i = rho*dc*sin(beta_L_i)/2 % Eqn.21
% k = (-rho^2+2*dc*rho-dc^2+r^2)/(2*rho) % Eqn.20
% c = sqrt(4*r^2-((rho^2-dc^2+r^2)/rho)^2) % Eqn.19
% 
% if k>dc
%     A_H = dc^2*asin(c/(2*dc))-(c*dc/2)*cos(asin(c/(2*dc)))-...
%         r^2*asin(c/(2*r))+(c*r/2)*cos(asin(c/(2*r))) 
% else
%     A_H = pi*dc^2-r^2*asin(c/(2*r))+(c*r/2)*cos(asin(c/(2*r)))-...
%         dc^2*asin(c/(2*dc))+(c*dc/2)*cos(asin(c/(2*dc)))
% end % Eqn.18
% 
% A_C_i = double(e/2*((subs(xi,phi,phi_L_i)-subs(xi,phi,phi_F_i))*sin(theta_i)+...
%     (subs(yi,phi,phi_L_i)-subs(yi,phi,phi_F_i))*cos(theta_i)))% Eqn.32
% 
% A_O_i = A_s-A_t+A_H-A_TL_i+A_SL_i+A_TF_i-A_SF_i % Eqn.15
% 
% A_i = A_O_i-A_I_i+A_C_i

