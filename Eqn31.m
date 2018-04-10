function [ theta_s ] = Eqn31(rho,m,d,phi_s)
% Eqn.31

dxo = -rho*sin(phi_s) - d*m*sin(m*phi_s);
dyo = rho*cos(phi_s) + d*m*cos(m*phi_s);
theta_s = atan(dyo/dxo);

end

