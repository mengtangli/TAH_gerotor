function [ R ] = Eqn6(rho,m,d,phi)
% Eqn.6

R_num = (rho^2+m^2*d^2+2*m*d*rho*cos((m-1)*phi))^1.5;
R_den = (rho^2+m^3*d^2+m*d*rho*(m+1)*cos((m-1)*phi));
R = R_num/R_den; 

end

