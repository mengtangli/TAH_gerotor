function [ beta_F_i ] = Eqn26(rho,m,e,theta_i)
% Eqn.26

beta_F_i = asin((m*e*sin((m-1)/m*theta_i-pi/m))/...
        sqrt(m^2*e^2+rho^2-2*m*e*rho*cos((m-1)/m*theta_i-pi/m))); 

end

