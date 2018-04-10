function [ beta_L_i ] = Eqn25(rho,m,e,theta_i)
% Eqn.25

beta_L_i = asin((m*e*sin((m-1)/m*theta_i+pi/m))/...
        sqrt(m^2*e^2+rho^2-2*m*e*rho*cos((m-1)/m*theta_i+pi/m))); % Eqn.25

end

