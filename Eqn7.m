function [ phi_m ] = Eqn7(rho,m,d)
% Eqn.7

phi_m = 1/(m-1)*...
    acos((rho^2+m^2*d^2)/(m*d*rho)-3*(rho^2+m^3*d^2)/(m*d*rho*(m+1)));

end

