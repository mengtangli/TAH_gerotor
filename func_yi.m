function [ yi ] = func_yi(rho,m,d,dc,phi)
% Eqn for yi

yi = d*sin(m*phi)+rho*sin(phi)-(dc*(rho*sin(phi)+d*m*sin(m*phi)))/ ...
    ((rho*cos(phi)+d*m*cos(m*phi))^2+(rho*sin(phi)+d*m*sin(m*phi))^2)^(1/2);

end

