function [ xi ] = func_xi(rho,m,d,dc,phi)
% Eqn for xi

xi = d*cos(m*phi)+rho*cos(phi)-(dc*(rho*cos(phi)+d*m*cos(m*phi)))/ ...
    ((rho*cos(phi)+d*m*cos(m*phi))^2+(rho*sin(phi)+d*m*sin(m*phi))^2)^(1/2);

end

