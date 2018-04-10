function [ A_C ] = Eqn32(rho,m,d,dc,e,phi_L,phi_F,theta_i)
% Eqn.32

part1 = func_xi(rho,m,d,dc,phi_L);
part2 = func_xi(rho,m,d,dc,phi_F);
part3 = func_yi(rho,m,d,dc,phi_L);
part4 = func_yi(rho,m,d,dc,phi_F);

A_C = e/2*((part1-part2)*sin(theta_i)+(part3-part4)*cos(theta_i));
% Eqn.32

end

