% Try to calculate the total area A_i.
% But not sure whether all of those equations are correct.
% Use Func instead of script.
% mengtang li
% Nov 26 2017

function [angle,A_I_i,A_C_i,A_O_i,A_i] = func_test1(rho, d, m, dc, r)
e = d; % ???
% Primary Input: rho, m, d
% Secondary Input: dc, e, r

% Constraint1: md < rho % Eqn.3 or Eqn.47
% if (m*d > rho)
%     error(message('Constraint1: md < rho % Eqn.3 or Eqn.47'));
% end
%
% % Constraint2: dc < R_min % Eqn.6+7 or Eqn.48
% phi_m = Eqn7(rho,m,d) % Eqn.7
% R_min = Eqn6(rho,m,d,phi_m) % Eqn.6
% if (dc > R_min)
%     error(message('Constraint2: dc < R_min % Eqn.6+7 or Eqn.48'));
% end
%
% % Constraint3: e < dc/2 % Eqn.50 or Eqn.51
% if (e > dc/2)
%     error(message('% Constraint3: e < dc/2 % Eqn.50 or Eqn.51'));
% end
%
% % Constraint4: r_min < r < r_max
% r_min = rho+2*e-dc % Eqn.8
% r_max = rho+dc % Eqn.9
% if (r < r_min || r > r_max)
%     error(message('% Constraint4: r_min < r < r_max'));
% end

angle = 0:0.01:20; %0:0.01:2*pi*m/(m-1)
n = size(angle,2);
A_I_i = zeros(1,n);
A_C_i = zeros(1,n);
A_O_i = zeros(1,n);
A_i = zeros(1,n);
A_s = pi*r^2/m; % Eqn.16
A_t = pi*dc^2; % Eqn.17

for i = 1:1:n
    theta_i = angle(i);
    beta_F_i = Eqn26(rho,m,e,theta_i); % Eqn.26
    beta_L_i = Eqn25(rho,m,e,theta_i); % Eqn.25
    A_SF_i = beta_F_i*dc^2/2; % Eqn.24
    A_TF_i = rho*dc*sin(beta_F_i)/2; % Eqn.23
    A_SL_i = beta_L_i*dc^2/2; % Eqn.22
    A_TL_i = rho*dc*sin(beta_L_i)/2; % Eqn.21
    k = (-rho^2+2*dc*rho-dc^2+r^2)/(2*rho); % Eqn.20
    c = sqrt(4*r^2-((rho^2-dc^2+r^2)/rho)^2); % Eqn.19
    
    if k>dc
        A_H = dc^2*asin(c/(2*dc))-(c*dc/2)*cos(asin(c/(2*dc)))-...
            r^2*asin(c/(2*r))+(c*r/2)*cos(asin(c/(2*r)));
    else
        A_H = pi*dc^2-r^2*asin(c/(2*r))+(c*r/2)*cos(asin(c/(2*r)))-...
            dc^2*asin(c/(2*dc))+(c*dc/2)*cos(asin(c/(2*dc)));
    end % Eqn.18
    
    phi_s = acos(-m*e/rho)/(m-1); % Eqn.30
    theta_s = Eqn31(rho,m,d,phi_s); % Eqn.31
    phi_L_i = -phi_s/(pi+theta_s)*(theta_i+theta_s)+phi_s; % Eqn.28
    phi_F_i = -phi_s/(pi+theta_s)*(theta_i+theta_s+2*pi)+phi_s; % Eqn.29
    
    int_element = @(phi)(d*cos(m*phi)+rho*cos(phi)-(dc*(rho*cos(phi)+ ...
        d*m*cos(m*phi)))./((rho*cos(phi)+d*m*cos(m*phi)).^2+(rho*sin(phi)+ ...
        d*m*sin(m*phi)).^2).^(1/2)).*(rho*cos(phi)-(dc*(rho*cos(phi)+ ...
        d*m^2*cos(m*phi)))./((rho*cos(phi)+d*m*cos(m*phi)).^2+(rho*sin(phi)+ ...
        d*m*sin(m*phi)).^2).^(1/2)+d*m*cos(m*phi)+(dc*(2*(rho*cos(phi)+ ...
        d*m^2*cos(m*phi)).*(rho*sin(phi)+d*m*sin(m*phi))-2*(rho*sin(phi)+ ...
        d*m^2*sin(m*phi)).*(rho*cos(phi)+d*m*cos(m*phi))).*(rho*sin(phi)+ ...
        d*m*sin(m*phi)))./(2*((rho*cos(phi)+d*m*cos(m*phi)).^2+(rho*sin(phi)+ ...
        d*m*sin(m*phi)).^2).^(3/2)))+(d*sin(m*phi)+rho*sin(phi)- ...
        (dc*(rho*sin(phi)+d*m*sin(m*phi)))./((rho*cos(phi)+d*m*cos(m*phi)).^2+ ...
        (rho*sin(phi)+d*m*sin(m*phi)).^2).^(1/2)).*(rho*sin(phi)- ...
        (dc*(rho*sin(phi)+d*m^2*sin(m*phi)))./((rho*cos(phi)+ ...
        d*m*cos(m*phi)).^2+(rho*sin(phi)+d*m*sin(m*phi)).^2).^(1/2)+ ...
        d*m*sin(m*phi)-(dc*(2*(rho*cos(phi)+d*m^2*cos(m*phi)).*(rho*sin(phi)+ ...
        d*m*sin(m*phi))-2*(rho*sin(phi)+d*m^2*sin(m*phi)).*(rho*cos(phi)+ ...
        d*m*cos(m*phi))).*(rho*cos(phi)+d*m*cos(m*phi)))./(2*((rho*cos(phi)+ ...
        d*m*cos(m*phi)).^2+(rho*sin(phi)+d*m*sin(m*phi)).^2).^(3/2)));
    
    
    A_I_i(i) = 0.5*integral(int_element,phi_L_i,phi_F_i); % Eqn.27
    A_C_i(i) = Eqn32(rho,m,d,dc,e,phi_L_i,phi_F_i,theta_i); % Eqn.32
    A_O_i(i) = A_s-A_t+A_H-A_TL_i+A_SL_i+A_TF_i-A_SF_i; % Eqn.15
    
    A_i(i) = A_O_i(i)+A_I_i(i)-A_C_i(i); % Eqn.14
end

clf;
figure(1);
plot(angle,A_O_i,'b','LineWidth',2);
hold on; grid minor; grid on;
plot(angle,A_I_i,'r','LineWidth',2);
plot(angle,A_C_i,'y','LineWidth',2);
legend('AOi', 'AIi', 'ACi');
ax = gca; % current axis handle
set(ax,'XTick',[0 pi/2 pi 3*pi/2 2*pi 5*pi/2 3*pi 7*pi/2 4*pi 9*pi/2 5*pi 11*pi/2 6*pi 13*pi/2 7*pi 15*pi/2 8*pi]);
set(ax,'XTickLabel',{'0','90^o','180^o','270^o','360^o', '450^o',...
    '540^o', '630^o', '720^o', '810^o', '900^o', '990^o', '1080^o',...
    '1170^o', '1260^o', '1350^o', '1440^o'});

figure(2);
plot(angle,A_i,'b','LineWidth',2);
grid minor; grid on;
xlim([0 5*pi/2]);
bx = gca; % current axis handle
set(bx,'XTick',[0 pi/2 pi 3*pi/2 2*pi 5*pi/2]);
set(bx,'XTickLabel',{'0','90^o','180^o','270^o','360^o', '450^o',...
    '540^o'});

Area = (max(A_i)-min(A_i))*(m-1)
end