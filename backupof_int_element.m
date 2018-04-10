% my hand calculation of Eqn.27
int_element = @(phi)(d*cos(m*phi)+rho*cos(phi)-(dc*(rho*cos(phi)+ ...
    d*m*cos(m*phi)))./((rho*cos(phi)+d*m*cos(m*phi)).^2+(rho*sin(phi)+ ...
    d*m*sin(m*phi)).^2).^(1/2)).*(rho*cos(phi)-(dc*((rho*cos(phi)+ ...
    d*m.^2*cos(m*phi)).*((rho*cos(phi)+d*m*cos(m*phi)).^2+(rho*sin(phi)+ ...
    d*m*sin(m*phi)).^2)-(rho*sin(phi)+d*m*sin(m*phi)).*(rho*cos(phi)- ...
    rho*sin(phi)+d*m*cos(m*phi)-d*m*sin(m*phi))))./((rho*cos(phi)+ ...
    d*m*cos(m*phi)).^2+(rho*sin(phi)+d*m*sin(m*phi)).^2).^(3/2)+ ...
    d*m*cos(m*phi))+(d*sin(m*phi)+rho*sin(phi)-(dc*(rho*sin(phi)+ ...
    d*m*sin(m*phi)))./((rho*cos(phi)+d*m*cos(m*phi)).^2+(rho*sin(phi)+ ...
    d*m*sin(m*phi)).^2).^(1/2)).*(rho*sin(phi)-(dc*((rho*sin(phi)+ ...
    d*m.^2*sin(m*phi)).*((rho*cos(phi)+d*m*cos(m*phi)).^2+(rho*sin(phi)+ ...
    d*m*sin(m*phi)).^2)+(rho*cos(phi)+d*m*cos(m*phi)).*(rho*cos(phi)- ...
    rho*sin(phi)+d*m*cos(m*phi)-d*m*sin(m*phi))))./((rho*cos(phi)+ ...
    d*m*cos(m*phi)).^2+(rho*sin(phi)+d*m*sin(m*phi)).^2).^(3/2)+ ...
    d*m*sin(m*phi)); 

% diff command calculation of Eqn.27
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
    (rho*sin(phi)+d*m*sin(m*phi)).^2).^(1/2))*(rho*sin(phi)- ...
    (dc*(rho*sin(phi)+d*m^2*sin(m*phi)))./((rho*cos(phi)+ ...
    d*m*cos(m*phi)).^2+(rho*sin(phi)+d*m*sin(m*phi)).^2).^(1/2)+ ...
    d*m*sin(m*phi)-(dc*(2*(rho*cos(phi)+d*m^2*cos(m*phi)).*(rho*sin(phi)+ ...
    d*m*sin(m*phi))-2*(rho*sin(phi)+d*m^2*sin(m*phi)).*(rho*cos(phi)+ ...
    d*m*cos(m*phi))).*(rho*cos(phi)+d*m*cos(m*phi)))./(2*((rho*cos(phi)+ ...
    d*m*cos(m*phi)).^2+(rho*sin(phi)+d*m*sin(m*phi)).^2).^(3/2)));