%% Goal: Function containing the ODEs of model extension 6.1
% function containing the ODEs of model extension 6.1, the parameters and an IF
% statement to change the formaldehyde concentration at a time point during the
% simulation (not used in this simulation)
%% ODEs of model extension 6.1
function dydt = model_e6_2(t,y,Fval,t1on,t1off,p)

% dydt(1) = d/dt⋅m_frmR
% dydt(2) = d/dt⋅frmR
% dydt(3) = d/dt⋅FfrmR
% dydt(4) = d/dt m_LacI
% dydt(5) = d/dt⋅LacI
% dydt(6) = d/dt m_sok
% dydt(7) = d/dt m_hoksok
% dydt(8) = d/dt m_hok
% dydt(9) = d/dt⋅m_x
% dydt(10)= d/dt⋅X

%a1 = p(1);
a2 = p(2);
a3 = p(3);
a4 = p(4);
a5 = p(5);
b1 = p(6);
b2 = p(7);
b3 = p(8);
b4 = p(9);
b5 = p(10);
b6 = p(11);
b7 = p(12);
g1 = p(13);
g2 = p(14);
k1 = p(15);
k2 = p(16);
k3 = p(17);
k4 = p(18);
k5 = p(19);
k6 = p(20);
n1 = p(21);
n2 = p(22);
% new parameters
a6 = p(23);
a7 = p(24);
b8 = p(25);
b9 = p(26);
g3 = p(27);
k7 = p(28);
n3 = p(29);

% IF statement to turn inputs on or off

if t > t1on && t <= t1off
    F = Fval;
else
    F = 9*(10^-6);
end

% Model equations
dydt = zeros(size(y));

dydt(1) = 3600/(10^-9)*(a6 + a7*(k7^n3)/((k7^n3)+((y(10)*(10^-9))^n3)) - b1*y(1)*(10^-9));
dydt(2) = 3600/(10^-9)*(g1*y(1)*(10^-9)-b2*y(2)*(10^-9)+k1*y(3)*(10^-9)-k2*F*y(2)*(10^-9));
dydt(3) = 3600/(10^-9)*(k2*F*y(2)*(10^-9)-k1*y(3)*(10^-9)- b2*y(3)*(10^-9));
dydt(4) = 3600/(10^-9)*(a2 + a3*(k3^n1)/((k3^n1)+((y(2)*(10^-9))^n1))-b3*y(4)*(10^-9));
dydt(5) = 3600/(10^-9)*(g2*y(4)*(10^-9)-b4*y(5)*(10^-9));
dydt(6) = 3600/(10^-9)*(a2 + a3*(k3^n1)/((k3^n1)+((y(2)*(10^-9))^n1))-k4*y(6)*(10^-9)*y(8)*(10^-9)+k5*y(7)*(10^-9)-b5*y(6)*(10^-9));
dydt(7) = 3600/(10^-9)*(k4*y(6)*(10^-9)*y(8)*(10^-9)-k5*y(7)*(10^-9)-b6*y(7)*(10^-9));
dydt(8) = 3600/(10^-9)*(a4 + a5*(k6^n2)/((k6^n2)+((y(5)*(10^-9))^n2))-k4*y(6)*(10^-9)*y(8)*(10^-9)+k5*y(7)*(10^-9)-b7*y(8)*(10^-9));
dydt(9) = 3600/(10^-9)*(a2 + a3*(k3^n1)/((k3^n1)+((y(2)*(10^-9))^n1))-b8*y(9)*(10^-9));
dydt(10)= 3600/(10^-9)*(g3*y(9)*(10^-9)-b9*y(10)*(10^-9));

end