%% Goal: Function containing the ODEs of the toy model
% function containing the ODEs of the toy model, the parameters and an IF
% statement to change the IPTG concentration at a time point during the
% simulation (not used in this simulation)
%% ODEs of the toy model
function dydt = hoksok_toymodel(t,y,Ival,t1on,t1off,p)

% dydt(1) = d/dt m_sok
% dydt(2) = d/dt m_hoksok
% dydt(3) = d/dt m_hok

% parameters 
a4 = p(1);
a5 = p(2);
a6 = p(3);
b5 = p(4);
b6 = p(5);
b7 = p(6);
k5 = p(7);
k6 = p(8);
k7 = p(9);
n3 = p(10);

% IF statement to turn IPTG input on or off
if t > t1on && t <= t1off
    I = Ival;
else
    I = 0;
end

% Model equations
dydt = zeros(size(y));

dydt(1) = 3600/(10^-9)*(a6 - k5*y(1)*(10^-9)*y(3)*(10^-9)+k6*y(2)*(10^-9)-b5*y(1)*(10^-9)); %sok
dydt(2) = 3600/(10^-9)*(k5*y(1)*(10^-9)*y(3)*(10^-9)-k6*y(2)*(10^-9)-b6*y(2)*(10^-9)); %hoksok
dydt(3) = 3600/(10^-9)*(a4 + a5*(I^n3)/(k7^n3+I^n3)-k5*y(1)*(10^-9)*y(3)*(10^-9)+k6*y(2)*(10^-9)-b7*y(3)*(10^-9)); %hok

end