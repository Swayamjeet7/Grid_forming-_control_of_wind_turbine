function [P, Q, L_1, L_2, C_f, R_f] = lcl_param()

P = 5000;
Q = 250;
f_g = 60;
f_sw = 15000;
V_dc = 400;
V_g = 120 * sqrt(3);
L = 2500;

x = 0.05; % maximum power factor variation seen by the grid
k = 0.2; % attenuation factor
m = 0.5; % modulation factor

% base values
Z_b = (V_g^2)/P;
C_b = 1/(2*pi*f_g*Z_b);
C_f = x*C_b;

% current ripple
I_m = P*sqrt(2)/(sqrt(3)*V_g);
dI_m = 0.1*I_m;

% inductor parameters
L_1 = (2*V_dc*m*(1-m))/(3*f_sw*dI_m);
L_2 = (1+sqrt(1/k^2))/(C_f*(2*pi*f_sw)^2);

% resonance frequency
w_res = sqrt((L_1+L_2+L)/(L_1*L_2*C_f));
f_res = w_res/(2*pi);

% damping resistor
R_f = 1000/(3*w_res*C_f);

% delta connection
R_f = 3*R_f;
C_f = C_f/3;

% check for resonance
if (10*f_g) < f_res && f_res < (0.5*f_sw) 
    disp('all good')
else
    disp('parameters should be rechosen !')
end

end
