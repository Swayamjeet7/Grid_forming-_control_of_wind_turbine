
% function f = LCL_filter_code()
Lg = 1e-3;
C = 50e-6;
Lc = 3e-3;
Rd = 0.1;
Cdc = 500e-6;
m = 0.8;

x0 = [Lg C Lc Rd Cdc m];
xL = [0.1e-3, 0.1e-6, 0.1e-3, 0.001, 10e-6, 0.5];
xU = [30e-3, 100e-6, 30e-3, 50, 2000e-6, 0.99];

% [Iga, Ica, Vfa, Vdc, h] = LCL_filter_equation(x0);
% disp(Ica)

options = optimoptions('fmincon',...
    'Algorithm','interior-point',...
    'Display','iter-detailed',...
    'StepTolerance',1e-9,...
    'ConstraintTolerance',1e-3,...
    'OptimalityTolerance',1e-9,...
    'HonorBounds',true,...
    'SpecifyObjectiveGradient',false,...
    'SpecifyConstraintGradient',false,...
    'ScaleProblem',true,...
    'SubproblemAlgorithm','cg',...
    'UseParallel',false,...
    'InitTrustRegionRadius',0.1,...
    'EnableFeasibilityMode',true, ...
    'MaxFunctionEvaluations', 1e3);  



[x_opt, fval] = fmincon(@objectiveFun, x0, [], [], [], [], xL, xU, @nonlinConstr, options);
disp(x_opt)


function F = objectiveFun(x)
[Iga, Ica, Vfa, Vdc, h] = LCL_filter_equation(x);
R = 0.1;
Pcu = sum(R .* abs(Iga).^2);
F = Pcu;
end

function [c, ceq] = nonlinConstr(x)
[Iga, Ica, Vfa, Vdc, h] = LCL_filter_equation(x);
% THD = sqrt(sum(abs(Ia(2:end)).^2)) / abs(Ia(1));
% L = x(1); C = x(2);
% c1 = 1e-6 - L;            
% c2 = 1e-8 - C;

I1 = Iga(h+1+1);      
num = 0;
for k = 2:h
    num = num + abs(Iga(h+1+k))^2;   
end
THD = sqrt(num)/abs(I1);

c3 = THD - 0.03;    
c = [c3];
ceq = [];
% disp(THD)
end

function [Iga, Ica, Vfa, Vdc, h] = LCL_filter_equation(x)

Ls = 0.5e-3;
Rs = 0.1;
Rlg = 0.1;
Rlc = 0.1;

Lg = x(1);
C = x(2);
Lc = x(3);
Rd = x(4);
Cdc = x(5);
m = x(3);

f_sw = 5000;
f0 = 60;
h = 2;

% fetching EHD matrices
[D, Ga, Za, Va, Il] = EHD_matrices(f0, f_sw, m, 2048, h);

% making constant matrices
O = zeros(2*h+1,2*h+1);
Id = eye(2*h+1);
a = -(Rs+Rlg+Rd)/(Ls+Lg);
b = -(Rd+Rlc)/Lc;

% state space matrices
A1 = [a*Id-D (Rd/(Ls+Lg))*Id (-1/(Ls+Lg))*Id O];
A2 = [(Rd/Lc)*Id b*Id-D (1/Lc)*Id -Ga/Lc];
A3 = [(1/C)*Id (-1/C)*Id -D O];
A4 = [O 3*Za/Cdc O -D];

A = [A1; A2; A3; A4];
B = [(1/(Ls+Lg))*Id O; O O; O O; O (1/Cdc)*Id];

U = [Va; Il];

% steady state solution
Xss = -inv(A)*B*U;

% voltage and current harmonics vectors
Iga = Xss(1:2*h+1);
Ica = Xss(2*h+1+1:2*(2*h+1));
Vfa = Xss(2*(2*h+1)+1:3*(2*h+1));
Vdc = Xss(3*(2*h+1)+1:4*(2*h+1));

% ia = inverseFFT(Ia_h, f0, h);
% vdc = inverseFFT(Vdc_h, f0, h);

end


function X = inverseFFT(X_h, f, h)

Nfft = 1+1/(f*0.0001);

FFTbins = zeros(Nfft,1);
FFTbins(1) = X_h(h+1); %dc term
for k = 1:h
    % for positive k
    FFTbins(k+1) = X_h(h+1+k);          
    % conjugate for negative k
    FFTbins(Nfft - k + 1) = conj(X_h(h+1-k)); 
end
X = real(ifft(FFTbins) * Nfft);

end
