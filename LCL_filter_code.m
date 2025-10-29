function x_opt = LCL_filter_code()

    h = 50;         % no of harmonics used
    f0 = 60;        % grid frequency
    f_sw = 5000;    % switching frequency

    % define fixed parameters
    Ls = 0.5e-3;
    Rs = 0.1;
    Rlg = 0.1;
    Rlc = 0.1;

    % Target/reference values
    THD_Vac_ref  = 0.06; 
    THD_Ig_ref = 0.02;
    THD_Ic_ref = 0.1;
    THD_Vdc_ref = 0.1;
    Vdc_ref = 700;
    Vac_peak = 310;
    P_ref = 100000;

    % define initial values of parameters to be estimated
    Lg0 = 1;
    C0 = 50;
    Lc0 = 5;
    Rd0 = 1;
    Cdc0 = 100;
    m0 = 0.8;
    x0 = [Lg0 C0 Lc0 Rd0 Cdc0 m0];

    % define limits 
    xL = [0.1, 0.1, 0.1, 0.001, 10, 0.5];
    xU = [30, 100, 30, 50, 2000, 0.99];

    scale = [1e-3, 1e-6, 1e-3, 1, 1e-6, 1];
    x0 = x0 .* scale;
    xL = xL .* scale;
    xU = xU .* scale;

    % LCL_filter_equation(x0);

    % opts = optimoptions('lsqnonlin','Algorithm', 'interior-point', 'Display','iter', 'FunctionTolerance', 1e-8, 'OptimalityTolerance', 1e-8, 'MaxFunctionEvaluations',5e3,'MaxIterations',1e3);
    % 
    % [x_opt, resnorm, residual, exitflag, output] = lsqnonlin(@objectiveFun, x0, xL, xU, [], [], [], [], @nonlinConstr, opts);

    opts = optimoptions('fmincon', ...
    'Algorithm','interior-point', ...   
    'Display','iter', ...
    'MaxFunctionEvaluations',5e3, ...
    'MaxIterations',1e3, ...
    'OptimalityTolerance',1e-8, ...
    'StepTolerance',1e-8, ...
    'ConstraintTolerance',1e-6);

    weight = [10; 10; 10; 10; 1; 1e8];
    obj = @(x) sum(weight .* (objectiveFun(x).^2));      
    [x_opt, fval] = fmincon(obj, x0, [], [], [], [], xL, xU, @nonlinConstr, opts);

    spect = [];
    P_spect = [];
    P_ref_init = P_ref;
    for i = -50000:5000:50000
        P_ref = P_ref_init + i;
        P_spect = [P_spect, P_ref];
        F = objectiveFun(x_opt);
        spect = [spect, F];
    end

    THD_Vac = spect(1,:) .*100;
    THD_Ig = spect(2,:) .*100;
    THD_Ic = spect(3,:) .*100;
    THD_Vdc = spect(4,:) .*100;
    % P_out = spect(5,:) + P_spect;

    subplot(2,2,1)
    plot(P_spect, THD_Vac)
    xlabel("Reference Power")
    ylabel("Capacitor voltage THD")
    subplot(2,2,2)
    plot(P_spect, THD_Ig)
    xlabel("Reference Power")
    ylabel("Grid side current THD")
    subplot(2,2,3)
    plot(P_spect, THD_Ic)
    xlabel("Reference Power")
    ylabel("Converter side current THD")
    subplot(2,2,4)
    plot(P_spect, THD_Vdc)
    xlabel("Reference Power")
    ylabel("DC source voltage THD")

    % [Iga, Ica, Vfa, Vdc, Idc] = LCL_filter_equation(x_opt);
    % disp(real(abs(Vdc(h+1))))
    % disp(real(abs(Idc(h+1))))
    % s = objectiveFun(x_opt);
    % disp(s)
    
    % x = x_opt ./ scale;
    % disp(x)

    function F = objectiveFun(x)

        [Iga, Ica, Vfa, Vdc, Idc] = LCL_filter_equation(x);
        
        rip1 = sum(abs(Vfa(h+3:end)).^2);
        rip2 = sum(abs(Iga(h+3:end)).^2);
        rip3 = sum(abs(Ica(h+3:end)).^2);
        rip4 = sum(abs(Vdc(h+2:end)).^2);
        
        Vac_rip = real(sqrt(rip1)/abs(Vfa(h+2)));
        Ig_rip = real(sqrt(rip2)/abs(Iga(h+2)));
        Ic_rip = real(sqrt(rip3)/abs(Ica(h+2)));
        Vdc_rip = real(sqrt(rip4)/abs(Vdc(h+1)));

        P = real(Idc(h+1) * Vdc(h+1));

        r1 = (Vac_rip);
        r2 = (Ig_rip);
        r3 = (Ic_rip);
        r4 = (Vdc_rip);
        r5 = (P - P_ref)/P_ref;
        r6 = (abs(Vdc(h+1)) - Vdc_ref) / Vdc_ref;
        F = [r1; r2; r3; r4; r5; r6];
        
    end


    function [c, ceq] = nonlinConstr(x)

        [Iga, Ica, Vfa, Vdc, Idc] = LCL_filter_equation(x);
        tol_vdc = 10; 
        % c = [abs(Vdc(h+1)) - (Vdc_ref + tol_vdc); (Vdc_ref - tol_vdc) - abs(Vdc(h+1))]; 
        % ceq = [];
        c = [];
        ceq = [(abs(Vdc(h+1)) - Vdc_ref) / Vdc_ref];
        % opts.ConstraintTolerance = 1e-8;

    end


    function [Iga, Ica, Vfa, Vdc, Idc] = LCL_filter_equation(x)
    
        [Lg, C, Lc, Rd, Cdc, m] = struct('s', num2cell(x)).s;
        
        % fetching EHD matrices
        [D, Ga, Za, Va, Idc] = EHD_matrices(f0, f_sw, m, 2048, h, P_ref, Vdc_ref, Vac_peak);
       
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
        
        U = [Va; Idc];
        
        % steady state solution
        Xss = -(A \ (B*U));
        
        % voltage and current harmonics vectors
        Iga = Xss(1:2*h+1);
        Ica = Xss(2*h+1+1:2*(2*h+1));
        Vfa = Xss(2*(2*h+1)+1:3*(2*h+1));
        Vdc = Xss(3*(2*h+1)+1:4*(2*h+1));

    end

end