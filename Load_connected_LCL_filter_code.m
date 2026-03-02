function [] = Load_connected_LCL_filter_code()

    % define values
    h = 5;          % no of harmonics used
    f0 = 60;        % grid frequency
    f_sw = 5040;    % switching frequency

    Vdc_ref = 30;   % source DC voltage
    Vac_peak = 0;   
    P_ref = 50;     % required load power
    Rload = 2;      % load resistance
    Ls = 0;         % load inductance
    Rloss = 0;      % power converter loss

    % define fixed filter parameters
    Rlg = 0.1;
    Rlc = 0.1;

    % define initial values of filter parameters to be estimated
    Lg0 = 1;
    C0 = 30;
    Lc0 = 2;
    Rd0 = 10;
    m0 = 0.8;
    theta0 = 0;
    x0 = [Lg0 C0 Lc0 Rd0 m0 theta0];

    % define limits 
    xL = [0.1, 0.1, 0.1, 0.001, 0.1, -pi];
    xU = [30, 200, 30, 50, 0.99, pi];

    % scale for optimization
    scale = [1e-3, 1e-6, 1e-3, 1, 1, 1];
    x0 = x0 .* scale;
    xL = xL .* scale;
    xU = xU .* scale;

    % options for fmincon
    opts = optimoptions('fmincon', ...
    'Algorithm','interior-point', ...
    'Display','none', ...
    'MaxFunctionEvaluations',1e5, ...
    'MaxIterations',5000, ...
    'OptimalityTolerance',1e-6, ...
    'StepTolerance',1e-6, ...
    'ConstraintTolerance',1e-3, ...
    'FiniteDifferenceType','forward', ...
    'FiniteDifferenceStepSize',1e-5, ...
    'TypicalX', [1, 1, 1, 10, 1, 0.1], ...
    'SpecifyObjectiveGradient',false);
  
    % calling fmincon
    [x_opt, fval, exitflag, output] = fmincon(@objectiveFun, x0, [], [], [], [], xL, xU, @nonlinConstr, opts);

    % diagnose()
    out()

    % diagnostic prints
    function [] = diagnose()
        
        disp('--- DIAGNOSTICS ---');
        disp(['exitflag: ', num2str(exitflag)]);
        if exist('output','var')
            if isfield(output,'constrviolation'); disp(['constrviolation: ', num2str(output.constrviolation)]); end
            if isfield(output,'message'); disp(['message: ', output.message]); end
        end

        % evaluate constraints at x_opt
        [c_opt, ceq_opt] = nonlinConstr(x_opt);
        disp('nonlinear inequality (c) at x_opt:'); disp(c_opt(:).');
        disp('nonlinear equality (ceq) at x_opt:'); disp(ceq_opt(:).');

        % evaluate steady state
        [Iga_opt, Ica_opt, Vfa_opt] = LCL_filter_equation(x_opt);
        P = 2*((abs(Iga_opt(h+2))^2) * Rload );
        fprintf("\nsteady state power is %.f watt", P)
   
    end
     
    % power and thd calculations
    function [Ig_thd, Ic_thd, Vf_thd, P_actual] = calc(x)
        
        [Iga, Ica, Vfa] = LCL_filter_equation(x);
        
        thd1 = sum(abs(Vfa(h+3:end)).^2);
        thd2 = sum(abs(Iga(h+3:end)).^2);
        thd3 = sum(abs(Ica(h+3:end)).^2);
        
        eps_small = 1e-12;
        Vf_thd = real(sqrt(thd1) / max(eps_small, abs(Vfa(h+2))));
        Ig_thd = real(sqrt(thd2) / max(eps_small, abs(Iga(h+2))));
        Ic_thd = real(sqrt(thd3) / max(eps_small, abs(Ica(h+2))));
        
        P_actual = 2*((abs(Iga(h+2))^2) * Rload );
        
    end

    % objective function for fmincon
    function F = objectiveFun(x)
        
        [Ig_thd, Ic_thd, Vf_thd, P_actual] = calc(x);

        P_error = (P_actual - P_ref) / P_ref;

        Fv = [Ig_thd; Ic_thd; Vf_thd; P_error];

        weight = [10; 10; 1; 1e4]; 
        F = sum(weight .* (Fv.^2));
    
    end

    % constraints for fmincon
    function [c, ceq] = nonlinConstr(x)

        Lg = x(1);
        C  = x(2);
        Lc = x(3);
        Rd = x(4);
        m  = x(5);
        
        % Base Values for limits
        V_rms_nom = Vdc_ref / sqrt(2); 
        Z_base = (V_rms_nom^2) / P_ref;
        C_base = 1 / (2 * pi * f0 * Z_base);
        L_base = Z_base / (2 * pi * f0);

        % resonance frequency constraints
        f_res = (1 / (2*pi)) * sqrt((Lg + Lc) / (Lg * Lc * C));
        f_res_min = 10 * f0;
        f_res_max = 0.5 * f_sw;
        
        % normalizing the constraints
        c1 = (f_res_min - f_res) / f_res_min;  % c1 <= 0 means f_res >= f_res_min
        c2 = (f_res - f_res_max) / f_res_max;  % c2 <= 0 means f_res <= f_res_max
        
        % reactive power limit (max 10% of base capacitance)
        C_max = 0.1 * C_base;
        c3 = (C - C_max) / C_max;        % c3 <= 0 means C <= C_max  
        
        % total inductance limit (max 20% voltage drop)
        L_max = 0.2 * L_base;
        c4 = ((Lg + Lc) - L_max) / L_max;  
        
        % damping resistor limit
        w_res = 2 * pi * f_res;
        Rd_opt = 1 / (3 * w_res * C);
        % bound Rd to not exceed 3x the theoretical value
        c5 = (Rd - (3 * Rd_opt)) / (3 * Rd_opt); 

        % power matching
        [~, Iga, ~] = LCL_filter_equation(x);
        P_actual = 2 * ((abs(Iga(h+2))^2) * Rload);
        
        c = [c1, c2, c3, c4, c5]; 
        ceq = (P_actual - P_ref) / P_ref;

    end

    % steady state equations in EHD
    function [Iga, Ica, Vfa] = LCL_filter_equation(x)
    
        % unpacking parameters
        Lg    = x(1);
        C     = x(2);
        Lc    = x(3);
        Rd    = x(4);
        m     = x(5);
        theta = x(6);

        % fetching EHD matrices
        [D, Ga, ~, ~, Vdc] = EHD_matrices(f0, f_sw, m, theta, 2048, h, Vac_peak, Vdc_ref);
       
        % making constant matrices
        O = zeros(2*h+1,2*h+1);
        Id = eye(2*h+1);
        a = -(Rlg+Rd+Rload)/(Ls+Lg);
        b = -(Rd+Rlc+Rloss)/Lc;
        
        % state space matrices
        A1 = [a*Id-D (Rd/(Ls+Lg))*Id (1/(Ls+Lg))*Id];
        A2 = [(Rd/Lc)*Id b*Id-D -(1/Lc)*Id];
        A3 = [-(1/C)*Id (1/C)*Id -D];
        
        A = [A1; A2; A3];
        B = [O; (1/(Lc))*Ga; O];
        
        U = [Vdc];
       
        % check if A is singular, rcond checks reciprocal condition number
        if rcond(A) < 1e-12 
            Iga = NaN * ones(2*h+1, 1);
            Ica = NaN * ones(2*h+1, 1);
            Vfa = NaN * ones(2*h+1, 1);
            return
        end
        
        % steady state solution
        Xss = -(A \ (B*U));
        
        % voltage and current harmonics vectors
        Iga = Xss(1:2*h+1);
        Ica = Xss(2*h+1+1:2*(2*h+1));
        Vfa = Xss(2*(2*h+1)+1:3*(2*h+1));

    end

    % final outputs and plots
    function [] = out()
        [Iga, Ica, Vfa] = LCL_filter_equation(x_opt);

        [Ig_thd, Ic_thd, Vf_thd, P_actual] = calc(x_opt);

        % optimized paramter values
        ansp = round((x_opt ./ scale), 3);
        fprintf('\nFilter parameters for %.f V supply and %.f W reference power are \n', Vdc_ref, P_actual);
        disp('----------------------------------------------------------------------------');
        disp(['converter side inductance:   ', num2str(ansp(3)), ' mH']);
        disp(['filter capacitance:          ', num2str(ansp(2)), ' μF']);
        disp(['load side inductance:        ', num2str(ansp(1)), ' mH']);
        disp(['damping resitance:           ', num2str(ansp(4)), ' Ω']);
        disp(['modulation index:            ', num2str(ansp(5))]);
        disp(['theta:                       ', num2str(ansp(6)), ' rad']);
        disp(' ');
        
        % THD values at reference power        
        anst = [Ig_thd, Ic_thd, Vf_thd];
        ansthd = round((anst*100), 2);
        fprintf("\nTHD values at reference power are\n");
        disp('-----------------------------------');
        disp(['Load side current THD:        ', num2str(ansthd(1)), ' %']);
        disp(['Converter side current THD:   ', num2str(ansthd(2)), ' %']);
        disp(['Filter capacitor voltage THD: ', num2str(ansthd(3)), ' %']);
        disp(' ');

        % % graphs for performance indices w.r.t. reference power
        % spect = [];
        % P_spect = [];
        % P_ref_init = P_ref;
        % for i = -(P_ref_init * 0.8):10:(P_ref_init * 0.8)
        %     P_ref = P_ref_init + i;
        %     x1 = [P_ref/Vdc_ref];
        %     [x1_opt] = fmincon(@obj2, x1, [], [], [], [], 0, [], [], opts2);
        %     x_opt(7) = x1_opt;
        %     P_spect = [P_spect, P_ref];
        %     F = objectiveFun(x_opt);
        %     spect = [spect, F];
        % end
        % 
        % THD_Vac = spect(1,:) .*100;
        % THD_Ig = spect(2,:) .*100;
        % THD_Ic = spect(3,:) .*100;
        % THD_Vdc = spect(4,:) .*100;
        % 
        % subplot(2,2,1)
        % plot(P_spect, THD_Vac)
        % xlabel("Reference Power")
        % ylabel("Capacitor voltage THD")
        % subplot(2,2,4)
        % plot(P_spect, THD_Ig)
        % xlabel("Reference Power")
        % ylabel("Grid side current THD")
        % subplot(2,2,3)
        % plot(P_spect, THD_Ic)
        % xlabel("Reference Power")
        % ylabel("Converter side current THD")
        % subplot(2,2,2)
        % plot(P_spect, THD_Vdc)
        % xlabel("Reference Power")
        % ylabel("DC source voltage THD")

        Idc = P_actual / Vdc_ref;
        Ig_thd_spect = [];
        Ic_thd_spect = [];
        Vf_thd_spect = [];
        P_spect = [];

        volt = [];
        volt_init = Vdc_ref;
        for i = -(volt_init * 0.33):2:(volt_init * 0.33)
            Vdc_ref = volt_init + i;
            volt = [volt, Vdc_ref];
            [Ig_thd, Ic_thd, Vf_thd, P_actual] = calc(x_opt);
            Ig_thd_spect = [Ig_thd_spect, Ig_thd];
            Ic_thd_spect = [Ic_thd_spect, Ic_thd];
            Vf_thd_spect = [Vf_thd_spect, Vf_thd];
            P_spect = [P_spect, P_actual];
        end
        
        % plot(volt, P_spect);

    end

end