function [] = LCL_filter_code()

    h = 50;         % no of harmonics used
    f0 = 60;        % grid frequency
    f_sw = 5000;    % switching frequency

    % define fixed parameters
    Ls = 0.5e-3;
    Rs = 0.1;
    Rlg = 0.1;
    Rlc = 0.1;

    % Target/reference values
    Vdc_ref = 700;
    Vac_peak = 310;
    P_ref = 5000;

    % define initial values of parameters to be estimated
    Lg0 = 1;
    C0 = 30;
    Lc0 = 2;
    Rd0 = 10;
    Cdc0 = 50;
    m0 = 0.8;
    theta0 = 0;
    Idc0 = P_ref/Vdc_ref;
    x0 = [Lg0 C0 Lc0 Rd0 Cdc0 m0 theta0 Idc0];

    % define limits 
    xL = [0.1, 0.1, 0.1, 0.001, 10, 0.5, -pi, 0];
    xU = [30, 200, 30, 50, 2000, 0.99, pi, 5*Idc0];

    % scale for optimization
    scale = [1e-3, 1e-6, 1e-3, 1, 1e-6, 1, 1, 1];
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
    'ConstraintTolerance',1e-1, ...
    'FiniteDifferenceType','forward', ...
    'FiniteDifferenceStepSize',1e-5, ...
    'TypicalX', [1, 1, 1, 10, 1, 1, 0.1, 10], ...
    'SpecifyObjectiveGradient',false);

    % objective function provided to fmincon
    obj = @(x) obj_wrapper(x);

    function J = obj_wrapper(x)
        
        weight = [1; 10; 10; 1; 1e6];
        penalty_vdc = 1e6;

        Fvec = objectiveFun(x);
        [~, ~, ~, Vdc] = LCL_filter_equation(x);
        v_err = (abs(Vdc(h+1)) - Vdc_ref) / Vdc_ref;
        J = sum(weight .* (Fvec.^2)) + penalty_vdc * (v_err^2);

    end
    
    % calling fmincon
    [x_opt, fval, exitflag, output] = fmincon(obj, x0, [], [], [], [], xL, xU, @nonlinConstr, opts);

    % output optimized paramter values
    ansp = round((x_opt ./ scale), 3);
    fprintf('\nFilter parameters for %.f V supply and %.f W reference power are \n', Vdc_ref, P_ref);
    disp('----------------------------------------------------------------------------');
    disp(['grid side inductance:      ', num2str(ansp(1)), ' mH']);
    disp(['filter capacitance:        ', num2str(ansp(2)), ' uF']);
    disp(['converter side inductance: ', num2str(ansp(3)), ' mH']);
    disp(['filter resitance:          ', num2str(ansp(4)), ' ohm']);
    disp(['source capacitance:        ', num2str(ansp(5)), ' uF']);
    disp(['modulation index:          ', num2str(ansp(6))]);
    disp(['theta:                     ', num2str(ansp(7))]);
    disp(' ')

    % graphs for performance indeces w.r.t. reference power
    spect = [];
    P_spect = [];
    P_ref_init = P_ref;
    for i = -(P_ref_init * 0.04):100:(P_ref_init * 0.04)
        P_ref = P_ref_init + i;
        x1 = [P_ref/Vdc_ref];
        [x1_opt] = fmincon(@obj2, x1, [], [], [], [], 0, [], [], opts2);
        x_opt(7) = x1_opt;
        P_spect = [P_spect, P_ref];
        F = objectiveFun(x_opt);
        spect = [spect, F];
    end

    THD_Vac = spect(1,:) .*100;
    THD_Ig = spect(2,:) .*100;
    THD_Ic = spect(3,:) .*100;
    THD_Vdc = spect(4,:) .*100;

    subplot(2,2,1)
    plot(P_spect, THD_Vac)
    xlabel("Reference Power")
    ylabel("Capacitor voltage THD")
    subplot(2,2,4)
    plot(P_spect, THD_Ig)
    xlabel("Reference Power")
    ylabel("Grid side current THD")
    subplot(2,2,3)
    plot(P_spect, THD_Ic)
    xlabel("Reference Power")
    ylabel("Converter side current THD")
    subplot(2,2,2)
    plot(P_spect, THD_Vdc)
    xlabel("Reference Power")
    ylabel("DC source voltage THD")

    function F = objectiveFun(x)

        [Iga, Ica, Vfa, Vdc] = LCL_filter_equation(x);
        
        rip1 = sum(abs(Vfa(h+3:end)).^2);
        rip2 = sum(abs(Iga(h+3:end)).^2);
        rip3 = sum(abs(Ica(h+3:end)).^2);
        rip4 = sum(abs(Vdc(h+2:end)).^2);

        eps_small = 1e-12;
        Vac_rip = real(sqrt(rip1) / max(eps_small, abs(Vfa(h+2))));
        Ig_rip  = real(sqrt(rip2) / max(eps_small, abs(Iga(h+2))));
        Ic_rip  = real(sqrt(rip3) / max(eps_small, abs(Ica(h+2))));
        Vdc_rip = real(sqrt(rip4) / max(eps_small, abs(Vdc(h+1))));

        Idc_f = complex(x(8),0);
        P = (real( Vdc(h+1) * conj(Idc_f) ) - P_ref ) / P_ref;

        F = [Vac_rip; Ig_rip; Ic_rip; Vdc_rip; P];
        
    end

    function F2 = obj2(x)
        y = x_opt;
        y(8) = x;

        [~, ~, ~, Vdc] = LCL_filter_equation(y);
        Idc_f = complex(y(8),0);
        Vdc_f = Vdc(h+1);

        V = (abs(Vdc_f) - (Vdc_ref)) / Vdc_ref;
        P = (real(Vdc_f * conj(Idc_f)) - P_ref ) / P_ref;

        F2 = [V; P];
        weight = [1; 1];
        F2 = sum(weight .* (F2.^2));

    end

    % constraints for fmincon
    function [c, ceq] = nonlinConstr(x)
       
        [~, ~, ~, Vdc] = LCL_filter_equation(x);
        Idc_f = complex(x(8),0);
        Vdc_f = Vdc(h+1);
        eq1 = (abs(Vdc_f) - (Vdc_ref)) / Vdc_ref;
        eq2 = (real(Vdc_f * conj(Idc_f)) - P_ref) / P_ref;
        c = [];
        ceq = [0.1*eq1; 0.1*eq2];

    end

    % steady state equations in EHD
    function [Iga, Ica, Vfa, Vdc] = LCL_filter_equation(x)
    
        % unpacking parameters
        Lg    = x(1);
        C     = x(2);
        Lc    = x(3);
        Rd    = x(4);
        Cdc   = x(5);
        m     = x(6);
        theta = x(7);
        Idc_f = x(8); 

        % fetching EHD matrices
        [D, Ga, Za, Va, Idc] = EHD_matrices(f0, f_sw, m, theta, 2048, h, Vac_peak, Idc_f);
       
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
       
        % check if A is singular, rcond checks reciprocal condition number
        if rcond(A) < 1e-12 
            Iga = NaN * ones(2*h+1, 1);
            Ica = NaN * ones(2*h+1, 1);
            Vfa = NaN * ones(2*h+1, 1);
            Vdc = NaN * ones(2*h+1, 1);
            return
        end
        
        % steady state solution
        Xss = -(A \ (B*U));
        
        % voltage and current harmonics vectors
        Iga = Xss(1:2*h+1);
        Ica = Xss(2*h+1+1:2*(2*h+1));
        Vfa = Xss(2*(2*h+1)+1:3*(2*h+1));
        Vdc = Xss(3*(2*h+1)+1:4*(2*h+1));

    end

end

function [D, Ga, Za, Va, Idc] = EHD_matrices(f0, f_sw, M, theta, N, h, Vac_peak, Idc_f)

    % Inputs:
    %   f0      -   fundamental frequency (Hz)
    %   f_sw    -   switching (carrier) frequency (Hz)
    %   M       -   modulation index (0..1)
    %   theta   -   phase angle in rad
    %   N       -   samples per fundamental period (suggest 2048 or 4096)
    %   h       -   max harmonic index to return (positive integer)
    % Vac_peak  -   peak ac voltage of the grid
    % Idc_f     -   source side DC curent
    
    % time vector of one period
    T = 1/f0;
    Fs = N/T;
    t = (0:N-1)/Fs;
    
    % reference
    w = 2*pi*f0;
    va = Vac_peak*sin(w*t); % phase voltage
    
    % modulating signals
    m_a = M * sin(w*t + theta);
    m_b = M * sin(w*t - (2*pi/3) + theta);
    m_c = M * sin(w*t + (2*pi/3) + theta);
    
    % sawtooth carrier 
    carrier = 2*abs(2*mod(f_sw*t,1)-1)-1;  % triangle in range [-1,1]
    
    % pwm
    s_a = double(m_a > carrier);
    s_b = double(m_b > carrier);
    s_c = double(m_c > carrier);
    
    % FFT 
    Sfull_a = fftshift(fft(s_a)/N);
    Sfull_b = fftshift(fft(s_b)/N);
    Sfull_c = fftshift(fft(s_c)/N);
    Vfull_a = fftshift(fft(va)/N);
    Ifull_l = fftshift(fft(Idc_f)/N);
    
    % extract -h..h (must have N/2 > h)
    center = N/2 + 1;
    idx = center + (-h:h);
    freqs = ( -h:h ) * f0;  % harmonic frequencies 
    
    % make switching matrices for desired h range
    Sa_h = Sfull_a(idx);
    Sb_h = Sfull_b(idx);
    Sc_h = Sfull_c(idx);
    Va = transpose(Vfull_a(idx));
    Idc = complex(zeros(2*h+1,1));
    Idc(h+1) = complex(Idc_f,0);
    
    % generate toeplitz matrices
    Sa = buildToeplitz(Sa_h,h);
    Sb = buildToeplitz(Sb_h,h);
    Sc = buildToeplitz(Sc_h,h);
    
    Ga = Sa - (Sa + Sb + Sc)/3;
    % Gb = Sb - (Sa + Sb + Sc)/3;
    % Gc = Sc - (Sa + Sb + Sc)/3;

    % operational matrix of differentiation
    diagElements = 1j*(-h:h)*w;
    D = diag(diagElements);
    
    
    Za = buildSimplificationMatrix(Sa,h);
    % Zb = buildSimplificationMatrix(Sb,h);
    % Zc = buildSimplificationMatrix(Sc,h);
    
    end
      
    
    function Zx = buildSimplificationMatrix(Sx,h)
    
        m = 6*(-1*ceil(h/3):ceil(h/3));
        check_row = h+1+m;
        Zx = zeros(2*h+1,2*h+1);
        
        % put value when row equal to h+1+m
        for row = (1:2*h+1)
            for col = (1:2*h+1)
                if ~isempty(find(check_row == row, 1))
                    Zx(row,col) = Sx(row, col);
                end
            end
        end
    
    end
    
    
    function Sx = buildToeplitz(Sx_h,h)
        H = -h:h;
        
        % full (2h+1)x(2h+1) by convolution rule
        Nh = 2*h+1;
        Sx = zeros(Nh,Nh);
        for r = (1:Nh)
            for c = (1:Nh)
                shift = H(r)- H(c);    % harmonic difference
                pos = find(H == shift, 1);
                if ~isempty(pos)
                    Sx(r,c) = Sx_h(pos);
                end
            end
        end
    end