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
    % THD_Vac_ref  = 0.06; 
    % THD_Ig_ref = 0.02;
    % THD_Ic_ref = 0.1;
    % THD_Vdc_ref = 0.1;
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

    opts2 = optimoptions('fmincon', ...
    'Algorithm','interior-point', ...
    'Display','none', ...
    'MaxFunctionEvaluations',1e5, ...
    'MaxIterations',5000, ...
    'OptimalityTolerance',1e-6, ...
    'StepTolerance',1e-6, ...
    'ConstraintTolerance',1e-6, ...
    'FiniteDifferenceType','forward', ...
    'FiniteDifferenceStepSize',1e-5, ...
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
    
    % [~, ceq_test] = nonlinConstr(x0);
    % fprintf('Initial constraint errors: Vdc=%.3f, Power=%.3f\n', ceq_test(1), ceq_test(2));
    
    % calling fmincon
    [x_opt, fval, exitflag, output] = fmincon(obj, x0, [], [], [], [], xL, xU, @nonlinConstr, opts);

    % % diagnostic prints
    % disp('--- DIAGNOSTICS ---');
    % disp(['exitflag: ', num2str(exitflag)]);
    % if exist('output','var')
    %     if isfield(output,'constrviolation'); disp(['constrviolation: ', num2str(output.constrviolation)]); end
    %     if isfield(output,'message'); disp(['message: ', output.message]); end
    % end
    % 
    % % evaluate constraints at x_opt
    % [c_opt, ceq_opt] = nonlinConstr(x_opt);
    % disp('nonlinear inequality (c) at x_opt:'); disp(c_opt(:).');
    % disp('nonlinear equality (ceq) at x_opt:'); disp(ceq_opt(:).');
    % 
    % % evaluate steady state
    % [Iga_opt, Ica_opt, Vfa_opt, Vdc_opt] = LCL_filter_equation(x_opt);
    % fprintf('|Vdc(fund)| = %.6f V\n', abs(Vdc_opt(h+1)));
    % fprintf('x_opt(7) (Idc variable) = %.6f A (scaled)\n', x_opt(7));
    % % compute actual P using phasor conj for safety
    % Pfund = real( Vdc_opt(h+1) * conj( complex(x_opt(7),0) ) );
    % fprintf('P_fund (using x_opt(7)) = %.6f W\n', Pfund);
    % 
    % [~, ceq_test] = nonlinConstr(x_opt);
    % fprintf('Final constraint errors: Vdc=%.3f, Power=%.3f\n', ceq_test(1), ceq_test(2));

    % [Iga, Ica, Vfa, Vdc] = LCL_filter_equation(x_opt);
    % disp(real(abs(Vdc(h+1))));
    % Ix = x_opt(7);
    % Px = real(Vdc(h+1) * Ix);
    % disp(Px);
    % s = objectiveFun(x_opt);
    % disp(s);


    [Iga, Ica, Vfa, Vdc] = LCL_filter_equation(x_opt);
    
%% ------------------ Robust two-period reconstruction & plots ------------------
% Assumes X_h vectors are ordered from -h .. 0 .. +h (length 2*h+1; center index h+1).
% Replace prior reconstruction calls with this block.

% get harmonic vectors (already computed)
% [Iga, Ica, Vfa, Vdc] = LCL_filter_equation(x_opt);

% Sampling and time vector (two periods)
numPeriods = 400;
Tfund = 1 / f0;
% pick Fs so that Fs >> switching freq (to show waveform detail). 10 kHz is fine.
Fs = 10000;                    
t = 0 : 1/Fs : numPeriods*Tfund - 1/Fs;   % exactly two periods
t = t*100;

% Helper: vectorized inverse EHD (harmonic sum)
% NOTE: this uses direct sum: x(t) = real( sum_{k=-h}^{h} X_h[k]*exp(j*2*pi*k*f0*t) )
invEHD = @(X_h) real( (reshape(X_h,1,[]) * exp(1j*2*pi*((-h:h)')*t)) );

% Reconstruct
ig_t  = invEHD(Iga);    % grid-side current time waveform
ic_t  = invEHD(Ica);    % converter-side current waveform
vf_t  = invEHD(Vfa);    % filter capacitor voltage
vdc_t = invEHD(Vdc);    % DC-link voltage (should be mostly DC + tiny ripple)

% Basic plots: two periods
figure('Units','normalized','Position',[0.05 0.05 0.85 0.8],'Color','w');

subplot(4,1,1);
plot(t, ig_t, 'LineWidth', 1.0);
xlabel('Time (s)'); ylabel('i_g (A)');
title('Grid-side current — 2 periods'); xlim([0 numPeriods*Tfund]);

subplot(4,1,2);
plot(t, ic_t, 'LineWidth', 1.0); hold on;
plot(t(1:round(length(t)/2)), ic_t(1:round(length(t)/2)), 'k:');
plot(t(round(length(t)/2)+1:end), ic_t(round(length(t)/2)+1:end), 'r--');
hold off;
xlabel('Time (s)'); ylabel('i_c (A)');
title('Converter-side current — 2 periods'); xlim([0 numPeriods*Tfund]);

subplot(4,1,3);
plot(t, vf_t, 'LineWidth', 1.0); hold on;
plot(t(1:round(length(t)/2)), vf_t(1:round(length(t)/2)), 'k:');
plot(t(round(length(t)/2)+1:end), vf_t(round(length(t)/2)+1:end), 'r--');
hold off;
xlabel('Time (s)'); ylabel('v_f (V)');
title('Filter capacitor voltage — 2 periods'); xlim([0 numPeriods*Tfund]);

subplot(4,1,4);
plot(t, vdc_t, 'LineWidth', 1.0); hold on;
plot(t(1:round(length(t)/2)), vdc_t(1:round(length(t)/2)), 'k:');
plot(t(round(length(t)/2)+1:end), vdc_t(round(length(t)/2)+1:end), 'r--');
hold off;
xlabel('Time (s)'); ylabel('V_{dc} (V)');
title('DC-link voltage — 2 periods'); xlim([0 numPeriods*Tfund]);

% Optional: zoom into a small window to inspect switching ripple (e.g., first 2 ms)
zoomWindow = 0.002; % 2 ms
figure('Units','normalized','Position',[0.1 0.1 0.6 0.35]);
plot(t, ig_t, 'LineWidth', 1.0); xlim([0 zoomWindow]);
xlabel('Time (s)'); ylabel('i_g (A)');
title(sprintf('Grid current zoom (first %.1f ms)', zoomWindow*1e3));

% Optional: compute and print basic stats to sanity-check
% fprintf('Signal stats (grid current): mean=%.4f A, pk2pk=%.4f A\n', mean(ig_t), max(ig_t)-min(ig_t));
% fprintf('DC-link: mean=%.4f V, ripple-pk2pk=%.6f V\n', mean(vdc_t), max(vdc_t)-min(vdc_t));


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


    % function X = inverseFFT(X_h, f, h)
    % 
    %     Nfft = floor(1+1/(f*0.0001));
    % 
    %     FFTbins = zeros(2*Nfft,1);
    %     FFTbins(1) = X_h(h+1); %dc term
    %     for k = 1:h
    %         % for positive k
    %         FFTbins(k+1) = X_h(h+1+k);          
    %         % conjugate for negative k
    %         FFTbins(Nfft - k + 1) = conj(X_h(h+1-k)); 
    %     end
    %     X = real(ifft(FFTbins) * Nfft);
    % 
    % end

    % function x_t = inverseEHD_autodetect(X_h, f0, h, t)
    % % Reconstruct time waveform from EHD harmonic vector X_h (length 2*h+1).
    % % This routine tests two common conventions and picks the one that yields
    % % a dominant fundamental magnitude.
    % %
    % % Inputs:
    % %   X_h : (2*h+1)x1 harmonic vector (complex). Commonly ordered -h..+h with center at h+1.
    % %   f0  : fundamental frequency (Hz)
    % %   h   : maximum harmonic index
    % %   t   : time vector (1xNt)
    % % Output:
    % %   x_t : real time waveform (1xNt)
    % %
    % % It will also print which mapping was chosen.
    % 
    %     X_h = X_h(:);            % ensure column
    %     % Try mapping A: assume X_h indices correspond to k = -h:h (center index = h+1)
    %     ksA = (-h:h).';          % harmonic indices for mapping A
    %     % Try mapping B: sometimes data is shifted by +1 (rare). We'll try also k = -h+1 : h+1
    %     ksB = (-h+1:h+1).';      % shifted indices (center at h+2)
    % 
    %     % build time exponentials and compute candidate waveforms
    %     E_A = exp(1j*2*pi*(ksA)*t);  % (2h+1) x Nt
    %     try
    %         xA_complex = X_h.' * E_A;
    %         xA = real(xA_complex);
    %     catch
    %         xA = -inf(1,length(t)); % fail-safe
    %     end
    % 
    %     % For mapping B we need to align sizes: if ksB has same length, shift X_h accordingly
    %     if length(ksB) == length(X_h)
    %         E_B = exp(1j*2*pi*(ksB)*t);
    %         xB_complex = X_h.' * E_B;
    %         xB = real(xB_complex);
    %     else
    %         % if lengths don't match (shouldn't happen) make xB empty
    %         xB = [];
    %     end
    % 
    %     % Compute FFT magnitude near f0 for both candidates to find which has dominant fundamental
    %     % FFT of one period (use 1/f0 seconds)
    %     NtFund = round(1/(f0) * (1/(t(2)-t(1))));
    %     if NtFund < 8, NtFund = 1024; end
    %     % window a single period from the start
    %     idx1 = 1: min(NtFund, length(t));
    %     magFundA = dominantFundamentalMag(xA(idx1), f0, t(idx1));
    %     magFundB = dominantFundamentalMag(xB(idx1), f0, t(idx1));
    % 
    %     fprintf('Detected fundamental magnitudes -> A: %.6e , B: %.6e\n', magFundA, magFundB);
    % 
    %     if magFundA >= 10*magFundB  % A much better
    %         chosen = 'A (center = h+1)';
    %         x_t = xA;
    %     elseif magFundB > 10*magFundA
    %         chosen = 'B (shifted center)';
    %         x_t = xB;
    %     else
    %         % ambiguous: pick the larger magnitude
    %         if magFundA >= magFundB
    %             chosen = 'A (fallback)';
    %             x_t = xA;
    %         else
    %             chosen = 'B (fallback)';
    %             x_t = xB;
    %         end
    %     end
    %     fprintf('inverseEHD_autodetect selected mapping: %s\n', chosen);
    % end
    
    % function mag = dominantFundamentalMag(x_segment, f0, t_segment)
    %     % return magnitude of fundamental component using FFT peak around f0
    %     N = length(x_segment);
    %     Xfft = fft(x_segment);
    %     f = (0:N-1)/ (t_segment(end)-t_segment(1) + (t_segment(2)-t_segment(1))) ; % simple freq vector
    %     % find index nearest to f0
    %     [~, idx] = min(abs(f - f0));
    %     mag = abs(Xfft(idx));
    % end


end
