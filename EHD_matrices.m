function [D, Ga, Za, Va, Il] = EHD_matrices(f0, f_sw, M, N, h)

% Inputs:
%   f0      -   fundamental frequency (Hz)
%   f_sw    -   switching (carrier) frequency (Hz)
%   M       -   modulation index (0..1)
%   N       -   samples per fundamental period (suggest 2048 or 4096)
%   h       -   max harmonic index to return (positive integer)

% time vector of one period
T = 1/f0;
Fs = N/T;
t = (0:N-1)/Fs;

% reference
w = 2*pi*f0;
% va = 120*sin(w*t); % phase voltage given for L_filter
va = 310*sin(w*t); % phase voltage given for LCL_filter
il = 50*cos(w*t);

m_a = M*sin(w*t);
m_b = M*sin(w*t - (2*pi/3));
m_c = M*sin(w*t + (2*pi/3));

% sawtooth carrier 
carrier = 2*abs(2*mod(f_sw*t,1)-1)-1;  % triangle in range [-1,1]

% pwm
s_a = double(m_a > carrier);
s_b = double(m_b > carrier);
s_c = double(m_c > carrier);
% s = [s_a; s_b; s_c];

% FFT 
Sfull_a = fftshift(fft(s_a)/N);
Sfull_b = fftshift(fft(s_b)/N);
Sfull_c = fftshift(fft(s_c)/N);
Vfull_a = fftshift(fft(va)/N);
Ifull_l = fftshift(fft(il)/N);

% map bins to harmonic index -N/2 .. N/2-1
% k = (-N/2):(N/2-1);

% extract -h..h (must have N/2 > h)
center = N/2 + 1;
idx = center + (-h:h);
freqs = ( -h:h ) * f0;  % harmonic frequencies 

% make switching matrices for desired h range
Sa_h = Sfull_a(idx);
Sb_h = Sfull_b(idx);
Sc_h = Sfull_c(idx);
Va = transpose(Vfull_a(idx));
Il = transpose(Ifull_l(idx));

% S_h = [Sfull_a(idx); Sfull_b(idx); Sfull_c(idx)];
% Sa = S_h(1,:);
% Sb = S_h(2,:);
% Sc = S_h(3,:);

% generate toeplitz matrices
Sa = buildToeplitz(Sa_h,h);
Sb = buildToeplitz(Sb_h,h);
Sc = buildToeplitz(Sc_h,h);

Ga = Sa - (Sa + Sb + Sc)/3;
% Gb = Sb - (Sa + Sb + Sc)/3;
% Gc = Sc - (Sa + Sb + Sc)/3;

% identity matrix
% I = eye(2*h+1);

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
% build Toeplitz 
% c = Sx_h(h+1:end);
% r = Sx_h(h+1:-1:1);
% Sx = toeplitz(c, r);

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

