function [t, s, S_h, freqs] = switchh(f0, f_sw, M, N, h)
% generate binary switching signals for 3-phase SPWM

% Inputs:
%   f0      -   fundamental frequency (Hz)
%   f_sw    -   switching (carrier) frequency (Hz)
%   M       -   modulation index (0..1)
%   N       -   samples per fundamental period (suggest 2048 or 4096)
%   h       -   max harmonic index to return (positive integer)

% Outputs:
%   t       - time vector (1 x N)
%   s       - 3 x N matrix of switching signals (rows: a,b,c) values 0/1
%   S_h     - 3 x (2*h+1) complex Fourier coefficients for harmonics [-h..h]
%   freqs   - harmonic frequencies vector (Hz) for S_h columns

% time vector of one period
T = 1/f0;
Fs = N/T;
t = (0:N-1)/Fs;

% reference
w = 2*pi*f0;
m_a = M*sin(w*t);
m_b = M*sin(w*t - (2*pi/3));
m_c = M*sin(w*t + (2*pi/3));

% sawtooth carrier 
carrier = 2*abs(2*mod(f_sw*t,1)-1)-1;  % triangle in range [-1,1]

% pwm
s_a = double(m_a > carrier);
s_b = double(m_b > carrier);
s_c = double(m_c > carrier);
s = [s_a; s_b; s_c];

% FFT 
Sfull_a = fftshift(fft(s_a)/N);
Sfull_b = fftshift(fft(s_b)/N);
Sfull_c = fftshift(fft(s_c)/N);

% map bins to harmonic index -N/2 .. N/2-1
k = (-N/2):(N/2-1);
% extract -h..h (must have N/2 > h)
center = N/2 + 1;
idx = center + (-h:h);

S_h = [Sfull_a(idx); Sfull_b(idx); Sfull_c(idx)];
freqs = ( -h:h ) * f0;  % harmonic frequencies 

% compute gi
Sa = S_h(1,:);
Sb = S_h(2,:);
Sc = S_h(3,:);
Ga = Sa - (Sa + Sb + Sc)/3;
Gb = Sb - (Sa + Sb + Sc)/3;
Gc = Sc - (Sa + Sb + Sc)/3;

col = Ga(h+1:end);
row = Ga(h+1:-1:1);
G_top = toeplitz(col,row);

disp(G_top)

end

f = switchh(50, 5000, 0.8, 2048, 2);