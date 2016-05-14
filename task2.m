clear all;
close all;
clc;
%% ==============================================Task 1=========================================================

%% Initialize parameters
fs = 8000;
Ts = 1/fs;
L = 256; % so muc luong tu
Nb = log2(L); % so bit can cho ma hoa
fm = [200 800 100];
Am = [1 2 3];
A = 87.6; % compression parameter
phi = [0 pi/2 pi/4];
t = 0:Ts:0.1; % 10s
Ns = length(t);
f = (-Ns/2:1:Ns/2-1)/(Ns*Ts);
%% Sampling
s = zeros(1, length(t));
for kk = 1:length(fm)
    s = s + Am(kk)*cos(2*pi*fm(kk)*t+phi(kk));
end
Amax = max(abs(s));
Sf = fft(s);
Sf = fftshift(Sf);
%% Compress
x1 = compand(s, A, Amax, 'A/compressor');
%% Quantization
delta = 2*Amax/(L - 1);   % buoc luong tu
Mq = -Amax:delta:Amax;  % Gia tri muc luong tu
Ml = 0:L-1; % Cac muc luong tu
xcode = zeros(size(x1));
for k = 1:L
    index = find(x1 > Mq(k)-delta/2 & x1 <= Mq(k)+delta/2);
    xcode(index) = Ml(k);
end
%% Encoding
code = de2bi(xcode);  % chuyen tu decimal sang binary

%% ==============================================Task 2=========================================================

%% Convert binary matrix to data
data = code';
data = data(:);
%% Initialize parameters
Rb = fs*log2(L);
T_sym = 1/Rb;
Es = 2/3*T_sym;
M = 8; % so muc dieu che 8-PSK
SNR = [5 8 12]; % [dB]
%% Modulation
% Create a 8-PSK modulator System object with bits as inputs and Gray-coded signal constellation
hModulator = comm.PSKModulator(M, 'BitInput', true);
% Modulate and plot the data
d = step(hModulator, data);

% Mo hinh tuong duong bang goc
tt = 0:T_sym/50:T_sym;
p = sqrt(2*Es/T_sym)*(1-cos(2*pi*tt/T_sym));
plot(tt, p);

% Xac dinh dang tin hieu phat voi 10 ky hieu dau tien
pp = repmat(p,10,1);    % clone vector p
s_mod = diag(d(1:10))*pp;
% convert matrix to vector
s_mod = s_mod';
s_mod = s_mod(:)';

%% Transmit signal through an AWGN channel.
d_noise = zeros(size(d), size(SNR));
for k = 1:length(SNR)
    d_noise(:, k) = awgn(d, SNR(k), 'measured');
end
y_noise = awgn(s_mod, SNR(2), 'measured');

%% Demodulation
d_demod = zeros(length(data), length(SNR));
hDemod = comm.PSKDemodulator(M, 'BitOutput', true);
for k = 1:length(SNR)
    d_demod(:, k) = step(hDemod, d_noise(:, k));
end

%% Tinh BER
BER = zeros(length(SNR), 1);
for k = 1:length(SNR)
    BER(k) = sum(abs(d_demod(:, k) - data))/length(data);
end

%% Decoding
y_index = zeros(length(SNR), length(xcode));
for k = 1: length(SNR)
    tmp = vec2mat(d_demod(:, k)', Nb);
    y_index(k, :) = bi2de(tmp)';
end

%% De-quantization
yq = Mq(y_index + 1);
%% Expand
y = zeros(length(SNR), length(s));
Yf = y;
for k = 1:length(SNR)
    y(k, :) = compand(yq(k, :), A, Amax, 'A/expander');
    Yf(k, :) = fft(y(k, :));
    Yf(k, :) = fftshift(Yf(k, :));
end

%% Plotting
h = scatterplot(d_noise(2, :),1,0,'xb');
hold on;
scatterplot(d,1,0,'or',h);
grid on;

figure(3)
plot(SNR, BER);

figure(4)
subplot(2,1,1)
stem([real(s_mod) imag(s_mod)])
title('Dang xung phat')
subplot(2,1,2)
stem([real(y_noise) imag(y_noise)])
hold on
plot([real(s_mod) imag(s_mod)],'color','r')
title('Dang tin hieu tai bo thu')

% Eye diagram
eyediagram([real(s_mod) imag(s_mod)],length(p)*2)
title('Bieu do mat tin hieu phat')
eyediagram([real(y_noise) imag(y_noise)],length(p)*2)
title('Bieu do mat tin hieu tai bo thu')

% Task 2d
figure(7)
subplot(211)
plot(t, s);
title('Original signal');
xlabel('t');
ylabel('s(t)');
grid on;
subplot(212);
plot(f, Sf);
xlabel('f');
ylabel('S(f)')
grid on;

figure(8)
subplot(211);
plot(t, y(1,:));
title('Reconstruct signal');
xlabel('t');
ylabel('y(t)')
grid on;
subplot(212);
plot (f, Yf(1, :));
xlabel('f');
ylabel('Y(f)')
grid on;

figure(9)
subplot(length(SNR)+1, 1, 1)
plot(t, s);
title('Original signal');
xlabel('t');
ylabel('s(t)');
grid on;
for k = 1:length(SNR)
    subplot(length(SNR)+1, 1, k+1)
    plot(t, y(k, :));
    title(['Reconstruct signal SNR = ', num2str(SNR(k)), 'dB']);
    xlabel('t');
    ylabel('y(t)');
    grid on;
end
