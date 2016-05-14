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
t = 0:Ts:3/min(fm); % 3 chu ky
Ns = length(t);
f = (-Ns/2:1:Ns/2-1)/(Ns*Ts);
%% Sampling
s = zeros(1, length(t));
for i = 1:length(fm)
    s = s + Am(i)*cos(2*pi*fm(i)*t+phi(i));
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
code = dec2bin(xcode);  % chuyen tu decimal sang string binary
%% Decoding
y_index = (bin2dec(code))'; % chuyen tu string binary sang decimal
%% De-quantization
yq = Mq(y_index + 1);
%% Expand
y = compand(yq, A, Amax, 'A/expander');
Yf = fft(y);
Yf = fftshift(Yf);
%% --------------------------------Plotting--------------------------------
figure(1)
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

figure(2)
subplot(211);
plot(t, y);
title('Reconstruct signal');
xlabel('t');
ylabel('y(t)')
grid on;
subplot(212);
plot (f, Yf);
xlabel('f');
ylabel('Y(f)')
grid on;