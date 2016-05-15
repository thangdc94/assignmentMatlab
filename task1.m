clear all;
close all;
clc;
%% =================================Task 1=================================

%% Initialize parameters
fs = 8000;                      % Tan so lay mau
Ts = 1/fs;                      % Chu ky lay mau
L = 256;                        % so muc luong tu
Nb = log2(L);                   % so bit can cho ma hoa
fm = [200 800 100];             % Tan so tuong ung voi MSV
Am = [1 2 3];                   % Bien do tin hieu
A = 87.6;                       % compression parameter
phi = [0 pi/2 pi/4];            % pha tin hieu
t = 0:Ts:3/min(fm);             % 3 chu ky
Ns = length(t);                 % number of samples
f = (-Ns/2:1:Ns/2-1)/(Ns*Ts);   % Bien doi sang mien tan so  
%% Sampling
s = zeros(1, length(t));                    % intialize
for i = 1:length(fm)
    s = s + Am(i)*cos(2*pi*fm(i)*t+phi(i)); % tong hop tin hieu
end
Amax = max(abs(s));                         % bien do cuc dai cua tin hieu
Sf = fft(s);                                % FFT
Sf = fftshift(Sf);                          % Dich pho tan ve trung tam
%% Compress
x1 = compand(s, A, Amax, 'A/compressor');   % A-law compressor
%% Quantization
delta = 2*Amax/(L - 1);     % buoc luong tu
Mq = -Amax:delta:Amax;      % Gia tri muc luong tu
Ml = 0:L-1;                 % Cac muc luong tu
xcode = zeros(size(x1));
for k = 1:L
    index = find(x1 > Mq(k)-delta/2 & x1 <= Mq(k)+delta/2); % luong tu hoa
    xcode(index) = Ml(k);
end
%% Encoding
code = dec2bin(xcode);          % convert decimal to string binary
code2 = de2bi(xcode);           % convert decimal to binary matrix
data = code2';                  % transform matrix to vector
data = data(:);
%% Decoding
y_index = (bin2dec(code))';     % convert string binary to decimal
%% De-quantization
yq = Mq(y_index + 1);
%% Expand
y = compand(yq, A, Amax, 'A/expander'); % A-law expander
Yf = fft(y);                            % FFT 
Yf = fftshift(Yf);                      % Do something
%% --------------------------------Plotting--------------------------------
figure(1)
subplot(211)
plot(t, s, 'r');
hold on;
stem(t, s);
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
subplot(211)
plot(t, x1, 'r');
hold on;
stem(t, x1);
grid on;
title('Compressed Signal')
xlabel('t');
ylabel('x1(t)');
subplot(212)
stairs(t, xcode);
grid on;
title('Quantized Signal');

figure(3)
stairs(data);
grid on;
title('Encoded signal');
axis([0 100 -2 3]);

figure(4)
subplot(211)
stairs(t, yq);
title('Decoded signal');
grid on;
subplot(212)
stem(t, y);
grid on;
title('Expanded signal')

figure(5)
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