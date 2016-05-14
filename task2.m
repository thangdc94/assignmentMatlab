clear all;
close all;
clc;
%% ==============================================Task 1=========================================================

%% Initialize parameters
fs = 8000;
Ts = 1/fs;
L = 16; % so muc luong tu
Nb = log2(L); % so bit can cho ma hoa
fm = [200 800 100];
Am = [1 2 3];
A = 87.6; % compression parameter
phi = [0 pi/2 pi/4];
t = 0:Ts:10; % 10s
Ns = length(t);
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
data = code';
data = data(:);

%% Convert binary matrix to data
% data = cell2mat(cellstr(code)');

%% ==============================================Task 2=========================================================

%% Initialize parameters
Rb = fs*log2(L);
T_sym = 1/Rb;
Es = 2/3*T_sym;
M = 8; % so muc dieu che 8-PSK
SNR = 15; % [dB]

% Create a 8-PSK modulator System object with bits as inputs and Gray-coded signal constellation
hModulator = comm.PSKModulator(M,'BitInput',true);
% Modulate and plot the data
dk = step(hModulator, data);
    

% symbols = {'000', '001', '011', '010', '110', '111', '101', '100'};
% value = [exp(0), exp(j*pi/4), exp(j*pi/2), exp(j*3*pi/4), exp(j*pi), exp(j*5*pi/4), exp(j*3*pi/2), exp(j*7*pi/4)];
% %% Modulation
% data_mod = strread(data,'%3s')';    % split data into 3 bit bundles
% dk = zeros(size(data_mod));
% for kk = 1:M
%     ind = find(strcmp(data_mod,symbols(kk)));
%     dk(ind) = value(kk);
% end
% 
%% Mo hinh tuong duong bang goc
t = 0:T_sym/50:T_sym;
p = sqrt(2*Es/T_sym)*(1-cos(2*pi*t/T_sym));
plot(t, p);
%% Transmit signal through an AWGN channel.
ynoisy = awgn(dk, SNR, 'measured');
%% Plotting
h = scatterplot(ynoisy,1,0,'xb');
hold on;
scatterplot(dk,1,0,'or',h);
grid on;

