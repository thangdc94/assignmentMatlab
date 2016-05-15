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
t = 0:Ts:10;                    % 10s
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
code = de2bi(xcode);  % convert decimal to binary matrix
%% Convert binary matrix to data
data = code';
data = data(:);

%% =================================Task 2=================================

%% Initialize parameters
Rb = fs*log2(L);            % bit rate
T_sym = 1/Rb;               % chu ky ky hieu
Es = 2/3*T_sym;             % nang luong ky hieu
M = 8;                      % so muc dieu che 8-PSK
SNR = [5 8 12];             % SNR[dB]
%% Modulation
% Create a 8-PSK modulator System object with bits as inputs 
% and Gray-coded signal constellation
hModulator = comm.PSKModulator(M, 'BitInput', true);
% Modulate
d = step(hModulator, data);         % su dung modulator

% Baseband-equivalent Model
tt = 0:T_sym/50:T_sym;
p = sqrt(2*Es/T_sym)*(1-cos(2*pi*tt/T_sym));

% Xac dinh dang tin hieu phat voi 1000 ky hieu dau tien
numofSym = 1000;                    % so symbol dung cho mo phong
tt = 0:T_sym/50:T_sym*numofSym;     % time domain
pp = repmat(p, numofSym, 1);        % clone vector p
s_mod = diag(d(1:numofSym))*pp;     % d*p
% convert matrix to vector
s_mod = s_mod';
s_mod = s_mod(:)';

% Pho tin hieu
Ns2 = length(s_mod);                    
f2 = (-Ns2/2:1:Ns2/2-1)/(Ns2*T_sym/50); % frequency domain
Sf_mod = fft(s_mod);                    % FFT
Sf_mod = fftshift(Sf_mod);

%% Transmit signal through an AWGN channel.
d_noise = zeros(length(d), length(SNR));
for k = 1:length(SNR)                   % vong lap cho cac SNR khac nhau
    d_noise(:, k) = awgn(d, SNR(k));    % nhieu tac dong vao cac ky hieu
end
y_noise = awgn(s_mod, SNR(2));          % nhieu tac dong vao dang tin hieu
Yf_noise = fft(y_noise);                % FFT
Yf_noise = fftshift(Yf_noise);

%% Demodulation
d_demod = zeros(length(data), length(SNR));
% Create a 8-PSK demodulator System object with bits as outputs 
% and Gray-coded signal constellation
hDemod = comm.PSKDemodulator(M, 'BitOutput', true);
for k = 1:length(SNR)
    d_demod(:, k) = step(hDemod, d_noise(:, k));    % su dung demodulator
end

%% Decoding
y_index = zeros(length(SNR), length(xcode));
for k = 1: length(SNR)                          % a loop
    tmp = vec2mat(d_demod(:, k)', Nb);
    y_index(k, :) = bi2de(tmp)';                % binary to decimal
end

%% De-quantization
yq = Mq(y_index + 1);
%% Expand
y = zeros(length(SNR), length(s));
Yf = y;
for k = 1:length(SNR)                                   % loop again
    y(k, :) = compand(yq(k, :), A, Amax, 'A/expander'); % A-law expander
    Yf(k, :) = fft(y(k, :));                            % FFT
    Yf(k, :) = fftshift(Yf(k, :));
end

%% --------------------------------Plotting--------------------------------

%% Task 2b: BER
BER = zeros(length(SNR), 1);                                % initialize
for k = 1:length(SNR)
    BER(k) = sum(abs(d_demod(:, k) - data))/length(data);   % tinh BER
end
figure(1)
plot(SNR, BER);
hold on;
plot(SNR, BER, 'ro');
grid on;
title(['Uoc tinh xac suat loi voi SNR = ', num2str(SNR)]);
xlabel('SNR[dB]');
ylabel('BER');

%% Task 2c

% Constellation
scatterplot(d,1,0,'or');
grid on;

h = scatterplot(d_noise(:, 2),1,0,'xb');
hold on;
scatterplot(d,1,0,'or', h);
grid on;

% Dang tin hieu va pho
figure(4)
subplot(2,1,1)
stem([real(s_mod) imag(s_mod)])
title('Dang xung phat')
xlim([0 350]);
subplot(2,1,2);
plot(f2, Sf_mod);
title('Pho tin hieu sau dieu che')
grid on;

figure(5)
subplot(2,1,1)
stem([real(y_noise) imag(y_noise)])
hold on
plot([real(s_mod) imag(s_mod)],'color','r')
xlim([0 350]);
title('Dang tin hieu tai bo thu')
subplot(2,1,2)
plot(f2, Yf_noise);
title('Pho tin hieu tai bo thu')
grid on;

% Eye diagram
eyediagram([real(s_mod) imag(s_mod)],length(p)*2)
title('Bieu do mat tin hieu phat')
grid on;
eyediagram([real(y_noise) imag(y_noise)],length(p)*2)
title('Bieu do mat tin hieu tai bo thu')
grid on;

%% Task 2d
figure(8)
subplot(length(SNR)+1, 1, 1)
plot(t, s);
xlim([0 0.08]);
title('Original signal');
xlabel('t');
ylabel('s(t)');
grid on;
for k = 1:length(SNR)
    subplot(length(SNR)+1, 1, k+1)
    plot(t, y(k, :));
    xlim([0 0.08]);
    title(['Reconstruct signal with SNR = ', num2str(SNR(k)), 'dB']);
    xlabel('t');
    ylabel('y(t)');
    grid on;
end
