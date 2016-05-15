close all;
clear all;
clc;
data = randi([0 1],10e6,1);
SNR = 8:12; %khoang gia tri cua SNR 

hMod = comm.PSKModulator(16, 'BitInput',true);
hDemod = comm.PSKDemodulator(16, 'BitOutput',true);
modSignal = step(hMod, data);
ber = zeros(1, length(SNR));
for i= 1:length(SNR)
    s_awgn = awgn(modSignal, SNR(i));
    receivedData = step(hDemod, s_awgn);
    ber(i) = (sum(abs(data-receivedData)))/length(data);
end

plot(SNR, ber)
