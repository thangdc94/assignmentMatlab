close all;
clear all;
clc;
data = randi([0 1],96,1);

hMod = comm.PSKModulator(16, 'BitInput',true);
hDemod = comm.PSKDemodulator(16, 'BitOutput',true);

modSignal = step(hMod, data);
receivedData = step(hDemod, modSignal);

err = sum(abs(data-receivedData))
