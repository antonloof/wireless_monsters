clc, close all, clear all

data = [0.02 -12.5
0.05 -4.4
0.07 -1.5
0.1  1.4
0.15 4.8
0.2  8.0
0.25 9.8
0.3  10.6
0.35 11.1
0.4  11.4
0.45 11.65
0.5  11.9];


in_ampl_dbm = 20*log10(data(:,1))+2;
out_ampl_dbm = data(:,2);
plot(in_ampl_dbm, out_ampl_dbm, 'o')
xlabel("input power [dBm]")
ylabel("output power [dBm]")
title("Transmitter transfer function")
p = polyfit(in_ampl_dbm(1:7), out_ampl_dbm(1:7),1);
hold on
plot(in_ampl_dbm, polyval(p, in_ampl_dbm))
legend("measured transfter", "linear fit")
