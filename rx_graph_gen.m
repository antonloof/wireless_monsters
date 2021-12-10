clc, close all, clear all

data = [-130,-74
-120,-63.6
-110,-53.4
-100,-44.9
-90,-34.9
-80,-24.9
-70,-14.9
-60,-4.2
-50,5.7
-48,7.7
-46,9.7
-44,11.9
-42,14
-40,15.7
-39,16.3
-38,16.8
-37,17.25
-36,17.55
-35,17.6];

input_dbm = data(:,1);
output_dbm = data(:,2);
hold on
plot(input_dbm, output_dbm, 'o');
xlabel("input power [dBm]")
ylabel("output power [dBm]")
title("Receiver transfer function")
p = polyfit(input_dbm(1:end-4), output_dbm(1:end-4), 1);
plot(input_dbm, polyval(p, input_dbm))
legend("measured transfter", "linear fit")
