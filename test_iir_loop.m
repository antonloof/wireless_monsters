clc, close all, clear all
sig = sin((1:1000)/10);
a = awgn(sig, 10);
a_filt = zeros(size(a));
a_filt(1) = a(1);


new_factor = 0.1;
for i = 2:length(a)
    a_filt(i) = new_factor * a(i) + (1 - new_factor) * a_filt(i - 1);
end

hold on
grid on
plot(a)
plot(a_filt)

figure
hold on

fs = 1e3;
fft_a = fft(a);
f = (0:1:length(fft_a)-1)*fs/length(fft_a);
plot(f, abs(fft_a))
plot(f, abs(fft(a_filt)))