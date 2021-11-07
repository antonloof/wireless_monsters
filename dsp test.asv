clc, close all, clear all

flo = 20;
oversampl = 6;
rrcos_filt_len = 6;
rrcos_filt_beta = 0.25;
rrcos_filt = rcosdesign(rrcos_filt_beta, rrcos_filt_len, oversampl);

bits_per_symb = 4;
bits_per_item = 8;
data = randi(2^bits_per_item - 1, 1, 10);
data = [10 20 30 40 50 60 70 80 90];
barker_code13 = [+1 +1 +1 +1 +1 -1 -1 +1 +1 -1 +1 -1 +1];
barker_symbs = barker_code13 + i*barker_code13;

bin_symbs = symbolify(data, bits_per_item, bits_per_symb);

points_per_axis = 2^(bits_per_symb/2);
qam16_matrix = (-1:2/(points_per_axis-1):1) + (i*(-1:2/(points_per_axis-1):1)');
qam16_lookup = reshape(qam16_matrix.',1,[]);

symbs = qam16_lookup(bin_symbs+1);
symbs_with_barker = [zeros(1, 100) barker_symbs symbs];
[data_to_send, t_tx] = transmitter(symbs_with_barker, oversampl, flo, rrcos_filt);
some_random_offset = zeros(1, (randi(100) + 10) * flo);
some_random_offset = zeros(1, 100*oversampl*flo);
data_recieved = awgn([some_random_offset data_to_send], 20); % do I even need HW?
data_recieved = [some_random_offset data_to_send];

df = 1;
some_phase = 0;
t_rx = linspace(0, length(data_recieved)/flo/oversampl,length(data_recieved));
lo_i = sin(2*pi*flo*df*t_rx+some_phase);
lo_q = cos(2*pi*flo*df*t_rx+some_phase);

b = fir1(2*flo,1/2, 'low');
rx_with_lo = data_recieved .* lo_i + 1j*data_recieved .* lo_q;

rx_high_sample = upfirdn(b, rx_with_lo, 1, flo);

baker_symbs_upsampled = upfirdn(barker_symbs, rrcos_filt, oversampl, 1);
all_symbs_upsampled = upfirdn(symbs_with_barker, rrcos_filt, oversampl, 1);

baker_filter = flip(barker_symbs);

hold on
plot(abs(filtfilt(baker_filter, 1, symbs_with_barker)))
plot(abs(xcorr(symbs_with_barker, barker_symbs)))

rx = upfirdn(rrcos_filt, rx_high_sample, 1, oversampl);
rx_shift = rx((rrcos_filt_len + 1):(end - rrcos_filt_len));

scale_factor = 0.9;
symbs_rx = rx_shift/max(abs(rx_shift))/scale_factor*sqrt(2);

figure
hold on
plot(symbs_rx, 'o');
plot(symbs, 'o')
grid on
axis([-1.5 1.5 -1.5 1.5])

rx_bin_symbs = zeros(size(symbs));
for k = 1:length(symbs_rx)
   [dont_care, qam_i] = min(abs(symbs_rx(k)-qam16_lookup));
   rx_bin_symbs(k) = qam_i - 1;
end

rx_data = symbolify(rx_bin_symbs, bits_per_symb, bits_per_item);
rx_data = rx_data(1:length(data));
[bit_errors, ber] = biterr(rx_data,data);

