clc, close all, clear all

fs = 5;
flo = 10;
oversampl = 8;
rrcos_filt_len = 8;
rrcos_filt_beta = 0.25;
rrcos_filt = rcosdesign(rrcos_filt_beta, rrcos_filt_len, oversampl);

bits_per_symb = 4;
bits_per_item = 8;
data = randi(2^bits_per_item, 1, 100) - 1;
barker_code13 = [+1 +1 +1 +1 +1 -1 -1 +1 +1 -1 +1 -1 +1];
barker_symbs = (1+i)*barker_code13;
barker_symbs_up = upfirdn(barker_symbs, rrcos_filt, oversampl);

bin_symbs = symbolify(data, bits_per_item, bits_per_symb);

points_per_axis = 2^(bits_per_symb/2);
qam_matrix = (-1:2/(points_per_axis-1):1) + (i*(-1:2/(points_per_axis-1):1)');
qam_lookup = qam_matrix(:)';

symbs = qam_lookup(bin_symbs+1);
packet_header = [barker_symbs zeros(1, oversampl - 1)];
packet = [packet_header symbs];
[data_to_send, t_tx] = transmitter(packet, oversampl, flo, rrcos_filt, fs);

% the channel is HERE
% do I even need HW?
the_num = randi(100) + 10;
some_random_offset = zeros(1, the_num * flo * fs);
data_recieved = awgn([some_random_offset data_to_send], 30);
df = 1+1e-5;
some_phase = 2*pi*(rand()-0.5);

t_rx = 0:1/fs:((length(data_recieved)-1)/fs);

lo_i = sin(df*t_rx+some_phase);
lo_q = cos(df*t_rx+some_phase);

b = fir1(2*flo*fs,2/fs/flo);
rx_with_lo = data_recieved .* lo_i + 1j*data_recieved .* lo_q;

rx_high_sample = upfirdn(b, rx_with_lo, 1, flo*fs);
rx_high_sample_shift = rx_high_sample(2:end-1);

abs_corr = abs(xcorr(rx_high_sample_shift, barker_symbs_up));
[v, corr_max_index] = max(abs_corr);
packet_offset = corr_max_index-length(rx_high_sample_shift);

rx_high_sample_shift_packet = rx_high_sample_shift((packet_offset+1):end);

rx = upfirdn(rrcos_filt, rx_high_sample_shift_packet, 1, oversampl);

rx_shift = rx((rrcos_filt_len + 1):(end - rrcos_filt_len));
rx_without_header = rx_shift(length(packet_header)+1:end);
rx_header = rx_shift(1:length(barker_symbs));
phase_error = mean(angle(rx_header.*barker_symbs));

rx_phase_corrected = rx_without_header*exp(i*(pi/2-phase_error));
gain_comp = mean(abs(barker_symbs./rx_header));

symbs_rx = rx_phase_corrected * gain_comp;

figure
hold on
plot(symbs_rx, 'o');
plot(symbs, 'o')
grid on
axis([-1.5 1.5 -1.5 1.5])

rx_bin_symbs = zeros(size(symbs));
for k = 1:length(symbs_rx)
   [dont_care, qam_i] = min(abs(symbs_rx(k)-qam_lookup));
   resolved_symb = qam_lookup(qam_i);
   
   rx_bin_symbs(k) = qam_i - 1;
end

rx_data = symbolify(rx_bin_symbs, bits_per_symb, bits_per_item);
rx_data = rx_data(1:length(data));
[bit_errors, ber] = biterr(rx_data,data);

