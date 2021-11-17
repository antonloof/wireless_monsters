clc, close all, clear all

fs = 5;
lo_up = 10;
flo = 1;
oversampl = 2;

rrcos_filt_len = 8;
rrcos_filt_beta = 0.5;
rrcos_filt = rcosdesign(rrcos_filt_beta, rrcos_filt_len, oversampl);

bits_per_symb = 6;
bits_per_item = 8;
data_count = 20000;
data = randi(2^bits_per_item, 1, data_count) - 1;
barker_code13 = [+1 +1 +1 +1 +1 -1 -1 +1 +1 -1 +1 -1 +1];
barker_symbs = (1+i)*barker_code13;
barker_symbs_up = upfirdn(barker_symbs, rrcos_filt, oversampl);

% data to complex numbers
bin_symbs = symbolify(data, bits_per_item, bits_per_symb);
points_per_axis = 2^(bits_per_symb/2);
qam_matrix = (-1:2/(points_per_axis-1):1) + (i*(-1:2/(points_per_axis-1):1)');
qam_lookup = qam_matrix(:)';
symbs = qam_lookup(bin_symbs + 1);

packet_header = [barker_symbs zeros(1, oversampl - 1)];
packet = [packet_header symbs];
[data_to_send, t_tx] = transmitter(packet, oversampl, lo_up, rrcos_filt, fs, flo);

% the channel is HERE
% do I even need HW?
the_num = randi(1000) + 10;
some_random_offset = zeros(1, the_num * lo_up * fs);
data_recieved = awgn([some_random_offset data_to_send], 30);
df = 1+1e-5;
some_phase = 2*pi*(rand()-0.5);

% plot fft
l = length(data_recieved);
freqdomain = fft(data_recieved)/l;
frequency_axis = linspace(0,fs/2,l/2+1);
PSD = abs(freqdomain(1:l/2+1));
plot(frequency_axis, PSD)
frequency_width = (1+rrcos_filt_beta)/lo_up/oversampl;
frequency_width_corr = ones(1, floor(l*frequency_width/fs/2));
slope_len = floor(frequency_width/2);
frequency_width_corr(1:slope_len) = linspace(0,1,slope_len);
frequency_width_corr(end-slope_len+1:end) = linspace(1,0,slope_len);
the_corr_in_f = xcorr(frequency_width_corr, PSD);
figure
plot(the_corr_in_f)
f_inds = find(the_corr_in_f>max(the_corr_in_f)/5);
lo_guess_frequency_index = floor((f_inds(1) + f_inds(end))/2);
flo_guess = frequency_axis(length(PSD)-lo_guess_frequency_index+floor(length(frequency_width_corr)/2));
lo_diff_f = flo_guess-flo

t_rx = 0:1/fs:((length(data_recieved)-1)/fs);
lo_i = sin(2*pi*flo_guess*t_rx+some_phase);
lo_q = cos(2*pi*flo_guess*t_rx+some_phase);

lo_remove_lp_filter = fir1(2*lo_up*fs,2/fs/lo_up);
rx_with_lo = data_recieved .* lo_i + 1j*data_recieved .* lo_q;

rx_high_sample = upfirdn(lo_remove_lp_filter, rx_with_lo, 1, lo_up*fs);
rx_high_sample_shift = rx_high_sample(2:end-1);

abs_corr = abs(xcorr(rx_high_sample_shift, barker_symbs_up));
[v, corr_max_index] = max(abs_corr);
packet_offset = corr_max_index-length(rx_high_sample_shift);
rx_high_sample_shift_packet = rx_high_sample_shift((packet_offset+1):end);

rx = upfirdn(rrcos_filt, rx_high_sample_shift_packet, 1, oversampl);
rx_shift = rx((rrcos_filt_len + 1):(end - rrcos_filt_len));

rx_without_header = rx_shift(length(packet_header)+1:end);
rx_header = rx_shift(1:length(barker_symbs));

symbs_error = barker_symbs./rx_header;
gain_comp = mean(abs(symbs_error));

header_phase_error = angle(symbs_error);
header_len = length(header_phase_error);
header_phase_error_fit = polyfit(1:header_len, header_phase_error, 1);
phase_error = header_phase_error_fit(2) + header_phase_error_fit(1) * header_len;

symbs_rx = rx_without_header * gain_comp;

rx_bin_symbs = zeros(size(symbs));
phase_error_coeff = 0.2;
measured_symbs = zeros(size(symbs_rx));
phase_errors = zeros(size(symbs_rx));
last_phase_error = header_phase_error_fit(2) + header_phase_error_fit(1) * (header_len - 1);

for k = 1:length(symbs_rx)
   measured_symb = symbs_rx(k)*exp(1i*phase_error);
   [dont_care, qam_i] = min(abs(measured_symb - qam_lookup));
   resolved_symb = qam_lookup(qam_i);
   new_phase_error = angle(resolved_symb / measured_symb);
   phase_error = phase_error_coeff * new_phase_error + phase_error;

   phase_errors(k) = phase_error - last_phase_error;
   last_phase_error = phase_error;
   measured_symbs(k) = measured_symb;
   rx_bin_symbs(k) = qam_i - 1;
   
end

rx_data = symbolify(rx_bin_symbs, bits_per_symb, bits_per_item);
rx_data = rx_data(1:length(data));
[bit_errors, ber] = biterr(rx_data,data)

figure
hold on
plot(measured_symbs, 'o');
plot(symbs, 'o')
grid on
axis([-1.5 1.5 -1.5 1.5])

