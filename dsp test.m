clc, close all, clear all


bits_per_symb = 6;
bits_per_item = 8;
data = randi(2^bits_per_item - 1, 1, 10000);

bin_symbs = symbolify(data, bits_per_item, bits_per_symb);

points_per_axis = 2^(bits_per_symb/2);
qam16_matrix = (-1:2/(points_per_axis-1):1) + (i*(-1:2/(points_per_axis-1):1)');
qam16_lookup = reshape(qam16_matrix.',1,[]);

symbs = qam16_lookup(bin_symbs+1);
oversampl = 9;
fs = 10;
[data_to_send, t_tx] = transmitter(symbs, oversampl, fs);
data_recieved = data_to_send; % awgn(data_to_send, 10); % do I even need HW?

df = 1;
fs = fs*df; % we have a frequency diff

lo_i = sin(2*pi*fs*t_tx);
lo_q = cos(2*pi*fs*t_tx);

b = fir1(10,1/fs);

rxi_high_sample = filtfilt(b, 1, data_recieved .* lo_i);
rxq_high_sample = filtfilt(b, 1, data_recieved .* lo_q);

raised_cos_filter = rcosdesign(0.1,6,oversampl);

rxi = resample(filtfilt(raised_cos_filter, 1, rxi_high_sample), 1, fs);
rxq = resample(filtfilt(raised_cos_filter, 1, rxq_high_sample), 1, fs);

sample_points = floor((1+oversampl)/2):oversampl:length(rxi);
i = rxi(sample_points);
q = rxq(sample_points);
scale_factor = 0.8;
symbs_rx = (i/max(i) + 1j*q/max(q))/scale_factor;

eyediagram(rxi, 3*oversampl);

figure
hold on
plot(symbs_rx, 'o');
plot(symbs, 'o')
grid on
axis([-1.5 1.5 -1.5 1.5])

rx_bin_symbs = zeros(size(symbs_rx));
for k = 1:length(symbs)
   [dont_care, qam_i] = min(abs(symbs_rx(k)-qam16_lookup));
   rx_bin_symbs(k) = qam_i - 1;
end

symberrors = sum(bin_symbs-rx_bin_symbs ~= 0)
rx_data = symbolify(rx_bin_symbs, bits_per_symb, bits_per_item);
figure
hold on
%plot(rx_data-data)
bit_errors = biterr(rx_data,data)

