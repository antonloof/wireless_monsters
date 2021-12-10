clc, close all


N = 1000;
samples = (randi(4, 1, N) + 1i*randi(4, 1, N) - 1 - 1i) * 2 / 3 - 1 - 1i;
%samples = 2*randi(2, 1, N) + 2i*randi(2, 1, N) - 3 - 3i;
t = 1:N;
high_samp_up = 100;
t_up = (1:0.5:N-0.5);
t_up = t_up - 0.1;

sampled = interp1(t, samples, t_up, 'spline');

resolved_samples = zeros(size(samples));

gardner_loop_ki = 0.1;
gardner_loop_kp = 0.05;
gardner_loop_i = 0;

dk = 0;
last_dk = 0;
iter = 4:N-3;
gardnered_samples = zeros(1, length(iter));
gardner_loop = zeros(1, length(iter));

for i = iter
    prev_sample = interp1([-1 0 1], sampled(2*i-3:2*i-1), last_dk);
    half_sample = interp1([-1 0 1], sampled(2*i-2:2*i), last_dk);
    this_sample = interp1([-1 0 1], sampled(2*i-1:2*i+1), dk);
    
    gardner_error_r = real(prev_sample - this_sample) * real(half_sample);
    gardner_error_i = imag(prev_sample - this_sample) * imag(half_sample);
    gardner_error = gardner_error_i + gardner_error_r;   
    gardnered_samples(i) = this_sample;
    
    gardner_loop_p = gardner_loop_kp * gardner_error;
    gardner_loop_i = gardner_loop_ki * gardner_error + gardner_loop_i;
    gardner_loop_out = gardner_loop_p + gardner_loop_i;
    gardner_loop(i) = gardner_loop_out;
    last_dk = dk;
    dk = gardner_loop_out;
end

figure
hold on
plot(gardner_loop)


figure
hold on
gardner_locked_at = 500;
gard = gardnered_samples(gardner_locked_at:end);
expected = samples(gardner_locked_at:end-1);
% rude gain control
% gard = mean(abs(expected(2:14) ./ gard(1:13))) * gard;
stupid = sampled((gardner_locked_at*2)+1:2:end);
plot(gard, 'o')
plot(stupid, '*')
plot(expected, 'o')

%gard_score = sum(abs(gard-expected))/length(gard)
%stupid_score = sum(abs(stupid-expected))/length(stupid)