function [txdata, t_up] = transmitter(symbs, oversampl, fs)
%TRANSMITTER Summary of this function goes here
%   Detailed explanation goes here
    
    symbs_up = zeros(1, oversampl * length(symbs));
    for i = 1:length(symbs)
        symbs_up((i-1)*oversampl+1:i*oversampl) = symbs(i);
    end

    symbs_q = imag(symbs_up);
    symbs_i = real(symbs_up);
    %eyediagram(symbs_i, 3*oversampl);
    t = linspace(0, length(symbs), length(symbs_up));

    raised_cos_filter = rcosdesign(0.1,6,oversampl);

    symbs_i_pulse_shaped = filtfilt(raised_cos_filter, 1, symbs_i);
    symbs_q_pulse_shaped = filtfilt(raised_cos_filter, 1, symbs_q);
    %eyediagram(symbs_i_pulse_shaped, 3*oversampl);
    
    symbs_i_pulse_shaped_up = resample(symbs_i_pulse_shaped, fs, 1);
    symbs_q_pulse_shaped_up = resample(symbs_q_pulse_shaped, fs, 1);
    t_up = linspace(min(t), max(t), length(symbs_q_pulse_shaped_up));
    
    lo_in = sin(2*pi*fs*t_up);
    lo_out = cos(2*pi*fs*t_up);

    txdata = symbs_i_pulse_shaped_up.*lo_in + symbs_q_pulse_shaped_up.*lo_out;
end

