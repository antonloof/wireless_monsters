function [txdata, t_up] = transmitter(symbs, oversampl, fs, rrcos_filt)
%TRANSMITTER Summary of this function goes here
%   Detailed explanation goes here
    symbs_i_pulse_shaped = upfirdn(real(symbs), rrcos_filt, oversampl, 1);
    symbs_q_pulse_shaped = upfirdn(imag(symbs), rrcos_filt, oversampl, 1);
    t = linspace(0, length(symbs), length(symbs_i_pulse_shaped));
    symbs_i_pulse_shaped_up = resample(symbs_i_pulse_shaped, fs, 1);
    symbs_q_pulse_shaped_up = resample(symbs_q_pulse_shaped, fs, 1);
    t_up = linspace(min(t), max(t), length(symbs_q_pulse_shaped_up));
    lo_in = sin(2*pi*fs*t_up);
    lo_out = cos(2*pi*fs*t_up);

    txdata = symbs_i_pulse_shaped_up.*lo_in + symbs_q_pulse_shaped_up.*lo_out;
end

