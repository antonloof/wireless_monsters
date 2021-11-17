function [txdata, t_up] = transmitter(symbs, oversampl, lo_up, rrcos_filt, fs, flo)
%TRANSMITTER Summary of this function goes here
%   Detailed explanation goes here
    symbs_pulse_shaped = upfirdn(symbs, rrcos_filt, oversampl, 1);
    symbs_pulse_shaped_up = resample(symbs_pulse_shaped, fs*lo_up, 1);
    t_up = 0:1/fs:((length(symbs_pulse_shaped_up)-1)/fs);
    lo_in = sin(2*pi*flo*t_up);
    lo_out = cos(2*pi*flo*t_up);

    txdata = real(symbs_pulse_shaped_up).*lo_in + imag(symbs_pulse_shaped_up).*lo_out;
end

