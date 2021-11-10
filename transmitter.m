function [txdata, t_up] = transmitter(symbs, oversampl, flo, rrcos_filt)
%TRANSMITTER Summary of this function goes here
%   Detailed explanation goes here
    symbs_pulse_shaped = upfirdn(symbs, rrcos_filt, oversampl, 1);
    symbs_pulse_shaped_up = resample(symbs_pulse_shaped, flo, 1);
    
    t_up = linspace(1, length(symbs), length(symbs_pulse_shaped_up));
    lo_in = sin(2*pi*flo*t_up);
    lo_out = cos(2*pi*flo*t_up);
    hold on
    %plot(linspace(1, length(symbs), length(symbs_pulse_shaped)), real(symbs_pulse_shaped))
    plot(t_up, lo_in);
    txdata = real(symbs_pulse_shaped_up).*lo_in + imag(symbs_pulse_shaped_up).*lo_out;
end

