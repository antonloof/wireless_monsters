function [symbs] = symbolify(data, bits_per_item, bits_per_symb)
%SYMBOLIFY Summary of this function goes here
%   Detailed explanation goes here
    total_bits = ceil(bits_per_item * length(data) / bits_per_symb) * bits_per_symb;
    symbs = zeros(1, ceil(total_bits / bits_per_symb));

    carry = data(1) - 0;
    next_symb_i = 1;
    next_word_i = 2;
    bits_long = bits_per_item;
    while next_word_i <= length(data) + 1
        if bits_long < bits_per_symb
            if next_word_i <= length(data)
                carry = bitshift(carry, bits_per_item) + data(next_word_i);
                bits_long = bits_long + bits_per_item;
            end
            next_word_i = next_word_i + 1;
        end
        if bits_long >= bits_per_symb
            symbs(next_symb_i) = bitshift(carry, -(bits_long - bits_per_symb));
            carry = carry - bitshift(symbs(next_symb_i), bits_long - bits_per_symb);
            next_symb_i = next_symb_i + 1;
            bits_long = bits_long - bits_per_symb;
        end
    end
end

