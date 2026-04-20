function bits = qpsk_demod(symbols)

symbols = symbols * sqrt(2);  % undo normalization

bits = zeros(2*length(symbols),1);

% Decision
bits(1:2:end) = real(symbols) < 0;
bits(2:2:end) = imag(symbols) < 0;

end