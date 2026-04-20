function bits = qpsk_demod(symbols)
% QPSK_DEMOD  Hard-decision QPSK demodulator.
%
%   bits = qpsk_demod(symbols)
%
%   Input:
%     symbols - Complex received symbols (any length).
%
%   Output:
%     bits    - Recovered bit vector (length = 2 * numel(symbols)).

symbols       = symbols * sqrt(2);          % undo power normalisation
bits          = zeros(2*numel(symbols), 1);
bits(1:2:end) = real(symbols) < 0;
bits(2:2:end) = imag(symbols) < 0;

end
