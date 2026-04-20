function symbols = qpsk_mod(bits)
% QPSK_MOD  Gray-coded QPSK modulation with unit average power.
%
%   symbols = qpsk_mod(bits)
%
%   Input:
%     bits    - Column vector of bits (length must be even).
%
%   Output:
%     symbols - Complex QPSK symbols, normalised to unit power.

bits    = reshape(bits, 2, []).';
symbols = (1 - 2*bits(:,1)) + 1j*(1 - 2*bits(:,2));
symbols = symbols / sqrt(2);

end
