function symbols = qpsk_mod(bits)

bits = reshape(bits, 2, []).';

% Gray-coded QPSK
symbols = (1 - 2*bits(:,1)) + 1j*(1 - 2*bits(:,2));

% Normalize power
symbols = symbols / sqrt(2);

end