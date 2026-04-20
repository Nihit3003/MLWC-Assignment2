function tx = ofdm_tx(X, N, cp_len)
% OFDM_TX  OFDM modulator: IFFT + optional cyclic-prefix insertion.
%
%   tx = ofdm_tx(X, N, cp_len)
%
%   Inputs:
%     X      - Frequency-domain symbol vector (N x 1).
%     N      - FFT size.
%     cp_len - Cyclic-prefix length (set to 0 to disable).
%
%   Output:
%     tx     - Time-domain signal with CP prepended.

tx = ifft(X) * sqrt(N);          % normalised IFFT

if cp_len > 0
    tx = [tx(end-cp_len+1:end); tx];
end

end
