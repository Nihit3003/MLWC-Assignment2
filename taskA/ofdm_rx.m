function Y = ofdm_rx(rx, N, cp_len)
% OFDM_RX  OFDM demodulator: CP removal + FFT.
%
%   Y = ofdm_rx(rx, N, cp_len)
%
%   Inputs:
%     rx     - Received time-domain signal.
%     N      - FFT size.
%     cp_len - Cyclic-prefix length that was added at the transmitter
%              (set to 0 if no CP was used).
%
%   Output:
%     Y      - Frequency-domain received vector (N x 1).

if cp_len > 0
    rx = rx(cp_len+1:end);
end

Y = fft(rx) / sqrt(N);           % normalised FFT

end
