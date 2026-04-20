function [H_LS, H_MMSE, H_DNN] = estimate_channel(Y, X, pilot_idx, N, h, noise_var, dnn)
% ESTIMATE_CHANNEL  LS, MMSE, and DNN channel estimators.
%
%   [H_LS, H_MMSE, H_DNN] = estimate_channel(Y, X, pilot_idx, N, h, noise_var, dnn)
%
%   Inputs:
%     Y         - Received frequency-domain vector (N x 1).
%     X         - Transmitted frequency-domain vector (N x 1) — used only
%                 at pilot positions.
%     pilot_idx - Indices of pilot subcarriers.
%     N         - FFT size (= number of subcarriers).
%     h         - Channel impulse response (used to build MMSE covariance).
%     noise_var - Noise variance (linear scale).
%     dnn       - Struct with fields W1,b1,W2,b2,W3,b3 for the DNN
%                 estimator.  Pass [] to skip DNN.
%
%   Outputs:
%     H_LS   - LS channel estimate  (N x 1).
%     H_MMSE - MMSE channel estimate (N x 1).
%     H_DNN  - DNN channel estimate  (N x 1), or [] if dnn is empty.

%% ---- LS ----------------------------------------------------------------
H_LS_pilot = Y(pilot_idx) ./ X(pilot_idx);

if length(pilot_idx) == N
    % All subcarriers are pilots — no interpolation needed
    H_LS = H_LS_pilot;
else
    H_LS = interp1(pilot_idx, H_LS_pilot, 1:N, 'linear', 'extrap').';
end

%% ---- MMSE --------------------------------------------------------------
L   = length(h);
pdp = abs(h).^2;
pdp = pdp / sum(pdp);

% Build full N×N channel correlation matrix
R = zeros(N, N);
for k = 1:N
    for l = 1:N
        for n = 1:L
            R(k,l) = R(k,l) + pdp(n) * exp(-1j*2*pi*(k-l)*(n-1)/N);
        end
    end
end

if length(pilot_idx) == N
    % Full-pilot case: standard Wiener filter
    H_MMSE = R * ((R + noise_var*eye(N)) \ H_LS);
else
    % Sparse-pilot case: cross-covariance between all and pilot subcarriers
    R_PP   = R(pilot_idx, pilot_idx);
    R_HP   = R(:, pilot_idx);
    H_MMSE = R_HP * ((R_PP + noise_var*eye(length(pilot_idx))) \ H_LS_pilot);
end

%% ---- DNN ---------------------------------------------------------------
if isempty(dnn)
    H_DNN = [];
    return;
end

dnn_in = [real(H_LS).' imag(H_LS).'];

a1    = max(0, dnn_in * dnn.W1 + dnn.b1);
a2    = max(0, a1     * dnn.W2 + dnn.b2);
z3    = a2 * dnn.W3 + dnn.b3;

H_DNN = z3(1:N).' + 1j*z3(N+1:end).';

end
