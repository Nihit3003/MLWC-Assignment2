function [BER_LS, BER_MMSE, BER_DNN] = run_sim(cfg)
% RUN_SIM  OFDM BER simulation with LS, MMSE, and (optionally) DNN channel
%          estimation over a range of SNR values.
%
%   [BER_LS, BER_MMSE, BER_DNN] = run_sim(cfg)
%
%   cfg fields (all required unless noted):
%     .N          - FFT size / number of subcarriers.
%     .cp_len     - Cyclic-prefix length (0 = no CP).
%     .pilot_idx  - Pilot subcarrier indices (row or column vector).
%     .SNR_dB_vec - Vector of SNR values in dB.
%     .numIter    - Monte-Carlo iterations per SNR point.
%     .dnn        - DNN weight struct (fields W1,b1,W2,b2,W3,b3).
%                   Pass [] to skip DNN.
%     .use_winner2- (optional) true/false; default true.
%
%   Outputs:
%     BER_LS, BER_MMSE, BER_DNN - BER vectors (length = numel(SNR_dB_vec)).

bitsPerSym   = 2;                                   % QPSK
N            = cfg.N;
cp_len       = cfg.cp_len;
pilot_idx    = cfg.pilot_idx(:);
SNR_dB_vec   = cfg.SNR_dB_vec;
numIter      = cfg.numIter;
dnn          = cfg.dnn;
use_winner2  = isfield(cfg,'use_winner2') && cfg.use_winner2;

data_idx  = setdiff(1:N, pilot_idx);
all_pilot = isempty(data_idx);                      % 64-pilot mode

nSNR   = numel(SNR_dB_vec);
BER_LS   = zeros(nSNR, 1);
BER_MMSE = zeros(nSNR, 1);
BER_DNN  = zeros(nSNR, 1);

for si = 1:nSNR

    SNR_dB    = SNR_dB_vec(si);
    noise_var = 1 / 10^(SNR_dB/10);

    err_ls = 0;  err_mmse = 0;  err_dnn = 0;
    total_bits = 0;

    for iter = 1:numIter

        %% ---- Channel --------------------------------------------------
        h = generate_channel(N, use_winner2);

        %% ---- Build frequency-domain frame -----------------------------
        X = zeros(N, 1);
        X(pilot_idx) = 1;                           % unit pilots

        if all_pilot
            % Pilots fill all subcarriers; data is re-encoded over them
            data_bits = randi([0 1], N*bitsPerSym, 1);
            X = qpsk_mod(data_bits);                % overwrite pilots
            X_data_tx = X;
        else
            data_bits  = randi([0 1], numel(data_idx)*bitsPerSym, 1);
            X(data_idx) = qpsk_mod(data_bits);
        end

        %% ---- Transmit -------------------------------------------------
        tx = ofdm_tx(X, N, cp_len);

        %% ---- Multipath channel + noise --------------------------------
        rx = conv(tx, h);
        rx = rx(1:length(tx));
        rx = rx + sqrt(noise_var/2) * (randn(size(rx)) + 1j*randn(size(rx)));

        %% ---- Receive --------------------------------------------------
        Y = ofdm_rx(rx, N, cp_len);

        %% ---- Channel estimation ---------------------------------------
        [H_LS, H_MMSE, H_DNN] = estimate_channel(Y, X, pilot_idx, N, h, noise_var, dnn);

        %% ---- Equalize + demod ----------------------------------------
        if all_pilot
            % Full-pilot frame: equalize all subcarriers
            bits_ls   = qpsk_demod(Y ./ H_LS);
            bits_mmse = qpsk_demod(Y ./ H_MMSE);
            if ~isempty(H_DNN)
                bits_dnn = qpsk_demod(Y ./ H_DNN);
            end
        else
            bits_ls   = qpsk_demod(Y(data_idx) ./ H_LS(data_idx));
            bits_mmse = qpsk_demod(Y(data_idx) ./ H_MMSE(data_idx));
            if ~isempty(H_DNN)
                bits_dnn = qpsk_demod(Y(data_idx) ./ H_DNN(data_idx));
            end
        end

        %% ---- Accumulate errors ----------------------------------------
        err_ls   = err_ls   + sum(bits_ls   ~= data_bits);
        err_mmse = err_mmse + sum(bits_mmse ~= data_bits);
        if ~isempty(H_DNN)
            err_dnn = err_dnn + sum(bits_dnn ~= data_bits);
        end
        total_bits = total_bits + numel(data_bits);

    end

    BER_LS(si)   = err_ls   / total_bits;
    BER_MMSE(si) = err_mmse / total_bits;
    if ~isempty(dnn)
        BER_DNN(si) = err_dnn / total_bits;
    end

    fprintf('SNR=%3d dB | LS=%.4f | MMSE=%.4f | DNN=%.4f\n', ...
        SNR_dB, BER_LS(si), BER_MMSE(si), BER_DNN(si));
end

end
