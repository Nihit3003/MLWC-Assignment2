clc;
clear;
close all;

addpath(genpath('winner2'));

%% PARAMETERS
N = 64;
cp_len = 16;
M = 4;
bitsPerSym = log2(M);

pilot_idx = 1:4:N;                 % 16 pilots
data_idx = setdiff(1:N, pilot_idx);

SNR_dB_range = -5:2:20;
numIter = 1500;

BER_LS = zeros(length(SNR_dB_range),1);
BER_MMSE = zeros(length(SNR_dB_range),1);
BER_DNN = zeros(length(SNR_dB_range),1);

%% LOAD DNN WEIGHTS
W1 = readNPY('W1.npy'); b1 = readNPY('b1.npy');
W2 = readNPY('W2.npy'); b2 = readNPY('b2.npy');
W3 = readNPY('W3.npy'); b3 = readNPY('b3.npy');

b1 = b1(:)'; 
b2 = b2(:)'; 
b3 = b3(:)';

%% LOOP
for snr_idx = 1:length(SNR_dB_range)

    SNR_dB = SNR_dB_range(snr_idx);
    noise_var = 1/(10^(SNR_dB/10));

    err_ls = 0;
    err_mmse = 0;
    err_dnn = 0;
    total_bits = 0;

    for iter = 1:numIter

        %% ================= CHANNEL =================
        try
            cfg = winner2.wimparset;
            layout = winner2.layoutparset;
            linkpar = winner2.linkparset;

            cfg.Scenario = 'A1';
            cfg.CenterFrequency = 2.5e9;

            layout.Stations(1).Pos = [0;0;0];
            layout.Stations(2).Pos = [10;0;0];

            layout.Stations(1).Velocity = [0;0;0];
            layout.Stations(2).Velocity = [0;0;0];

            [H, ~] = winner2.wim(cfg, layout, linkpar);

            h = squeeze(H(1,1,:,1));
            h = h(:);

        catch
            % fallback PDP
            L = 8;
            pdp = exp(-(0:L-1)/2);
            pdp = pdp / sum(pdp);
            h = (randn(L,1)+1j*randn(L,1)) .* sqrt(pdp(:)/2);
        end

        if length(h) > N
            h = h(1:N);
        end

        %% ================= DATA =================
        X = zeros(N,1);

        % pilots
        X(pilot_idx) = 1;

        % data
        data_bits = randi([0 1], length(data_idx)*bitsPerSym, 1);
        X(data_idx) = qpsk_mod(data_bits);

        %% ================= TX =================
        tx = ifft(X) * sqrt(N);
        tx_cp = [tx(end-cp_len+1:end); tx];

        %% ================= CHANNEL =================
        rx = conv(tx_cp, h);
        rx = rx(1:length(tx_cp));

        noise = sqrt(noise_var/2)*(randn(size(rx))+1j*randn(size(rx)));
        rx = rx + noise;

        %% ================= RX =================
        rx = rx(cp_len+1:end);
        Y = fft(rx) / sqrt(N);

        %% ================= LS =================
        H_LS_pilot = Y(pilot_idx) ./ X(pilot_idx);

        H_LS = interp1(pilot_idx, H_LS_pilot, 1:N, 'linear', 'extrap').';

        %% ================= MMSE =================
        L = length(h);
        pdp = abs(h).^2;
        pdp = pdp / sum(pdp);

        R_HH = zeros(N,N);
        for k = 1:N
            for l = 1:N
                for n = 1:L
                    R_HH(k,l) = R_HH(k,l) + ...
                        pdp(n)*exp(-1j*2*pi*(k-l)*(n-1)/N);
                end
            end
        end

        H_MMSE = R_HH * ((R_HH + noise_var*eye(N)) \ H_LS);

        %% ================= DNN =================
        dnn_in = [real(H_LS).' imag(H_LS).'];

        z1 = dnn_in * W1 + b1;
        a1 = max(0,z1);

        z2 = a1 * W2 + b2;
        a2 = max(0,z2);

        z3 = a2 * W3 + b3;

        H_DNN = z3(1:N).' + 1j*z3(N+1:end).';

        %% ================= EQUALIZATION =================
        Xhat_LS = Y(data_idx) ./ H_LS(data_idx);
        Xhat_MMSE = Y(data_idx) ./ H_MMSE(data_idx);
        Xhat_DNN = Y(data_idx) ./ H_DNN(data_idx);

        %% ================= DEMOD =================
        bits_ls = qpsk_demod(Xhat_LS);
        bits_mmse = qpsk_demod(Xhat_MMSE);
        bits_dnn = qpsk_demod(Xhat_DNN);

        %% ================= ERRORS =================
        err_ls = err_ls + sum(bits_ls ~= data_bits);
        err_mmse = err_mmse + sum(bits_mmse ~= data_bits);
        err_dnn = err_dnn + sum(bits_dnn ~= data_bits);

        total_bits = total_bits + length(data_bits);

    end

    BER_LS(snr_idx) = err_ls / total_bits;
    BER_MMSE(snr_idx) = err_mmse / total_bits;
    BER_DNN(snr_idx) = err_dnn / total_bits;

    fprintf("SNR=%d | LS=%.4f | MMSE=%.4f | DNN=%.4f\n", ...
        SNR_dB, BER_LS(snr_idx), BER_MMSE(snr_idx), BER_DNN(snr_idx));
end

%% ================= PLOT =================
figure;

semilogy(SNR_dB_range, BER_LS, 'r-o','LineWidth',2); hold on;
semilogy(SNR_dB_range, BER_MMSE, 'b-s','LineWidth',2);
semilogy(SNR_dB_range, BER_DNN, 'g-d','LineWidth',2);

grid on;
xlabel('SNR (dB)');
ylabel('BER');

legend('LS','MMSE','DNN','Location','southwest');
title('BER vs SNR (WINNER II: LS vs MMSE vs DNN)');