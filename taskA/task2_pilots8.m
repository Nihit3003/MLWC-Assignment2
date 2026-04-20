clc;
clear;
close all;

%% PARAMETERS
N = 64;
cp_len = 16;
M = 4;
bitsPerSym = log2(M);

pilot_idx = 1:8:N;
SNR_dB_range = -5:2:20;
numIter = 4000;

BER_LS = zeros(length(SNR_dB_range),1);
BER_MMSE = zeros(length(SNR_dB_range),1);

%% LOOP
for snr_idx = 1:length(SNR_dB_range)

    SNR_dB = SNR_dB_range(snr_idx);
    noise_var = 1/(10^(SNR_dB/10));   % ✅ correct SNR scaling

    err_ls = 0;
    err_mmse = 0;
    total_bits = 0;

    for iter = 1:numIter

        %% CHANNEL (MATCH PAPER)
        L = 4;   % 🔥 KEY FIX
        pdp = exp(-(0:L-1)/1.5);
        pdp = pdp / sum(pdp);

        h = (randn(L,1)+1j*randn(L,1)) .* sqrt(pdp(:)/2);

        %% DATA
        data_bits = randi([0 1], N*bitsPerSym, 1);
        X = qpsk_mod(data_bits);

        %% PILOTS
        X(pilot_idx) = 1;

        %% TX
        tx = ifft(X) * sqrt(N);   % ✅ IMPORTANT FIX
        tx_cp = [tx(end-cp_len+1:end); tx];

        %% CHANNEL
        rx = conv(tx_cp, h);
        rx = rx(1:length(tx_cp));

        %% NOISE
        noise = sqrt(noise_var/2)*(randn(size(rx))+1j*randn(size(rx)));
        rx = rx + noise;

        %% RX
        rx = rx(cp_len+1:end);
        Y = fft(rx) / sqrt(N);   % normalize back

        %% ================= LS =================
        H_LS_pilot = Y(pilot_idx) ./ X(pilot_idx);
        H_LS = interp1(pilot_idx, H_LS_pilot, 1:N, 'linear', 'extrap').';

        %% ================= MMSE =================
        R_HH = zeros(N,N);
        for k = 1:N
            for l = 1:N
                for n = 1:L
                    R_HH(k,l) = R_HH(k,l) + pdp(n)*exp(-1j*2*pi*(k-l)*(n-1)/N);
                end
            end
        end

        R_PP = R_HH(pilot_idx, pilot_idx);
        R_HP = R_HH(:, pilot_idx);

        H_MMSE = R_HP * ((R_PP + noise_var*eye(length(pilot_idx))) \ H_LS_pilot);

        %% ================= EQUALIZATION =================
        Xhat_LS = Y ./ H_LS;
        Xhat_MMSE = Y ./ H_MMSE;

        %% DEMOD
        bits_ls = qpsk_demod(Xhat_LS);
        bits_mmse = qpsk_demod(Xhat_MMSE);

        %% ERRORS
        err_ls = err_ls + sum(bits_ls ~= data_bits);
        err_mmse = err_mmse + sum(bits_mmse ~= data_bits);
        total_bits = total_bits + length(data_bits);

    end

    BER_LS(snr_idx) = err_ls / total_bits;
    BER_MMSE(snr_idx) = err_mmse / total_bits;

    fprintf("SNR=%d | LS=%.5f | MMSE=%.5f\n", ...
        SNR_dB, BER_LS(snr_idx), BER_MMSE(snr_idx));
end

%% PLOT
figure;
semilogy(SNR_dB_range, BER_LS, 'r-o','LineWidth',2);
hold on;
semilogy(SNR_dB_range, BER_MMSE, 'b-s','LineWidth',2);

grid on;
xlabel('SNR (dB)');
ylabel('BER');
legend('LS (8 pilots)','MMSE (8 pilots)','Location','southwest');
title('BER vs SNR (8 Pilots - TRUE PAPER MATCH)');