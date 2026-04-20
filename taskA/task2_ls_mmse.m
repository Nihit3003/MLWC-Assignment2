clc;
clear;
close all;

%% PARAMETERS
N = 64;
cp_len = 16;
M = 4;
bitsPerSym = log2(M);

SNR_dB_range = -5:2:20;
numIter = 1500;

BER_LS = zeros(length(SNR_dB_range),1);
BER_MMSE = zeros(length(SNR_dB_range),1);

%% LOOP
for snr_idx = 1:length(SNR_dB_range)

    SNR_dB = SNR_dB_range(snr_idx);
    noise_var = 10^(-SNR_dB/10);

    err_ls = 0;
    err_mmse = 0;
    total_bits = 0;

    for iter = 1:numIter

        %% CHANNEL
        L = 12;
        pdp = exp(-(0:L-1)/2);
        pdp = pdp / sum(pdp);
        h = (randn(L,1)+1j*randn(L,1)) .* sqrt(pdp(:)/2);

        %% DATA
        data_bits = randi([0 1], N*bitsPerSym, 1);
        X_data = qpsk_mod(data_bits);

        %% PILOTS (ALL SUBCARRIERS)
        X_pilot = ones(N,1);   % known pilots

        %% TRANSMIT SIGNAL (pilots only for estimation)
        tx = ifft(X_pilot);
        tx = tx / sqrt(mean(abs(tx).^2));
        tx_cp = [tx(end-cp_len+1:end); tx];

        %% CHANNEL + NOISE
        rx = conv(tx_cp, h);
        rx = rx(1:length(tx_cp));
        noise = sqrt(noise_var/2)*(randn(size(rx))+1j*randn(size(rx)));
        rx = rx + noise;

        %% RX FFT
        rx_no_cp = rx(cp_len+1:end);
        Y_pilot = fft(rx_no_cp);

        %% TRUE CHANNEL
        H_true = fft([h; zeros(N-length(h),1)]);

        %% LS ESTIMATION (CORRECT)
        H_LS = Y_pilot ./ X_pilot;

        %% MMSE (Wiener)
        R_HH = zeros(N,N);
        for m=1:N
            for n=1:N
                R_HH(m,n) = exp(-abs(m-n)/2);
            end
        end

        H_MMSE = R_HH * ((R_HH + noise_var*eye(N)) \ H_LS);

        %% ================= DATA TRANSMISSION =================
        tx_data = ifft(X_data);
        tx_data = tx_data / sqrt(mean(abs(tx_data).^2));
        tx_data_cp = [tx_data(end-cp_len+1:end); tx_data];

        rx_data = conv(tx_data_cp, h);
        rx_data = rx_data(1:length(tx_data_cp));
        noise = sqrt(noise_var/2)*(randn(size(rx_data))+1j*randn(size(rx_data)));
        rx_data = rx_data + noise;

        rx_data = rx_data(cp_len+1:end);
        Y_data = fft(rx_data);

        %% EQUALIZATION
        Xhat_LS = Y_data ./ H_LS;
        Xhat_MMSE = Y_data ./ H_MMSE;

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
legend('LS','MMSE','Location','southwest');
title('BER vs SNR (64 Pilots - Correct Paper Match)');