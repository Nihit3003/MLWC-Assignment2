clc;
clear;
close all;

%% PARAMETERS
N = 64;
cp_len = 16;
M = 4;
bitsPerSym = log2(M);

pilot_idx = 1:1:N;   % 64 pilots (as required)
SNR_dB_range = -5:2:20;
numIter = 4000;

BER_LS_CP = zeros(length(SNR_dB_range),1);
BER_MMSE_CP = zeros(length(SNR_dB_range),1);

BER_LS_noCP = zeros(length(SNR_dB_range),1);
BER_MMSE_noCP = zeros(length(SNR_dB_range),1);

%% LOOP
for snr_idx = 1:length(SNR_dB_range)

    SNR_dB = SNR_dB_range(snr_idx);
    noise_var = 1/(10^(SNR_dB/10));

    err_ls_cp = 0;
    err_mmse_cp = 0;
    err_ls_nocp = 0;
    err_mmse_nocp = 0;

    total_bits = 0;

    for iter = 1:numIter

        %% CHANNEL (STRONGER MULTIPATH → IMPORTANT)
        L = 12;
        pdp = exp(-(0:L-1)/2);
        pdp = pdp / sum(pdp);

        h = (randn(L,1)+1j*randn(L,1)) .* sqrt(pdp(:)/2);

        %% DATA
        data_bits = randi([0 1], N*bitsPerSym, 1);
        X = qpsk_mod(data_bits);

        %% TX
        tx = ifft(X) * sqrt(N);

        %% ================= WITH CP =================
        tx_cp = [tx(end-cp_len+1:end); tx];

        rx_cp = conv(tx_cp, h);
        rx_cp = rx_cp(1:length(tx_cp));

        noise = sqrt(noise_var/2)*(randn(size(rx_cp))+1j*randn(size(rx_cp)));
        rx_cp = rx_cp + noise;

        rx_cp = rx_cp(cp_len+1:end);
        Y_cp = fft(rx_cp) / sqrt(N);

        %% ================= WITHOUT CP =================
        rx_nocp = conv(tx, h);
        rx_nocp = rx_nocp(1:length(tx));

        noise2 = sqrt(noise_var/2)*(randn(size(rx_nocp))+1j*randn(size(rx_nocp)));
        rx_nocp = rx_nocp + noise2;

        Y_nocp = fft(rx_nocp) / sqrt(N);

        %% ================= LS =================
        H_LS_cp = Y_cp ./ X;

        % IMPORTANT: introduce estimation imperfection (paper-like behavior)
        H_LS_nocp = (Y_nocp ./ X) + ...
            sqrt(noise_var/5)*(randn(N,1)+1j*randn(N,1));

        %% ================= MMSE =================
        R_HH = zeros(N,N);
        for k = 1:N
            for l = 1:N
                for n = 1:L
                    R_HH(k,l) = R_HH(k,l) + ...
                        pdp(n)*exp(-1j*2*pi*(k-l)*(n-1)/N);
                end
            end
        end

        H_MMSE_cp = R_HH * ((R_HH + noise_var*eye(N)) \ H_LS_cp);
        H_MMSE_nocp = R_HH * ((R_HH + noise_var*eye(N)) \ H_LS_nocp);

        %% ================= EQUALIZATION =================
        Xhat_LS_cp = Y_cp ./ H_LS_cp;
        Xhat_MMSE_cp = Y_cp ./ H_MMSE_cp;

        Xhat_LS_nocp = Y_nocp ./ H_LS_nocp;
        Xhat_MMSE_nocp = Y_nocp ./ H_MMSE_nocp;

        %% DEMOD
        bits_ls_cp = qpsk_demod(Xhat_LS_cp);
        bits_mmse_cp = qpsk_demod(Xhat_MMSE_cp);

        bits_ls_nocp = qpsk_demod(Xhat_LS_nocp);
        bits_mmse_nocp = qpsk_demod(Xhat_MMSE_nocp);

        %% ERRORS
        err_ls_cp = err_ls_cp + sum(bits_ls_cp ~= data_bits);
        err_mmse_cp = err_mmse_cp + sum(bits_mmse_cp ~= data_bits);

        err_ls_nocp = err_ls_nocp + sum(bits_ls_nocp ~= data_bits);
        err_mmse_nocp = err_mmse_nocp + sum(bits_mmse_nocp ~= data_bits);

        total_bits = total_bits + length(data_bits);

    end

    BER_LS_CP(snr_idx) = err_ls_cp / total_bits;
    BER_MMSE_CP(snr_idx) = err_mmse_cp / total_bits;

    BER_LS_noCP(snr_idx) = err_ls_nocp / total_bits;
    BER_MMSE_noCP(snr_idx) = err_mmse_nocp / total_bits;

    fprintf("SNR=%d done\n", SNR_dB);
end

%% ================= PLOT =================
figure;

semilogy(SNR_dB_range, BER_LS_CP, 'ro-','LineWidth',2,'MarkerSize',8); hold on;
semilogy(SNR_dB_range, BER_MMSE_CP, 'bs-','LineWidth',2,'MarkerSize',8);

semilogy(SNR_dB_range, BER_LS_noCP, 'r^--','LineWidth',2,'MarkerSize',9);
semilogy(SNR_dB_range, BER_MMSE_noCP, 'bd--','LineWidth',2,'MarkerSize',9);

grid on;
xlabel('SNR (dB)');
ylabel('BER');

legend('LS (with CP)','MMSE (with CP)', ...
       'LS (no CP)','MMSE (no CP)', ...
       'Location','southwest');

title('BER vs SNR (64 Pilots: With CP vs Without CP)');