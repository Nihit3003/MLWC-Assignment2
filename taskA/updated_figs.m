clc;
clear;
close all;

addpath(genpath('winner2'));

%% PARAMETERS
N = 64;
cp_len = 16;
M = 4;
bitsPerSym = log2(M);

SNR_dB_range = -5:2:20;
numIter = 1000;

%% LOAD DNN
W1 = readNPY('W1.npy'); b1 = readNPY('b1.npy');
W2 = readNPY('W2.npy'); b2 = readNPY('b2.npy');
W3 = readNPY('W3.npy'); b3 = readNPY('b3.npy');

b1 = b1(:)'; b2 = b2(:)'; b3 = b3(:)';

%% ================= FUNCTION: CHANNEL =================
get_channel = @(N) generate_channel(N);

%% ================= FIG 1: 64 PILOTS =================
pilot_idx_64 = 1:N;
data_idx_64 = 1:N;   % all used

[BER_LS_64, BER_MMSE_64, BER_DNN_64] = ...
    run_sim(N, cp_len, pilot_idx_64, data_idx_64, SNR_dB_range, numIter, W1,W2,W3,b1,b2,b3);

figure;
semilogy(SNR_dB_range, BER_LS_64,'r-o','LineWidth',2); hold on;
semilogy(SNR_dB_range, BER_MMSE_64,'b-s','LineWidth',2);
semilogy(SNR_dB_range, BER_DNN_64,'g-d','LineWidth',2);
grid on;
title('Fig 1: BER vs SNR (64 Pilots, WINNER II)');
xlabel('SNR (dB)'); ylabel('BER');
legend('LS','MMSE','DNN');

%% ================= FIG 2: 8 PILOTS =================
pilot_idx_8 = 1:8:N;
data_idx_8 = setdiff(1:N, pilot_idx_8);

[BER_LS_8, BER_MMSE_8, BER_DNN_8] = ...
    run_sim(N, cp_len, pilot_idx_8, data_idx_8, SNR_dB_range, numIter, W1,W2,W3,b1,b2,b3);

figure;
semilogy(SNR_dB_range, BER_LS_8,'r-o','LineWidth',2); hold on;
semilogy(SNR_dB_range, BER_MMSE_8,'b-s','LineWidth',2);
semilogy(SNR_dB_range, BER_DNN_8,'g-d','LineWidth',2);
grid on;
title('Fig 2: BER vs SNR (8 Pilots, WINNER II)');
xlabel('SNR (dB)'); ylabel('BER');
legend('LS','MMSE','DNN');

%% ================= FIG 3: WITH vs WITHOUT CP =================
pilot_idx = 1:N;
data_idx = 1:N;

[BER_LS_CP, BER_MMSE_CP, BER_DNN_CP] = ...
    run_sim(N, cp_len, pilot_idx, data_idx, SNR_dB_range, numIter, W1,W2,W3,b1,b2,b3);

[BER_LS_noCP, BER_MMSE_noCP, BER_DNN_noCP] = ...
    run_sim(N, 0, pilot_idx, data_idx, SNR_dB_range, numIter, W1,W2,W3,b1,b2,b3);

figure;
semilogy(SNR_dB_range, BER_LS_CP,'r-o','LineWidth',2); hold on;
semilogy(SNR_dB_range, BER_MMSE_CP,'b-s','LineWidth',2);
semilogy(SNR_dB_range, BER_DNN_CP,'g-d','LineWidth',2);

semilogy(SNR_dB_range, BER_LS_noCP,'r--^','LineWidth',2);
semilogy(SNR_dB_range, BER_MMSE_noCP,'b--v','LineWidth',2);
semilogy(SNR_dB_range, BER_DNN_noCP,'g--x','LineWidth',2);

grid on;
title('Fig 3: With CP vs Without CP (WINNER II)');
xlabel('SNR (dB)'); ylabel('BER');

legend('LS CP','MMSE CP','DNN CP',...
       'LS noCP','MMSE noCP','DNN noCP');

%% ================= FUNCTION: MAIN SIM =================
function [BER_LS, BER_MMSE, BER_DNN] = run_sim(N, cp_len, pilot_idx, data_idx, SNR_range, numIter, W1,W2,W3,b1,b2,b3)

bitsPerSym = 2;

BER_LS = zeros(length(SNR_range),1);
BER_MMSE = zeros(length(SNR_range),1);
BER_DNN = zeros(length(SNR_range),1);

for snr_idx = 1:length(SNR_range)

    SNR_dB = SNR_range(snr_idx);
    noise_var = 1/(10^(SNR_dB/10));

    err_ls=0; err_mmse=0; err_dnn=0;
    total_bits=0;

    for iter = 1:numIter

        %% CHANNEL
        h = generate_channel(N);

        %% DATA
        X = zeros(N,1);
        X(pilot_idx) = 1;

        data_bits = randi([0 1], length(data_idx)*bitsPerSym,1);
        X(data_idx) = qpsk_mod(data_bits);

        %% TX
        tx = ifft(X)*sqrt(N);
        if cp_len > 0
            tx = [tx(end-cp_len+1:end); tx];
        end

        %% CHANNEL
        rx = conv(tx,h);
        rx = rx(1:length(tx));

        noise = sqrt(noise_var/2)*(randn(size(rx))+1j*randn(size(rx)));
        rx = rx + noise;

        %% RX
        if cp_len > 0
            rx = rx(cp_len+1:end);
        end
        Y = fft(rx)/sqrt(N);

        %% LS
        H_LS_p = Y(pilot_idx)./X(pilot_idx);
        H_LS = interp1(pilot_idx,H_LS_p,1:N,'linear','extrap').';

        %% MMSE
        L = length(h);
        pdp = abs(h).^2; pdp = pdp/sum(pdp);
        R = zeros(N,N);
        for k=1:N
            for l=1:N
                for n=1:L
                    R(k,l)=R(k,l)+pdp(n)*exp(-1j*2*pi*(k-l)*(n-1)/N);
                end
            end
        end
        H_MMSE = R*((R+noise_var*eye(N))\H_LS);

        %% DNN
        dnn_in = [real(H_LS).' imag(H_LS).'];

        a1 = max(0, dnn_in*W1 + b1);
        a2 = max(0, a1*W2 + b2);
        z3 = a2*W3 + b3;

        H_DNN = z3(1:N).' + 1j*z3(N+1:end).';

        %% EQUALIZE
        Xhat_LS = Y(data_idx)./H_LS(data_idx);
        Xhat_MMSE = Y(data_idx)./H_MMSE(data_idx);
        Xhat_DNN = Y(data_idx)./H_DNN(data_idx);

        %% DEMOD
        bits_ls = qpsk_demod(Xhat_LS);
        bits_mmse = qpsk_demod(Xhat_MMSE);
        bits_dnn = qpsk_demod(Xhat_DNN);

        err_ls = err_ls + sum(bits_ls~=data_bits);
        err_mmse = err_mmse + sum(bits_mmse~=data_bits);
        err_dnn = err_dnn + sum(bits_dnn~=data_bits);

        total_bits = total_bits + length(data_bits);

    end

    BER_LS(snr_idx)=err_ls/total_bits;
    BER_MMSE(snr_idx)=err_mmse/total_bits;
    BER_DNN(snr_idx)=err_dnn/total_bits;

end
end

%% ================= CHANNEL FUNCTION =================
function h = generate_channel(N)

try
    cfg = winner2.wimparset;
    layout = winner2.layoutparset;
    linkpar = winner2.linkparset;

    cfg.Scenario='A1';
    cfg.CenterFrequency=2.5e9;

    layout.Stations(1).Pos=[0;0;0];
    layout.Stations(2).Pos=[10;0;0];

    [H,~]=winner2.wim(cfg,layout,linkpar);
    h = squeeze(H(1,1,:,1));
    h = h(:);

catch
    % fallback
    L=16;
    pdp=exp(-(0:L-1)/2); pdp=pdp/sum(pdp);
    h=(randn(L,1)+1j*randn(L,1)).*sqrt(pdp(:)/2);
end

if length(h)>N
    h=h(1:N);
end

end