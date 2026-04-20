clc;
clear;
close all;

%% ================= PARAMETERS =================
N = 64;              % Subcarriers
cp_len = 16;         % Cyclic Prefix
M = 4;               % QPSK
numSymbols = 1000;   % Number of OFDM symbols
SNR_dB = 10;         % Example SNR

%% ================= TRANSMITTER =================

% Generate random bits
bits = randi([0 1], numSymbols * N * log2(M), 1);

% QPSK Modulation
symbols = qpsk_mod(bits);

% Reshape into OFDM symbols
symbols_mat = reshape(symbols, N, []);

% IFFT (convert to time domain)
tx_ifft = ifft(symbols_mat, N);

% Add Cyclic Prefix
tx_cp = [tx_ifft(end-cp_len+1:end, :); tx_ifft];

% Serialize
tx_signal = tx_cp(:);

%% ================= CHANNEL =================

L = 5; % number of paths
h = (randn(L,1) + 1j*randn(L,1)) / sqrt(2*L);

% Convolution
rx_signal = conv(tx_signal, h);

% Trim to original length
rx_signal = rx_signal(1:length(tx_signal));

% Add noise
rx_signal = awgn(rx_signal, SNR_dB, 'measured');

%% ================= RECEIVER =================

% Reshape back
rx_cp = reshape(rx_signal, N + cp_len, []);

% Remove CP
rx_no_cp = rx_cp(cp_len+1:end, :);

% FFT (back to frequency domain)
rx_fft = fft(rx_no_cp, N);

%% ================= IDEAL CHANNEL (FOR TASK 1 ONLY) =================
% (We use perfect channel knowledge for now — estimation comes in Task 2)

H = fft(h, N);                % 64 x 1
H = repmat(H, 1, size(rx_fft,2));  % match dimensions

rx_eq = rx_fft ./ H;

%% ================= DEMODULATION =================

rx_symbols = rx_eq(:);

rx_bits = qpsk_demod(rx_symbols);

%% ================= BER =================

[min_len] = min(length(bits), length(rx_bits));
BER = sum(bits(1:min_len) ~= rx_bits(1:min_len)) / min_len;

fprintf("BER (Ideal Channel) = %f\n", BER);