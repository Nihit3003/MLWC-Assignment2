% task1_ofdm_ideal.m
%
% Task 1 — Basic OFDM system with ideal (perfect) channel knowledge.
% Demonstrates end-to-end OFDM: modulation, multipath channel, demodulation.
%
% Run:  task1_ofdm_ideal

clc; clear; close all;

%% ---- Parameters --------------------------------------------------------
N          = 64;        % subcarriers
cp_len     = 16;        % cyclic-prefix length
M          = 4;         % QPSK
numSymbols = 1000;      % OFDM symbols
SNR_dB     = 10;        % fixed SNR for this demo

%% ---- Generate & modulate bits ------------------------------------------
bits      = randi([0 1], numSymbols * N * log2(M), 1);
symbols   = qpsk_mod(bits);
sym_mat   = reshape(symbols, N, []);

%% ---- OFDM transmitter (column-by-column) --------------------------------
tx_signal = [];
for k = 1:size(sym_mat, 2)
    tx_signal = [tx_signal; ofdm_tx(sym_mat(:,k), N, cp_len)]; %#ok<AGROW>
end

%% ---- Multipath channel -------------------------------------------------
L  = 5;
h  = (randn(L,1) + 1j*randn(L,1)) / sqrt(2*L);
H  = fft(h, N);            % frequency-domain channel (for ideal equalization)

rx_signal = conv(tx_signal, h);
rx_signal = rx_signal(1:length(tx_signal));
rx_signal = awgn(rx_signal, SNR_dB, 'measured');

%% ---- OFDM receiver ------------------------------------------------------
frame_len = N + cp_len;
rx_mat    = reshape(rx_signal, frame_len, []);
Y_mat     = zeros(N, size(rx_mat, 2));
for k = 1:size(rx_mat, 2)
    Y_mat(:,k) = ofdm_rx(rx_mat(:,k), N, cp_len);
end

%% ---- Ideal equalization (perfect channel knowledge) --------------------
H_rep  = repmat(H, 1, size(Y_mat, 2));
rx_eq  = Y_mat ./ H_rep;

%% ---- Demodulate --------------------------------------------------------
rx_bits = qpsk_demod(rx_eq(:));

%% ---- BER ---------------------------------------------------------------
n   = min(numel(bits), numel(rx_bits));
BER = sum(bits(1:n) ~= rx_bits(1:n)) / n;
fprintf('BER (Ideal Channel, SNR=%d dB) = %.5f\n', SNR_dB, BER);
