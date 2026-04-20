% task2_channel_estimation.m
%
% Task 2 — Channel estimation: LS vs MMSE.
% Runs two scenarios back-to-back:
%   (a) 64 pilots  — every subcarrier is a pilot (best-case estimation).
%   (b)  8 pilots  — sparse pilots, requires interpolation for LS.
%
% Run:  task2_channel_estimation

clc; clear; close all;

%% ---- Shared parameters -------------------------------------------------
N           = 64;
cp_len      = 16;
SNR_dB_vec  = -5:2:20;
numIter     = 1500;

base_cfg = struct( ...
    'N',          N,          ...
    'cp_len',     cp_len,     ...
    'SNR_dB_vec', SNR_dB_vec, ...
    'numIter',    numIter,    ...
    'dnn',        [],         ...   % no DNN in this task
    'use_winner2',false);           % use fast fallback channel

%% ---- (a) 64 pilots -----------------------------------------------------
fprintf('\n=== 64-pilot simulation ===\n');
cfg_64 = base_cfg;
cfg_64.pilot_idx = 1:N;

[BER_LS_64, BER_MMSE_64, ~] = run_sim(cfg_64);

%% ---- (b) 8 pilots ------------------------------------------------------
fprintf('\n=== 8-pilot simulation ===\n');
cfg_8 = base_cfg;
cfg_8.pilot_idx  = 1:8:N;
cfg_8.numIter    = 4000;           % more iterations for sparse pilots

[BER_LS_8, BER_MMSE_8, ~] = run_sim(cfg_8);

%% ---- Plots -------------------------------------------------------------
curves_64(1) = struct('ber', BER_LS_64,   'label', 'LS',   'style', 'r-o');
curves_64(2) = struct('ber', BER_MMSE_64, 'label', 'MMSE', 'style', 'b-s');
plot_ber(SNR_dB_vec, curves_64, 'BER vs SNR — 64 Pilots (LS vs MMSE)');

curves_8(1) = struct('ber', BER_LS_8,   'label', 'LS (8 pilots)',   'style', 'r-o');
curves_8(2) = struct('ber', BER_MMSE_8, 'label', 'MMSE (8 pilots)', 'style', 'b-s');
plot_ber(SNR_dB_vec, curves_8, 'BER vs SNR — 8 Pilots (LS vs MMSE)');
