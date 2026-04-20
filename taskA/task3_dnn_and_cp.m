% task3_dnn_and_cp.m
%
% Task 3 — DNN channel estimation + cyclic-prefix study.
% Produces three figures using the WINNER II channel model (with fallback):
%
%   Fig 1 — 64 pilots:  LS vs MMSE vs DNN
%   Fig 2 —  8 pilots:  LS vs MMSE vs DNN
%   Fig 3 — With CP vs Without CP (64 pilots, all three estimators)
%
% DNN weights must be present in the same folder as this script:
%   W1.npy  W2.npy  W3.npy  b1.npy  b2.npy  b3.npy
%
% Run:  task3_dnn_and_cp

clc; clear; close all;

addpath(genpath('winner2'));         % add WINNER II if available

%% ---- Load DNN weights --------------------------------------------------
dnn = load_dnn();                   % reads *.npy from current directory

%% ---- Shared parameters -------------------------------------------------
N          = 64;
cp_len     = 16;
SNR_dB_vec = -5:2:20;
numIter    = 1000;

base_cfg = struct( ...
    'N',           N,          ...
    'cp_len',      cp_len,     ...
    'SNR_dB_vec',  SNR_dB_vec, ...
    'numIter',     numIter,    ...
    'dnn',         dnn,        ...
    'use_winner2', true);

%% ---- Fig 1: 64 pilots --------------------------------------------------
fprintf('\n=== Fig 1: 64 pilots ===\n');
cfg1 = base_cfg;
cfg1.pilot_idx = 1:N;

[BER_LS_64, BER_MMSE_64, BER_DNN_64] = run_sim(cfg1);

c1(1) = struct('ber', BER_LS_64,   'label', 'LS',   'style', 'r-o');
c1(2) = struct('ber', BER_MMSE_64, 'label', 'MMSE', 'style', 'b-s');
c1(3) = struct('ber', BER_DNN_64,  'label', 'DNN',  'style', 'g-d');
plot_ber(SNR_dB_vec, c1, 'Fig 1: BER vs SNR — 64 Pilots (WINNER II)');

%% ---- Fig 2: 8 pilots ---------------------------------------------------
fprintf('\n=== Fig 2: 8 pilots ===\n');
cfg2 = base_cfg;
cfg2.pilot_idx = 1:8:N;

[BER_LS_8, BER_MMSE_8, BER_DNN_8] = run_sim(cfg2);

c2(1) = struct('ber', BER_LS_8,   'label', 'LS',   'style', 'r-o');
c2(2) = struct('ber', BER_MMSE_8, 'label', 'MMSE', 'style', 'b-s');
c2(3) = struct('ber', BER_DNN_8,  'label', 'DNN',  'style', 'g-d');
plot_ber(SNR_dB_vec, c2, 'Fig 2: BER vs SNR — 8 Pilots (WINNER II)');

%% ---- Fig 3: With CP vs Without CP (64 pilots) --------------------------
fprintf('\n=== Fig 3: With CP ===\n');
cfg3_cp = base_cfg;
cfg3_cp.pilot_idx = 1:N;
cfg3_cp.cp_len    = cp_len;

[BER_LS_CP, BER_MMSE_CP, BER_DNN_CP] = run_sim(cfg3_cp);

fprintf('\n=== Fig 3: Without CP ===\n');
cfg3_nocp = base_cfg;
cfg3_nocp.pilot_idx = 1:N;
cfg3_nocp.cp_len    = 0;

[BER_LS_noCP, BER_MMSE_noCP, BER_DNN_noCP] = run_sim(cfg3_nocp);

c3(1) = struct('ber', BER_LS_CP,     'label', 'LS (CP)',     'style', 'r-o');
c3(2) = struct('ber', BER_MMSE_CP,   'label', 'MMSE (CP)',   'style', 'b-s');
c3(3) = struct('ber', BER_DNN_CP,    'label', 'DNN (CP)',    'style', 'g-d');
c3(4) = struct('ber', BER_LS_noCP,   'label', 'LS (no CP)',  'style', 'r--^');
c3(5) = struct('ber', BER_MMSE_noCP, 'label', 'MMSE (no CP)','style', 'b--v');
c3(6) = struct('ber', BER_DNN_noCP,  'label', 'DNN (no CP)', 'style', 'g--x');
plot_ber(SNR_dB_vec, c3, 'Fig 3: With CP vs Without CP (WINNER II)');
