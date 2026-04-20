import numpy as np
import scipy.io as sio


# ----------------------------
# QPSK Modulation
# ----------------------------
def qpsk_mod(bits):
    return (2*(bits % 2) - 1 + 1j*(2*(bits // 2) - 1)) / np.sqrt(2)


# ----------------------------
# Steering Vector (ULA)
# ----------------------------
def steering_vector(N, angle):
    n = np.arange(N)
    return np.exp(1j * np.pi * n * np.sin(angle))


# ----------------------------
# mmWave Channel
# ----------------------------
def generate_mmwave_channel(N_rx, N_tx, num_paths):
    H = np.zeros((N_rx, N_tx), dtype=complex)

    for _ in range(num_paths):
        alpha = (np.random.randn() + 1j*np.random.randn()) / np.sqrt(2)
        theta = np.random.uniform(-np.pi/2, np.pi/2)
        phi = np.random.uniform(-np.pi/2, np.pi/2)

        a_rx = steering_vector(N_rx, theta)
        a_tx = steering_vector(N_tx, phi)

        H += alpha * np.outer(a_rx, np.conj(a_tx))

    return H / np.sqrt(num_paths)


# ----------------------------
# RIS Dataset Generator
# ----------------------------
def generate_ris_dataset_paper():
    np.random.seed(42)

    # PARAMETERS
    N_sc = 64
    N_cp = 16
    M = 4
    num_frames = 40000
    SNR_dB = np.arange(0, 31, 5)

    num_BS = 4
    RIS_els = 16 * 16
    num_paths = 20

    # STORAGE
    DL_RIS_Y_pilot = np.zeros((num_frames, len(SNR_dB), N_sc), dtype=complex)
    DL_RIS_Y_data = np.zeros((num_frames, len(SNR_dB), N_sc), dtype=complex)
    DL_RIS_X_data = np.zeros((num_frames, N_sc), dtype=complex)

    SER_LS = np.zeros(len(SNR_dB))
    SER_MMSE = np.zeros(len(SNR_dB))

    # PILOT
    pilot_bits = np.random.randint(0, M, N_sc)
    X_p = qpsk_mod(pilot_bits)

    print("Generating RIS dataset (corrected)...")

    for frame in range(num_frames):

        # DATA
        data_bits = np.random.randint(0, M, N_sc)
        X_d = qpsk_mod(data_bits)
        DL_RIS_X_data[frame] = X_d

        # CHANNELS
        H_BR = generate_mmwave_channel(RIS_els, num_BS, num_paths)
        H_RU = generate_mmwave_channel(1, RIS_els, num_paths)

        theta = np.exp(1j * np.random.uniform(0, 2*np.pi, RIS_els))
        Theta = np.diag(theta)

        H_eff = H_RU @ Theta @ H_BR
        P = np.ones((num_BS, 1)) / np.sqrt(num_BS)

        H_eff_scalar = (H_eff @ P).flatten()[0]

        # ----------------------------
        # Frequency-selective channel
        # ----------------------------
        h_time = np.zeros(N_sc, dtype=complex)

        for p in range(num_paths):
            delay = np.random.randint(0, N_cp)
            gain = (np.random.randn() + 1j*np.random.randn()) / np.sqrt(2*num_paths)
            h_time[delay] += gain

        H_eff_k = np.fft.fft(h_time) * H_eff_scalar

        # ----------------------------
        # SNR LOOP
        # ----------------------------
        for snr_idx, snr in enumerate(SNR_dB):

            noise_var = 10**(-snr/10)
            signal_power = 1
            noise_var = signal_power * noise_var

            noise_p = (np.random.randn(N_sc) + 1j*np.random.randn(N_sc)) * np.sqrt(noise_var/2)
            noise_d = (np.random.randn(N_sc) + 1j*np.random.randn(N_sc)) * np.sqrt(noise_var/2)

            y_p = X_p * H_eff_k + noise_p
            y_d = X_d * H_eff_k + noise_d

            # SAVE DATA
            DL_RIS_Y_pilot[frame, snr_idx] = y_p
            DL_RIS_Y_data[frame, snr_idx] = y_d

            # ----------------------------
            # LS
            # ----------------------------
            H_LS = y_p / (X_p + 1e-8)
            X_hat_LS = y_d / (H_LS + 1e-8)

            # ----------------------------
            # MMSE (stable)
            # ----------------------------
            H_MMSE = (np.conj(X_p) * y_p) / (np.abs(X_p)**2 + noise_var)
            X_hat_MMSE = y_d / (H_MMSE + 1e-8)

            # ----------------------------
            # HARD DECISIONS
            # ----------------------------
            bits_LS_real = (np.real(X_hat_LS) > 0).astype(int)
            bits_LS_imag = (np.imag(X_hat_LS) > 0).astype(int)

            bits_MMSE_real = (np.real(X_hat_MMSE) > 0).astype(int)
            bits_MMSE_imag = (np.imag(X_hat_MMSE) > 0).astype(int)

            true_real = (np.real(X_d) > 0).astype(int)
            true_imag = (np.imag(X_d) > 0).astype(int)

            # ERROR COUNT
            SER_LS[snr_idx] += np.sum(bits_LS_real != true_real) + np.sum(bits_LS_imag != true_imag)
            SER_MMSE[snr_idx] += np.sum(bits_MMSE_real != true_real) + np.sum(bits_MMSE_imag != true_imag)

    # ----------------------------
    # NORMALIZE (IMPORTANT)
    # ----------------------------
    SER_LS /= (num_frames * N_sc * 2)
    SER_MMSE /= (num_frames * N_sc * 2)

    # ----------------------------
    # SAVE
    # ----------------------------
    sio.savemat(r'C:\Users\Nihit Reddy\OneDrive\Desktop\sem6\mlwc\assignment2\RIS_Dataset.mat', {
        'DL_RIS_Y_pilot': DL_RIS_Y_pilot,
        'DL_RIS_Y_data': DL_RIS_Y_data,
        'DL_RIS_X_data': DL_RIS_X_data,
        'SNR_dB': SNR_dB,
        'SER_LS': SER_LS,
        'SER_MMSE': SER_MMSE
    })

    print("RIS dataset generated successfully (final correct version).")


# ----------------------------
# RUN
# ----------------------------
if __name__ == "__main__":
    generate_ris_dataset_paper()