function plot_ber(SNR_dB_vec, curves, title_str)
% PLOT_BER  Semi-log BER vs SNR plot.
%
%   plot_ber(SNR_dB_vec, curves, title_str)
%
%   Inputs:
%     SNR_dB_vec - SNR axis vector (dB).
%     curves     - Struct array, each with fields:
%                    .ber    - BER vector (same length as SNR_dB_vec).
%                    .label  - Legend label string.
%                    .style  - Plot style string (e.g. 'r-o').
%     title_str  - Figure title string.

figure;
hold on;

for k = 1:numel(curves)
    semilogy(SNR_dB_vec, curves(k).ber, curves(k).style, ...
             'LineWidth', 2, 'MarkerSize', 7);
end

grid on;
xlabel('SNR (dB)');
ylabel('BER');
title(title_str);
legend({curves.label}, 'Location', 'southwest');

hold off;

end
