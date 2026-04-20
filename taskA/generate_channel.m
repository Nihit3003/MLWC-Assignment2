function h = generate_channel(N, use_winner2)
% GENERATE_CHANNEL  Return a multipath channel impulse response.
%
%   h = generate_channel(N)
%   h = generate_channel(N, use_winner2)
%
%   Tries to use the WINNER II model when use_winner2 is true (default).
%   Falls back to a 16-tap exponential PDP Rayleigh channel if WINNER II
%   is unavailable or throws an error.
%
%   The returned vector h always satisfies length(h) <= N.

if nargin < 2
    use_winner2 = true;
end

h = [];

if use_winner2
    try
        cfg  = winner2.wimparset;
        layout = winner2.layoutparset;
        linkpar = winner2.linkparset;

        cfg.Scenario        = 'A1';
        cfg.CenterFrequency = 2.5e9;

        layout.Stations(1).Pos      = [0; 0; 0];
        layout.Stations(2).Pos      = [10; 0; 0];
        layout.Stations(1).Velocity = [0; 0; 0];
        layout.Stations(2).Velocity = [0; 0; 0];

        [H, ~] = winner2.wim(cfg, layout, linkpar);
        h = squeeze(H(1,1,:,1));
        h = h(:);
    catch
        % WINNER II unavailable — fall through to fallback
    end
end

if isempty(h)
    % Fallback: 16-tap exponential power-delay profile
    L   = 16;
    pdp = exp(-(0:L-1) / 2);
    pdp = pdp / sum(pdp);
    h   = (randn(L,1) + 1j*randn(L,1)) .* sqrt(pdp(:) / 2);
end

if length(h) > N
    h = h(1:N);
end

end
