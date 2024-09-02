% FUNCTION FREQ_HZ = EQN2FREQ(N_ERB)
%
% Converts ERBs to Hz, using the inverse formula of freq2eqn


function freq_Hz = eqn2freq(n_erb)

% freq_Hz = exp((n_erb + 21.6413) / 6);
freq_Hz = exp(2 * (n_erb + 65.8578) / ((2 ^ (1/12) - 2 ^ (-1/12))));


