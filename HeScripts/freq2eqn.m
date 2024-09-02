% FUNCTION N_FILTER = FREQ2EQN(FREQ_HZ)
%
% Converts Hz to ERBs, each bandwidth is exactly 1/6 octave

function n_erb = freq2eqn(freq_Hz)

% n_erb = 6 * log(freq_Hz) - 21.6413;
n_erb = (2 / (2 ^ (1/12) - 2 ^ (-1/12))) * log(freq_Hz) - 65.8578;


