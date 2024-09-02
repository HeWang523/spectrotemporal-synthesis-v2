function [SpecdBOut, W_out] = MarmosetExcitationPatternv1(FreqIn, SpecAmpIn, FreqOut)
% MarmosetExPatnSteady caclulate the steady excitation pattern of marmoset
% auditory periphery. (Assume the input sound is a steady sound)
%   FreqIn,     the discrete INPUT frequency bins,
%               should be in linear scale
%   SpecAmpIn,  the Amp(rms) for each INPUT frequency bin, 
%   FreqOut,    the discrete OUTPUT frequency bins,
%               can be in linear/log/ERB scale
%   SpecdBOut,  the Amp(dB relative) for each OUTPUT frequency bin
%       assume Amp=1(rms) is X dB SPL
%       the number in SpecdBOut is dB number relative to X in each bin

%   Technical NOTES:
%       (1) interpolation direct on p, r generally fails. Especially for r,
%           when not 0, the side leakage is huge. Even w/ 0 input, there is
%           still a high baseline "level"
%       (2) interpolation on ERB 1st, and then use roex filter that only
%           has p but not r is much better.

% SpecAmpIn = SpecAmpIn ./ std(SpecAmpIn);

ERB_freq_Raw    =[  250     500     1000    7000    16000];
ERB_raw         =[  90.97   126.85  180.51  460.83  2282.71];
% data for 500 1000 7000 16000 are from 4 monkeys
% data for 250 is from 1 monkey and r is assumed to be 0;

% % ERB_Out =   interp1(ERB_freq_Raw,	ERB_raw, FreqOut,'spline'); 
% 
ERB_Out_marmosetERB =   interp1(ERB_freq_Raw,	ERB_raw, FreqOut,'pchip'); 
% ERB_Out_1over6 = FreqOut/6;
% ERB_Out_1over12 = FreqOut/12;
ERB_Out_1over6 = FreqOut * (2 ^ (1/12) - 2 ^ (-1/12));
ERB_Out_1over12 = FreqOut * (2 ^ (1/12) - 1);
ERB_1990 = 24.7*(4.37*FreqOut/1000+1);

% % ERB_p_Out = 4.322 * FreqOut./ERB_Out; 
ERB_p_Out = 4 * FreqOut./ERB_Out_1over12; 

ERB_r_Out = 0*ERB_p_Out;

% figure
% loglog(FreqOut, ERB_Out_1over6);
% hold on
% loglog(FreqOut, ERB_Out_1over12);
% loglog(FreqOut, ERB_1990);
% scatter(ERB_freq_Raw, ERB_raw)
% grid on
% legend('1/6 octave', '1/12 octave', 'Glasberg & Moore, 1990', 'Marmoset ERB');
% xlabel('Central Frequency (Hz)');
% ylabel('Bandwidth (Hz)');
% hold off

for i = 1:length(FreqOut)
    % for each output Freqbin
    g = abs( (FreqIn - FreqOut(i)) / FreqOut(i) );
        % In 
    W = ERB_r_Out(i) + (1-ERB_r_Out(i)) * (1+ERB_p_Out(i)*g) .*...
        exp(-ERB_p_Out(i)*g);
        % In
    SpecdBOut(i) = 10*log10( sum(SpecAmpIn.^2 .*W) );
        % Out(i)  

    W_out(:,i) = W;
end


% figure
% plot(FreqOut/1000, SpecdBOut);
% ERB_r_raw       =[	0       0.03735 0.00848 0.00193 0.00095];
% ERB_p_Out   = interp1(ERB_freq_Raw,	ERB_p_raw, FreqOut,'spline');
% ERB_r_Out   = interp1(ERB_freq_Raw,	ERB_r_raw, FreqOut,'spline');
end


