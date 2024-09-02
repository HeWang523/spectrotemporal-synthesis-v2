function SpecdBOut = MarmosetExcitationPatternv2(SpecAmpIn, Filterbandk, FreqOut)
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


for i = 1:length(FreqOut)
    W = Filterbandk(:,i);
    SpecdBOut(i) = 10*log10( sum(SpecAmpIn.^2 .*W) );
end


% figure
% plot(FreqOut/1000, SpecdBOut);
% ERB_r_raw       =[	0       0.03735 0.00848 0.00193 0.00095];
% ERB_p_Out   = interp1(ERB_freq_Raw,	ERB_p_raw, FreqOut,'spline');
% ERB_r_Out   = interp1(ERB_freq_Raw,	ERB_r_raw, FreqOut,'spline');
end




