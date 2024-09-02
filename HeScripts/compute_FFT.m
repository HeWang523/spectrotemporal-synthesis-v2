function [pxx,f] = compute_FFT(wav,fs)
% L = length(wav);
% raw_p = fft(wav);
% raw_p = abs(raw_p / L)'; % amp
% 
% f = fs/L*(0:L-1);
% pxx = raw_p(1:L/2);
% pxx = pxx;
% f = f(1:L/2);

[pxx,f] = pspectrum(wav, fs, 'FrequencyLimits', [50 fs/2],'FrequencyResolution',4);
pxx = movmean(pxx, 2048);

end

