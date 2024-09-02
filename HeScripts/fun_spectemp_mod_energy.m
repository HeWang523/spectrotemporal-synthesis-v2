function spectemp_matrix = fun_spectemp_mod_energy(path,P,pltfig)
%% Generate Cochleagram
[wav_orig, ~] = audioread(path);
duration_sec = length(wav_orig) / P.audio_sr;
[audio_filts, audio_low_cutoff] = ...
    make_erb_cos_filts_quadruple2(duration_sec*P.audio_sr, P.audio_sr, ...
    P.n_filts, P.lo_freq_hz, P.audio_sr/2);

% remove filters below and above desired cutoffs
xi = audio_low_cutoff > P.lo_freq_hz - 1e-3 ...
    & audio_low_cutoff < P.audio_sr/2 + 1e-3;
audio_filts = audio_filts(:,xi);
audio_low_cutoff = audio_low_cutoff(xi);

[coch_orig, P.f, P.t, ~] = ...
    wav2coch(wav_orig, audio_filts, audio_low_cutoff, ...
    P.audio_sr, P.env_sr, P.compression_factor, P.logf_spacing);

% plot cochleograms
if pltfig
    figure;
    plot_cochleogram(coch_orig,P.f,P.t);
end

%% Spectrotemporal modulation

% Pad coch
% dimensions
% time x frequency
pad_value = mean(coch_orig(:));
n_temp_pad_smps = round(P.env_sr * 2);
% matrix set to the mean of the cochleogram
T_pad = pad_value * ones(n_temp_pad_smps, size(coch_orig,2));
temp_and_freq_padded_coch = [T_pad; coch_orig];
coch_orig = temp_and_freq_padded_coch;

[T,F] = size(coch_orig);

% 2D fourier transforms
FT_coch_orig = fft2(coch_orig);


temp_n = length(P.temp_mod_rates) * 2;
spec_n = length(P.spec_mod_rates);

for i = 1:temp_n
    for j = 1:spec_n
        if i > length(P.temp_mod_rates)
            idx = i - length(P.temp_mod_rates);
            temp_rate = -1 * P.temp_mod_rates(idx);
        else
            idx = i;
            temp_rate = P.temp_mod_rates(idx);
        end

        spec_rate = P.spec_mod_rates(j);
        Hts = filt_spectemp_mod(...
            spec_rate, temp_rate, F, T, P);
        filter_response = abs(ifft2(FT_coch_orig .* Hts));
        for k = 1:T
            energy(k) = sum(filter_response(k,:).^2);
        end
        avg_fr(j,i) = mean(energy);
    end
end

%% Average magnitude of positive and negative temporal modulations rates

avg_fr_pn = (avg_fr(:,1:temp_n/2) + avg_fr(:,(temp_n/2)+1:end))/2;
spectemp_matrix = flipud(avg_fr_pn);

if pltfig
    figure
    imagesc(flipud(avg_fr_pn));
end
end

