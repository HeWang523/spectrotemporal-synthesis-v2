%
clc;
clear;
%
in_Dir = {'/MarmosetVocal', '/BirdSong'};
out_Dir = {'/MarmosetVocal_filtered', '/BirdSong_selected'};

marmoset_input_dir = [pwd in_Dir{1}];


marmoset_output_dir = [pwd out_Dir{1}];


fname_list_str_marmoset  = dir( fullfile(marmoset_input_dir, '*.wav'));
fname_list_marmoset = {fname_list_str_marmoset.name};




marmoset_sound_n = length(fname_list_marmoset);


% Resample marmoset vocal sounds back to 50 kHz

% Design a lowpass filter with cutoff frequency @ 45 kHz
fc = 45000;
fs = 100000;
fs_target = 50000;
[b,a] = cheby1(6,10,fc/(fs/2));

for i = 1:marmoset_sound_n
    [wav_orig, ~] = audioread([marmoset_input_dir '/'  fname_list_marmoset{i}]);

    wav_lowpass = filtfilt(b, a, wav_orig);

    % resample sound

    wav_resample = resample(wav_lowpass,fs_target, fs);
    wav_resample = 0.01*wav_resample / rms(wav_resample);

    % wav_resample = wav_resample ./ std(wav_resample);

    fname = [marmoset_output_dir '/'  fname_list_marmoset{i}];
    audiowrite(fname,wav_resample, fs_target);

    % % plot power spectrum and excitation pattern
    % 
    % [pxx_ori,f] = pspectrum(wav_resample, fs_target, 'FrequencyLimits', [0 fs_target/2],'FrequencyResolution',1024);
    % 
    % 
    % f_out = 50*2.^(0:log2(fs_target / (50 * 2))/160:log2(fs_target / (50 * 2)));
    % [ExP_ori, W_ori] = MarmosetExcitationPatternv1(f, pxx_ori, f_out);
    % 
    % figure
    % plot(f/1000,pow2db(pxx_ori), 'r', LineWidth=2);
    % ylim([-95, 10]);
    % xlim([0, fs_target / 2000]);
    % box off;
    % xlabel('Frequency (kHz)');
    % ylabel('Power Spectrum (dB)');
    % title('Power Spectrum');
    % set(gcf,'color','w')
    % 
    % 
    % figure
    % semilogx(f_out/1000, ExP_ori, 'r', LineWidth=2);
    % xlabel('Central Frequency (kHz)');
    % ylabel('Excitation Pattern (dB)');
    % title('Excitation Pattern');
    % ylim([-95, 5]);
    % xlim([50/1000, fs_target / 2000]);

end









