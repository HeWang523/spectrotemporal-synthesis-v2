%
clc;
clear;
%
in_Dir = {'/MarmosetVocal', '/BirdSong'};
out_Dir = {'/MarmosetVocal_filtered', '/Bird_song_resample'};

bird_input_dir = [pwd in_Dir{2}];


bird_output_dir = [pwd out_Dir{2}];


fname_list_str_bird  = dir( fullfile(bird_input_dir, '*.wav'));
fname_list_bird = {fname_list_str_bird.name};




bird_sound_n = length(fname_list_bird);


% Resample marmoset vocal sounds back to 50 kHz

% high pass bird song @ 3 kHz, design filter


fc = 3000;

fs_target = 50000;


    for i = 1:bird_sound_n
        [wav_orig, fs] = audioread([bird_input_dir '/'  fname_list_bird{i}]);
        [b,a] = butter(4, fc / (fs/ 2), 'high');
        

        % resample sound

        wav_resample = resample(wav_orig,fs_target, fs);
        wav_highpass = filtfilt(b, a, wav_resample);
        wav_highpass = 0.01*wav_highpass / rms(wav_highpass);


        fname = [bird_output_dir '/'  fname_list_bird{i}];
        audiowrite(fname,wav_highpass, fs_target);

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









