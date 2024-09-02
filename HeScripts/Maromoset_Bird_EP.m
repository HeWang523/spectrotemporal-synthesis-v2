clc;
clear;
close all;
%% Compare excitation pattern

in_Dir = {'/MarmosetVocal_filtered', '/BirdSong_selected', '/165_natural_sounds_441'};
marmoset_input_dir = [pwd in_Dir{1}];
bird_input_dir = [pwd in_Dir{2}];
ns_input_dir = [pwd in_Dir{3}];

fname_list_str_marmoset  = dir( fullfile(marmoset_input_dir, '*.wav'));
fname_list_marmoset = {fname_list_str_marmoset.name};

fname_list_str_bird  = dir( fullfile(bird_input_dir, '*.wav'));
fname_list_bird = {fname_list_str_bird.name};

fname_list_str_ns  = dir( fullfile(ns_input_dir, '*.wav'));
fname_list_ns = {fname_list_str_ns.name};

marmoset_sound_n = length(fname_list_marmoset);
bird_sound_n = length(fname_list_bird);
ns_sound_n = length(fname_list_ns);

%%
fs = 50000;

f_out = 50*2.^(0:log2(fs / (50 * 2))/160:log2(fs / (50 * 2)));

for i = 1:marmoset_sound_n
    [wav_marmoset, ~] = audioread([marmoset_input_dir '/'  fname_list_marmoset{i}]);
    wav_marmoset = wav_marmoset ./ std(wav_marmoset);
    % [pxx_marmoset(:,i),f_marmoset] = pspectrum(wav_marmoset, fs, 'FrequencyLimits', [0 fs/2],'FrequencyResolution',1024);
    [pxx_marmoset(i,:),f_marmoset] = compute_FFT(wav_marmoset, fs);
    [ExP_marmoset(:,i), ~] = MarmosetExcitationPatternv1(f_marmoset', pxx_marmoset(i,:), f_out);
end

for i = 1:bird_sound_n
    [wav_bird, ~] = audioread([bird_input_dir '/'  fname_list_bird{i}]);

    if size(wav_bird, 2) == 2
        continue;
    end

    wav_bird = wav_bird ./ std(wav_bird);
    % [pxx_bird(i,:),f_bird] = pspectrum(wav_bird, fs, 'FrequencyLimits', [0 fs/2],'FrequencyResolution',1024);
    [pxx_bird(i,:),f_bird] = compute_FFT(wav_bird, fs);
    [ExP_bird(:,i), ~] = MarmosetExcitationPatternv1(f_bird', pxx_bird(i,:), f_out);
end

for i = 1:ns_sound_n
    [wav_ns, ns_sr] = audioread([ns_input_dir '/'  fname_list_ns{i}]);

    wav_ns = wav_ns ./ std(wav_ns);
    % [pxx_ns(:,i),f_ns] = pspectrum(wav_ns, ns_sr, 'FrequencyLimits', [0 ns_sr/2],'FrequencyResolution',1024);
    [pxx_ns(i,:),f_ns] = compute_FFT(wav_ns, fs);
    [EXP_ns(:,i), ~] = MarmosetExcitationPatternv1(f_ns', pxx_ns(i,:), f_out);
end


average_EXP_marmoset = mean(ExP_marmoset,2);
average_EXP_Bird = mean(ExP_bird,2);
average_EXP_ns = mean(EXP_ns,2);

SD1=std(ExP_marmoset,[],2);
% SD1 = SD1 / 15;
NSDR1=(average_EXP_marmoset+SD1);
PSDR1=(average_EXP_marmoset-SD1);

SD2=std(ExP_bird,[], 2);
% SD2 = SD2 / 15;
NSDR2=(average_EXP_Bird+SD2);
PSDR2=(average_EXP_Bird-SD2);

SD3=std(EXP_ns,[], 2);
% SD2 = SD2 / 15;
NSDR3=(average_EXP_ns+SD3);
PSDR3=(average_EXP_ns-SD3);


figure(1)

semilogx(f_out/1000, average_EXP_marmoset, 'r', LineWidth=2);
hold on
semilogx(f_out/1000, average_EXP_Bird, 'b', LineWidth=2);
% semilogx(f_out/1000, average_EXP_ns, 'y', LineWidth=2);
legend('Marmoset Vocal', 'Bird Songs');
patch([f_out/1000 fliplr(f_out/1000)],[PSDR1' fliplr(NSDR1')], 'r', 'FaceAlpha',0.2, 'EdgeColor','none', 'HandleVisibility', 'off');

patch([f_out/1000 fliplr(f_out/1000)],[PSDR2' fliplr(NSDR2')], 'b', 'FaceAlpha',0.2, 'EdgeColor','none', 'HandleVisibility', 'off');
% patch([f_out/1000 fliplr(f_out/1000)],[PSDR3' fliplr(NSDR3')], 'y', 'FaceAlpha',0.2, 'EdgeColor','none', 'HandleVisibility', 'off');
ylim([-96, 20]);
xlim([0, fs/2000]);
box off;
xlabel('Frequency (kHz)');
ylabel('Level (dB)');
set(gcf,'color','w')
grid on;

hold off
%% Power spectrum
pxx_marmoset = pow2db(pxx_marmoset)';
pxx_bird = pow2db(pxx_bird)';
pxx_ns = pow2db(pxx_ns)';

average_PW_marmoset = mean(pxx_marmoset,2);
average_PW_Bird = mean(pxx_bird,2);
average_PW_ns = mean(pxx_ns,2);

SD1=std(pxx_marmoset,[],2);
% SD1 = SD1 / 15;
NSDR1=(average_PW_marmoset+SD1);
PSDR1=(average_PW_marmoset-SD1);

SD2=std(pxx_bird,[], 2);
% SD2 = SD2 / 15;
NSDR2=(average_PW_Bird+SD2);
PSDR2=(average_PW_Bird-SD2);

SD3=std(pxx_ns,[], 2);
% SD2 = SD2 / 15;
NSDR3=(average_PW_ns+SD3);
PSDR3=(average_PW_ns-SD3);


figure(2)

plot(f_marmoset/1000, average_PW_marmoset, 'r', LineWidth=2);
hold on
plot(f_marmoset/1000, average_PW_Bird, 'b', LineWidth=2);

legend('Marmoset Vocal', 'Bird Songs');
% patch([f_marmoset'/1000 fliplr(f_marmoset'/1000)],[PSDR1' fliplr(NSDR1')], 'r', 'FaceAlpha',0.2, 'EdgeColor','none', 'HandleVisibility', 'off');
% 
% patch([f_marmoset'/1000 fliplr(f_marmoset'/1000)],[PSDR2' fliplr(NSDR2')], 'b', 'FaceAlpha',0.2, 'EdgeColor','none', 'HandleVisibility', 'off');
xlim([0, fs/2000]);
ylim([-96, 20]);
box off
xlabel('Frequency (kHz)');
ylabel('Amplitude (dB)');
set(gcf,'color','w')
grid on;

hold off