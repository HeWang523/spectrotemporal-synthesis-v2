clc;
clear;
close all;
fs = 100000;
%% Plot pxx of marmoset vocal, bird song, 165 natural sounds in the same figure
f = load("Agg_sounds\f.mat").f;

pxx_165 = load("165_natural_sounds_Info\pxx.mat");
pxx_mar = load("MarmosetVocal_info\pxx.mat");
pxx_bird = load("BirdSong_Info\pxx.mat");

pxx_165_mm = load("165_natural_sounds_Info\pxx_mm.mat");
pxx_mar_mm = load("MarmosetVocal_info\pxx_mm.mat");
pxx_bird_mm = load("BirdSong_Info\pxx_mm.mat");

average_p_165 = mean(pow2db(pxx_165.pxx_ori),1);
average_p_mar = mean(pow2db(pxx_mar.pxx_ori),1);
average_p_bird = mean(pow2db(pxx_bird.pxx_ori),1);
% f = load("f_in.mat").f;
SD1=std(pow2db(pxx_165.pxx_ori),1);
NSDR1=(average_p_165+SD1);
PSDR1=(average_p_165-SD1);

SD2=std(pow2db(pxx_mar.pxx_ori),1);
NSDR2=(average_p_mar+SD2);
PSDR2=(average_p_mar-SD2);

SD3=std(pow2db(pxx_bird.pxx_ori),1);
NSDR3=(average_p_bird+SD3);
PSDR3=(average_p_bird-SD3);

figure

semilogx(f/1000,average_p_165, 'g','LineWidth',3);
hold on
semilogx(f/1000,average_p_mar,'b','LineWidth',3);
semilogx(f/1000,average_p_bird,'r','LineWidth',3);

patch([f'/1000 fliplr(f'/1000)],[NSDR1 fliplr(PSDR1)], 'g', 'FaceAlpha',0.2, 'EdgeColor','none', 'HandleVisibility', 'off');
patch([f'/1000 fliplr(f'/1000)],[NSDR2 fliplr(PSDR2)], 'b', 'FaceAlpha',0.2, 'EdgeColor','none', 'HandleVisibility', 'off');
patch([f'/1000 fliplr(f'/1000)],[NSDR3 fliplr(PSDR3)], 'r', 'FaceAlpha',0.2, 'EdgeColor','none', 'HandleVisibility', 'off');

ylim([-96, 0]);
xlim([0, fs/2000]);
legend('165 natural sounds', 'Marmoset Vocal', 'Bird Song', 'Location', 'southwest');
box off;
xlabel('Frequency (kHz)');
ylabel('Amplitude (dB)');
title('Power Spectrum of 165 Natural Sounds, Marmoset Vocal and Bird Songs');
set(gcf,'color','w')
grid on;
hold off
pxx_165 = pxx_165.pxx_ori;
pxx_mar = pxx_mar.pxx_ori;
pxx_bird = pxx_bird.pxx_ori;
pxx_165_mm = pxx_165_mm.pxx_mm;
pxx_mar_mm = pxx_mar_mm.pxx_mm;
pxx_bird_mm = pxx_bird_mm.pxx_mm;
%% Excitation pattern by FFT
in_Dir = {'/165_natural_sounds_Info/', '/MarmosetVocal_Info/', '/BirdSong_Info/'};
marmoset_input_dir = [pwd in_Dir{2}];
bird_input_dir = [pwd in_Dir{3}];
ns_input_dir = [pwd in_Dir{1}];

fname_list_str_marmoset  = dir( fullfile(marmoset_input_dir, '*orig.wav'));
fname_list_marmoset = {fname_list_str_marmoset.name};

fname_list_str_bird  = dir( fullfile(bird_input_dir, '*orig.wav'));
fname_list_bird = {fname_list_str_bird.name};

fname_list_str_ns  = dir( fullfile(ns_input_dir, '*orig.wav'));
fname_list_ns = {fname_list_str_ns.name};

marmoset_sound_n = length(fname_list_marmoset);
bird_sound_n = length(fname_list_bird);
ns_sound_n = length(fname_list_ns);

f_out = 50*2.^(0:log2(fs / (50 * 2))/160:log2(fs / (50 * 2)));


for i = 1:165
    [wav_ns, ~] = audioread([ns_input_dir '/'  fname_list_ns{i}]);
    wav_ns = wav_ns ./ std(wav_ns);
    % raw_p = fft(wav_ns);
    % raw_p = abs(raw_p)'; % amp
    % L = length(wav_ns);
    % f = fs/L*(0:L-1);
    % pxx_ns = raw_p(1:L/2);
    % pxx_ns = pxx_ns ./ max(pxx_ns);
    % f = f(1:L/2);
    % pxx_temp(i,:) = pxx_ns;
    [pxx_ns,f] = compute_FFT(wav_ns, fs);
    SpecdBOut_165(i,:) = MarmosetExcitationPatternv1(f, pxx_ns, f_out);
end

for i = 1:15
    [wav_marmoset, ~] = audioread([marmoset_input_dir '/'  fname_list_marmoset{i}]);
    wav_marmoset = wav_marmoset ./ std(wav_marmoset);
    % raw_p = fft(wav_marmoset);
    % raw_p = abs(raw_p)'; % amp
    % L = length(wav_marmoset);
    % f = fs/L*(0:L-1);
    % pxx_marmoset = raw_p(1:L/2);
    % pxx_marmoset = pxx_marmoset ./ max(pxx_marmoset);
    % f = f(1:L/2);
    [pxx_marmoset,f] = compute_FFT(wav_marmoset, fs);
    SpecdBOut_marvocal(i,:) = MarmosetExcitationPatternv1(f, pxx_marmoset, f_out);
end

for i = 1:15
    [wav_birdsong, ~] = audioread([bird_input_dir '/'  fname_list_bird{i}]);
    wav_birdsong = wav_birdsong ./ std(wav_birdsong);
    % raw_p = fft(wav_birdsong);
    % raw_p = abs(raw_p)'; % amp
    % L = length(wav_birdsong);
    % f = fs/L*(0:L-1);
    % pxx_bird = raw_p(1:L/2);
    % pxx_bird = pxx_bird ./ max(pxx_bird);
    % f = f(1:L/2);
    [pxx_bird,f] = compute_FFT(wav_birdsong, fs);
    SpecdBOut_bird(i,:) = MarmosetExcitationPatternv1(f, pxx_bird, f_out);
end

SD1=std(SpecdBOut_165,1);
NSDR1=(mean(SpecdBOut_165)+SD1);
PSDR1=(mean(SpecdBOut_165)-SD1);

SD2=std(SpecdBOut_marvocal,1);
NSDR2=(mean(SpecdBOut_marvocal)+SD2);
PSDR2=(mean(SpecdBOut_marvocal)-SD2);

SD3=std(SpecdBOut_bird,1);
NSDR3=(mean(SpecdBOut_bird)+SD3);
PSDR3=(mean(SpecdBOut_bird)-SD3);


figure
semilogx(f_out/1000, mean(SpecdBOut_165,1), 'g', LineWidth=2);
hold on
semilogx(f_out/1000, mean(SpecdBOut_marvocal,1), 'b',LineWidth=2);
semilogx(f_out/1000, mean(SpecdBOut_bird,1),'r', LineWidth=2);
patch([f_out/1000 fliplr(f_out/1000)],[NSDR1 fliplr(PSDR1)], 'g', 'FaceAlpha',0.4, 'EdgeColor','none', 'HandleVisibility', 'off');
patch([f_out/1000 fliplr(f_out/1000)],[NSDR2 fliplr(PSDR2)], 'b', 'FaceAlpha',0.4, 'EdgeColor','none', 'HandleVisibility', 'off');
patch([f_out/1000 fliplr(f_out/1000)],[NSDR3 fliplr(PSDR3)], 'r', 'FaceAlpha',0.4, 'EdgeColor','none', 'HandleVisibility', 'off');
legend('165 Natural Sounds', 'Marmoset Vocal', 'Bird Sounds', 'Location', 'southwest');
xlabel('Frequency (kHz)');
ylabel('Level (dB)');
title('1/12 octave weighted - 165 Natural Sounds, Marmoset Vocal and Bird Songs');
ylim([-96, 15]);
xlim([0, fs/2000]);
grid on;
set(gcf,'color','w')
box off;
hold off

%% 
in_Dir = {'/165_natural_sounds_Info/', '/MarmosetVocal_Info/', '/BirdSong_Info/'};
marmoset_input_dir = [pwd in_Dir{2}];
bird_input_dir = [pwd in_Dir{3}];
ns_input_dir = [pwd in_Dir{1}];

fname_list_str_marmoset  = dir( fullfile(marmoset_input_dir, '*synth.wav'));
fname_list_marmoset = {fname_list_str_marmoset.name};

fname_list_str_bird  = dir( fullfile(bird_input_dir, '*synth.wav'));
fname_list_bird = {fname_list_str_bird.name};

fname_list_str_ns  = dir( fullfile(ns_input_dir, '*synth.wav'));
fname_list_ns = {fname_list_str_ns.name};

marmoset_sound_n = length(fname_list_marmoset);
bird_sound_n = length(fname_list_bird);
ns_sound_n = length(fname_list_ns);

for i = 1:165
    [wav_ns, ~] = audioread([ns_input_dir '/'  fname_list_ns{i}]);
    wav_ns = wav_ns ./ std(wav_ns);
    % raw_p = fft(wav_ns);
    % raw_p = abs(raw_p)'; % amp
    % L = length(wav_ns);
    % f = fs/L*(0:L-1);
    % pxx_ns = raw_p(1:L/2);
    % pxx_ns = pxx_ns ./ max(pxx_ns);
    % f = f(1:L/2);
    [pxx_ns,f] = compute_FFT(wav_ns, fs);
    SpecdBOut_165_mm(i,:) = MarmosetExcitationPatternv1(f, pxx_ns, f_out);
end

for i = 1:15
    [wav_marmoset, ~] = audioread([marmoset_input_dir '/'  fname_list_marmoset{i}]);
    wav_marmoset = wav_marmoset ./ std(wav_marmoset);
    % raw_p = fft(wav_marmoset);
    % raw_p = abs(raw_p)'; % amp
    % L = length(wav_marmoset);
    % f = fs/L*(0:L-1);
    % pxx_marmoset = raw_p(1:L/2);
    % pxx_marmoset = pxx_marmoset ./ max(pxx_marmoset);
    % f = f(1:L/2);
    [pxx_marmoset,f] = compute_FFT(wav_marmoset, fs);
    SpecdBOut_marvocal_mm(i,:) = MarmosetExcitationPatternv1(f, pxx_marmoset, f_out);
end

for i = 1:15
    [wav_birdsong, ~] = audioread([bird_input_dir '/'  fname_list_bird{i}]);
    wav_birdsong = wav_birdsong ./ std(wav_birdsong);
    % raw_p = fft(wav_birdsong);
    % raw_p = abs(raw_p)'; % amp
    % L = length(wav_birdsong);
    % f = fs/L*(0:L-1);
    % pxx_bird = raw_p(1:L/2);
    % pxx_bird = pxx_bird ./ max(pxx_bird);
    % f = f(1:L/2);
    [pxx_bird,f] = compute_FFT(wav_birdsong, fs);
    SpecdBOut_bird_mm(i,:) = MarmosetExcitationPatternv1(f, pxx_bird, f_out);
end

SD4=std(SpecdBOut_165_mm,1);
NSDR4=(mean(SpecdBOut_165_mm)+SD4);
PSDR4=(mean(SpecdBOut_165_mm)-SD4);

SD5=std(SpecdBOut_marvocal_mm,1);
NSDR5=(mean(SpecdBOut_marvocal_mm)+SD5);
PSDR5=(mean(SpecdBOut_marvocal_mm)-SD5);

SD6=std(SpecdBOut_bird_mm,1);
NSDR6=(mean(SpecdBOut_bird_mm)+SD6);
PSDR6=(mean(SpecdBOut_bird_mm)-SD6);

figure
semilogx(f_out/1000, mean(SpecdBOut_165,1), 'r', LineWidth=2);
hold on
semilogx(f_out/1000, mean(SpecdBOut_165_mm,1), 'b',LineWidth=2);
patch([f_out/1000 fliplr(f_out/1000)],[NSDR1 fliplr(PSDR1)], 'r', 'FaceAlpha',0.2, 'EdgeColor','none', 'HandleVisibility', 'off');
patch([f_out/1000 fliplr(f_out/1000)],[NSDR4 fliplr(PSDR4)], 'b', 'FaceAlpha',0.2, 'EdgeColor','none', 'HandleVisibility', 'off');
legend('165 Natural Sounds', 'Model Matched Sounds', 'Location', 'northwest');
xlabel('Frequency (kHz)');
ylabel('Level (dB)');
title('1/12 octave weighted - 165 Natural Sounds');
ylim([-96, 15]);
xlim([0, fs/2000]);
grid on;
hold off

figure
semilogx(f_out/1000, mean(SpecdBOut_marvocal,1), 'r', LineWidth=2);
hold on
semilogx(f_out/1000, mean(SpecdBOut_marvocal_mm,1), 'b',LineWidth=2);
patch([f_out/1000 fliplr(f_out/1000)],[NSDR2 fliplr(PSDR2)], 'r', 'FaceAlpha',0.2, 'EdgeColor','none', 'HandleVisibility', 'off');
patch([f_out/1000 fliplr(f_out/1000)],[NSDR5 fliplr(PSDR5)], 'b', 'FaceAlpha',0.2, 'EdgeColor','none', 'HandleVisibility', 'off');
legend('Marmoset Vocal', 'Model Matched Sounds', 'Location', 'northwest');
xlabel('Frequency (kHz)');
ylabel('Level (dB)');
title('1/12 octave weighted - Marmoset Vocal');
ylim([-96, 15]);
xlim([0, fs/2000]);
grid on;
hold off

figure
semilogx(f_out/1000, mean(SpecdBOut_bird,1), 'r', LineWidth=2);
hold on
semilogx(f_out/1000, mean(SpecdBOut_bird_mm,1), 'b',LineWidth=2);
patch([f_out/1000 fliplr(f_out/1000)],[NSDR3 fliplr(PSDR3)], 'r', 'FaceAlpha',0.2, 'EdgeColor','none', 'HandleVisibility', 'off');
patch([f_out/1000 fliplr(f_out/1000)],[NSDR6 fliplr(PSDR6)], 'b', 'FaceAlpha',0.2, 'EdgeColor','none', 'HandleVisibility', 'off');
legend('Bird Songs', 'Model Matched Sounds', 'Location', 'northwest');
xlabel('Frequency (kHz)');
ylabel('Level (dB)');
title('1/12 octave weighted - Bird Songs');
ylim([-96, 15]);
xlim([0, fs/2000]);
grid on;
hold off
