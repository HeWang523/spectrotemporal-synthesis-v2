clc
clear

%% load in data, set output dir

in_Dir = {'/MarmosetVocal_filtered', '/Bird_song_resample'};
out_Dir = {'/MarmosetVocal_2_s', '/BirdSong_2_s'};
marmoset_input_dir = [pwd in_Dir{1}];
bird_input_dir = [pwd in_Dir{2}];

marmoset_output_dir = [pwd out_Dir{1}];
bird_output_dir = [pwd out_Dir{2}];

fname_list_str_marmoset  = dir( fullfile(marmoset_input_dir, '*.wav'));
fname_list_marmoset = {fname_list_str_marmoset.name};

fname_list_str_bird  = dir( fullfile(bird_input_dir, '*.wav'));
fname_list_bird = {fname_list_str_bird.name};


marmoset_sound_n = length(fname_list_marmoset);
bird_sound_n = length(fname_list_bird);


%% Generate 2 s marmoset vocal sounds
fs = 50000;
target_dur = fs * 2;
interval_zeros = zeros(0.2*fs,1);

marmoset_2_s = zeros(target_dur,1);

for i = 1:marmoset_sound_n
    [wav_marmoset, ~] = audioread([marmoset_input_dir '/'  fname_list_marmoset{i}]);
    
    idx = find(abs(wav_marmoset) > 10e-4); % to remove the startup and end blank, find the indeces of wav that are greater than 0 (10e-4)

    wav_marmoset_no_blank = wav_marmoset(idx(1):idx(end));

    marmoset_length_no_blank(i) = (idx(end) - idx(1) + 1) / fs; % store the real length of the audio

    % Cat the no-blank wav 3 times to make sure the audio is greater than 2
    % s

    marmoset_no_blank_cat = [wav_marmoset_no_blank; interval_zeros; wav_marmoset_no_blank; interval_zeros; wav_marmoset_no_blank];

    marmoset_2_s = marmoset_no_blank_cat(1:target_dur);
    marmoset_2_s = 0.01*marmoset_2_s / rms(marmoset_2_s);
    fname = [marmoset_output_dir '/'  fname_list_marmoset{i}];
    audiowrite(fname,marmoset_2_s, fs);

    % Now find the audio length after rearraging
    idx = find(abs(marmoset_2_s) > 10e-4); 

    marmoset_length_no_blank_after_rearrange(i) = (idx(end) - idx(1) + 1) / fs; % store the real length of the audio

end

%% Generate 2 s Bird song 

bird_2_s = zeros(target_dur,1);


for i = 1:bird_sound_n
    [wav_bird, ~] = audioread([bird_input_dir '/'  fname_list_bird{i}]);
    
    idx = find(abs(wav_bird) > 10e-4); % to remove the startup and end blank, find the indeces of wav that are greater than 0 (10e-4)

    wav_bird_no_blank = wav_bird(idx(1):idx(end));

    bird_length_no_blank(i) = (idx(end) - idx(1) + 1) / fs; % store the real length of the audio

    % Cat the no-blank wav 3 times to make sure the audio is greater than 2
    % s

    bird_no_blank_cat = [wav_bird_no_blank; interval_zeros; wav_bird_no_blank; interval_zeros; wav_bird_no_blank; interval_zeros; wav_bird_no_blank];

    bird_2_s = bird_no_blank_cat(1:target_dur);
    bird_2_s = 0.01*bird_2_s / rms(bird_2_s);
    fname = [bird_output_dir '/'  fname_list_bird{i}];
    audiowrite(fname,bird_2_s, fs);
    
    % Now find the audio length after rearraging
    idx = find(abs(bird_2_s) > 10e-4); 

    bird_length_no_blank_after_rearrange(i) = (idx(end) - idx(1) + 1) / fs; % store the real length of the audio

end

%% Plot histograms to visualize the distribution of audio length before and after the rearranging
figure(1);
h = histogram(marmoset_length_no_blank, 10);
h.FaceColor = [0 0.5 0.5];
h.EdgeColor = 'r';
xline(2,'-','Desired Duration of Stimulus');
xlim([1.6,2.1]);
title("Histogram of Audio of Marmoset Vocal before Rearranging");
xlabel("Length of Audio (sec)");
ylabel("Count");


figure(2);
h = histogram(marmoset_length_no_blank_after_rearrange);
h.FaceColor = [0 0.5 0.5];
h.EdgeColor = 'r';
title("Histogram of Audio of Marmoset Vocal After Rearranging");
xlabel("Length of Audio (sec)");
ylabel("Count");

figure(3);
h = histogram(bird_length_no_blank, 50);
h.FaceColor = [0 0.5 0.5];
h.EdgeColor = 'r';
xline(2,'-','Desired Duration of Stimulus');
title("Histogram of Audio of Bird Songs Before Rearranging");
xlabel("Length of Audio (sec)");
ylabel("Count");

figure(4);
h = histogram(bird_length_no_blank_after_rearrange);
h.FaceColor = [0 0.5 0.5];
h.EdgeColor = 'r';
title("Histogram of Audio of Bird Songs After Rearranging");
xlabel("Length of Audio (sec)");
ylabel("Count");
