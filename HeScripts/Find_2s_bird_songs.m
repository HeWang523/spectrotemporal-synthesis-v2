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


%% Find 1.6 s - 2 s Bird song, which shares a similar distribution of duration with the length of marmoset vocal sounds
fs = 50000;
target_dur = fs * 2;

bird_2_s = zeros(target_dur,1);
selected_idx = 1;


for i = 1:bird_sound_n
    [wav_bird, ~] = audioread([bird_input_dir '/'  fname_list_bird{i}]);
    
    idx = find(abs(wav_bird) > 10e-4); % to remove the startup and end blank, find the indeces of wav that are greater than 0 (10e-4)

    wav_bird_no_blank = wav_bird(idx(1):idx(end));

    bird_length_no_blank(i) = (idx(end) - idx(1) + 1) / fs; % store the real length of the audio

    % if the length is greater than 1.6 s and shorter than 2 s

    if bird_length_no_blank(i) >= 1.6 && bird_length_no_blank(i) <= 2
        
        bird_2_s = wav_bird_no_blank;
        bird_2_s = 0.01*bird_2_s / rms(bird_2_s);
        fname = [bird_output_dir '/'  fname_list_bird{i}];
        audiowrite(fname,bird_2_s, fs);
        bird_length_no_blank_after_rearrange(selected_idx) = bird_length_no_blank(i);
        selected_idx = selected_idx + 1;
    end


end

%% Plot histograms to visualize the distribution of audio length before and after the rearranging

figure(1);
h = histogram(bird_length_no_blank, 50);
h.FaceColor = [0 0.5 0.5];
h.EdgeColor = 'r';
xline(2,'-','Desired Duration of Stimulus');
title("Histogram of Audio of Bird Songs Before Rearranging");
xlabel("Length of Audio (sec)");
ylabel("Count");

figure(2);
h = histogram(bird_length_no_blank_after_rearrange,20);
h.FaceColor = [0 0.5 0.5];
h.EdgeColor = 'r';
title("Histogram of Audio of Bird Songs After Selecting the Length");
xlabel("Length of Audio (sec)");
ylabel("Count");

%% Generate 2 s bird song 

in_Dir =  '/BirdSong_2_s';
input_directory = [pwd in_Dir];

out_Dir =  '/BirdSong_2_s';
out_directory = [pwd out_Dir];
fname_list_str  = dir( fullfile(input_directory, '*.wav'));
fname_list = {fname_list_str.name};

wav_des_length = fs * 2;
wav_out = zeros(wav_des_length,1);
s_n = length(fname_list);

for i = 1:s_n
    [wav_orig, wav_sr] = audioread([input_directory '/' fname_list{i}]);

    wav_real_length = length(wav_orig);

    if wav_real_length == wav_des_length
        wav_out = wav_orig;
    elseif wav_real_length < wav_des_length
        if wav_real_length <= wav_des_length/2
            wav_orig = [wav_orig wav_orig];
            wav_real_length = wav_real_length * 2;
        end
        res_length = wav_des_length - wav_real_length;
        if mod(res_length,2) == 0
            head = res_length / 2;
            tail = head;
        else
            head = (res_length + 1) / 2;
            tail = head - 1;
        end
        wav_out(head + 1: end - tail) = wav_orig;
    elseif wav_real_length > wav_des_length
        res_length = abs(wav_des_length - wav_real_length);
        if mod(res_length,2) == 0
            head = res_length / 2;
            tail = head;
        else
            head = (res_length + 1) / 2;
            tail = head - 1;
        end
        wav_out = wav_orig(head + 1: end - tail);
    end
    audiowrite([out_directory '/' fname_list{i}],wav_out, fs);
    
end
