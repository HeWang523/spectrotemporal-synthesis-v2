clc
clear

in_Dir = {'/BirdSong_2_s'};
bird_input_dir = [pwd in_Dir{1}];
bird_output_dir = [pwd in_Dir{1}];

fname_list_str_bird  = dir( fullfile(bird_input_dir, '*.wav'));
fname_list_bird = {fname_list_str_bird.name};

bird_sound_n = length(fname_list_bird);

fs = 50000;

f_out = 50*2.^(0:log2(fs / (50 * 2))/160:log2(fs / (50 * 2)));

%% peak_filter
G_p = -100;
fc_p = 20750;
fb_p = 8500;
Q_p = fb_p / fc_p;

[b_p, a_p] = peak(G_p, fc_p, Q_p, fs);

%% notch filter

G_n = -10;
fc_n = 5000;
fb_n = 4000;
Q_n = fb_n / fc_n;

[b_n, a_n] = peak(G_n, fc_n, Q_n, fs);

for i = 1:bird_sound_n
    [wav_bird, ~] = audioread([bird_input_dir '/'  fname_list_bird{i}]);
    wav_bird = wav_bird ./ std(wav_bird);
    
    wav_bird = filtfilt(b_n, a_n, wav_bird);
    % wav_bird = filtfilt(b_p, a_p, wav_bird);
    

    wav_bird = 0.01*wav_bird / rms(wav_bird);
    fname = [bird_output_dir '/'  fname_list_bird{i}];
    audiowrite(fname,wav_bird, fs);

end