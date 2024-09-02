clc;
clear;

fs = 50000;
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






