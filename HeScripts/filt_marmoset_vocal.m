clc
clear

Dir_out ='/MarmosetVocal_filtered';
Dir_in = '/MarmosetVocal';
cutoff = 1000;
output_directory = [pwd Dir_out];
input_directory = [pwd Dir_in];
fs = 100000;
fname_list_str  = dir( fullfile(input_directory, '*.wav'));
fname_list = {fname_list_str.name};

[b,a] = butter(2,cutoff * 2/fs, 'high');

s_n = length(fname_list);
wav_orig = [];
for i = 1:s_n
    synth_mat_file = [output_directory fname_list{i}];
    [wav, wav_sr] = audioread([input_directory '/' fname_list{i}]);
    wav_temp= filtfilt(b, a, wav);
    
    [pxx_ori(i,:),f] = pspectrum(wav_temp, wav_sr, 'FrequencyLimits', [0 25000],'FrequencyResolution',512);
    figure(1)
    plot(f/1000,pow2db(pxx_ori(i,:)), 'r');
    audiowrite([output_directory '/' fname_list{i}],wav_temp, wav_sr);
    
    [wav, wav_sr] = audioread([output_directory '/' fname_list{i}]);
    % wav= filtfilt(b, a, wav);
    [pxx_ori(i,:),f] = pspectrum(wav, wav_sr, 'FrequencyLimits', [0 25000],'FrequencyResolution',512);
    figure(2)
    
    plot(f/1000,pow2db(pxx_ori(i,:)), 'r');

end
average_p_ori = mean(pxx_ori,1);
figure(3)

plot(f/1000,pow2db(average_p_ori), 'r');
