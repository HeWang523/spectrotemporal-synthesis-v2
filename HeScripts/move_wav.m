clc
clear

Dir_mm = {'/165_natural_sounds_model_matched', '/MarmosetVocal_model_matched', '/BirdSong_model_matched'};
save_dir = '/Agg_sounds';
cate = 3;
% directory containing the audio waveform
input_directory = [pwd Dir_mm{cate}];
fname_list_str  = dir( fullfile(input_directory, '*.wav'));
fname_list = {fname_list_str.name};
out_directory = [pwd save_dir];
s_n = length(fname_list);


for i = 1:s_n
    source = [input_directory '/' fname_list{i}];
    if cate == 2
        head = 'stim00';
        fname_split = regexp(fname_list{i}, '\_', 'split');
        new_fname = head;
        for j = 2:length(fname_split)
            new_fname = [new_fname '_' fname_split{j}];
        end
    elseif cate == 3
        head = 'stim01';
        new_fname = [head '_' fname_list{i}];
    elseif cate == 1
        new_fname = fname_list{i};
    end
    dest = [out_directory '/' new_fname];
    copyfile(source, dest, 'f');


end
