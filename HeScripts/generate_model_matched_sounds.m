clc
clear
% Demonstrates how to run the synthesis algorithm
% name of the audio waveform
load('C:\modelmatch\spectrotemporal-synthesis-v2-master\165_natural_sounds\category_labels.mat')
% fname = [stim_names{1} + ".wav"];
% fname = "000_TWM96.wav";

Dir = {'/165_natural_sounds_441', '/MarmosetVocal_filtered', '/BirdSong_selected'};
Dir_mm = {'/165_natural_sounds_model_matched', '/MarmosetVocal_model_matched', '/BirdSong_model_matched'};
% cate = 2;
for cate = 2:2
    % directory containing the audio waveform
    input_directory = [pwd Dir{cate}];
    % directory to save results of the synthesis process
    output_directory = [pwd Dir_mm{cate}];
    if ~exist(output_directory,'dir')
        mkdir(output_directory);
    end

    fname_list_str  = dir( fullfile(input_directory, '*.wav'));
    fname_list = {fname_list_str.name};

    s_n = length(fname_list);

    for i = 1:1
        if cate == 2 || cate == 3
            fname = fname_list{i};
        elseif cate == 1
            fname = [stim_names{i} + ".wav"];
        end
        % read parameters
        P = synthesis_parameters_default;
        if cate == 1
            P.audio_sr = 44100;
        elseif cate == 2 || cate == 3
            P.audio_sr = 50000;
        end
        P.n_filts = 160;
        P.lo_freq_hz = 50;
        P.n_iter = 100;
        P.overcomplete = 0;
        % run synthesis
        run_spectrotemporal_synthesis(P, fname, input_directory, output_directory);
    end
end