clc
clear
% Demonstrates how to run the synthesis algorithm
% name of the audio waveform
load('C:\modelmatch\spectrotemporal-synthesis-v2-master\165_natural_sounds\category_labels.mat')
% fname = [stim_names{1} + ".wav"];
fname = "000_TWM96.wav";

Dir = {'/165_natural_sounds', '/MarmosetVocal', '/BirdSong'};
Dir_mm = {'/165_natural_sounds_model_matched', '/MarmosetVocal_model_matched', '/BirdSong_model_matched'};
% directory containing the audio waveform
input_directory = [pwd Dir{2}];
% directory to save results of the synthesis process
output_directory = [pwd Dir_mm{2}];
if ~exist(output_directory,'dir')
    mkdir(output_directory);
end

fname_list_str  = dir( fullfile(input_directory, '*.wav'));
fname_list = {fname_list_str.name};

s_n = length(fname_list);

for i = 1:s_n
    fname = fname_list{i};    
    % read parameters
    P = synthesis_parameters_default;
    P.audio_sr = 100000;
    P.n_iter = 200;
    % run synthesis
    run_spectrotemporal_synthesis(P, fname, input_directory, output_directory);
end