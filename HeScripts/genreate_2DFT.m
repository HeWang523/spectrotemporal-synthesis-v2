clc
clear
% Demonstrates how to run the synthesis algorithm
% name of the audio waveform
load('C:\modelmatch\spectrotemporal-synthesis-v2-master\165_natural_sounds\category_labels.mat')
% fname = [stim_names{1} + ".wav"];
% fname = "000_TWM96.wav";

Dir = {'/165_natural_sounds', '/MarmosetVocal_filtered', '/BirdSong'};
Dir_mm = {'/165_natural_sounds_model_matched', '/MarmosetVocal_model_matched', '/BirdSong_model_matched'};
cate = 1;
% directory containing the audio waveform
input_directory = [pwd '/Agg_sounds'];
% directory to save results of the synthesis process
output_directory = [pwd '/Agg_sounds/2DFT'];

fname_list_str  = dir( fullfile(input_directory, '*.wav'));
fname_list = {fname_list_str.name};

s_n = length(fname_list);

for i = 1:s_n
    % fname = fname_list{i};   
    fname = fname_list{i};
    % read parameters
    P = synthesis_parameters_default;
    P.audio_sr = 100000;
    P.n_filts = 33;
    P.lo_freq_hz = 50;
    P.n_iter = 200;
    P.overcomplete = 0;
    % run synthesis
    func_plot_2DFT(P, fname, input_directory, output_directory);
end