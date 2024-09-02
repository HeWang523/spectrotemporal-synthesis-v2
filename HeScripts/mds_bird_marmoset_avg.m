%%
clc;
clear;
%%
in_Dir = {'/MarmosetVocal_filtered', '/Bird_song_resample'};
out_Dir = {'/BirdSong_selected'};

marmoset_input_dir = [pwd in_Dir{1}];
bird_input_dir = [pwd in_Dir{2}];

bird_output_dir = [pwd out_Dir{1}];

fname_list_str_marmoset  = dir( fullfile(marmoset_input_dir, '*.wav'));
fname_list_marmoset = {fname_list_str_marmoset.name};

fname_list_str_bird  = dir( fullfile(bird_input_dir, '*.wav'));
fname_list_bird = {fname_list_str_bird.name};


marmoset_sound_n = length(fname_list_marmoset);
bird_sound_n = length(fname_list_bird);

%% load in marmoset vocal sounds and compute excitation pattern
fs = 50000;
f_low = 50;

f_out = f_low*2.^(0:log2(fs / (f_low * 2))/160:log2(fs / (f_low * 2)));


for i = 1:marmoset_sound_n
    [wav_marmoset, ~] = audioread([marmoset_input_dir '/'  fname_list_marmoset{i}]);
    wav_marmoset = wav_marmoset ./ std(wav_marmoset);
    [pxx_marmoset(:,i),f_marmoset] = pspectrum(wav_marmoset, fs, 'FrequencyLimits', [f_low fs/2],'FrequencyResolution',128);

    
end
[ExP_marmoset(:,1), ~] = MarmosetExcitationPatternv1(f_marmoset, mean(pxx_marmoset,2), f_out);
%% load in bird sounds and compute excitation pattern

for i = 1:bird_sound_n
    [wav_bird, ~] = audioread([bird_input_dir '/'  fname_list_bird{i}]);

    if size(wav_bird, 2) == 2
        wav_bird(:,2) = [];
    end

    wav_bird = wav_bird ./ std(wav_bird);
    [pxx_bird(:,i),f_bird] = pspectrum(wav_bird, fs, 'FrequencyLimits', [f_low fs/2],'FrequencyResolution',128);

    [ExP_bird(:,i), ~] = MarmosetExcitationPatternv1(f_bird, pxx_bird(:,i), f_out);
end


%% 
EXP = [mean(ExP_marmoset,2) ExP_bird];

% generate RDM

eu_d = pdist(EXP');
rdm= squareform(eu_d);

%% Apply MDS on RDM

[Ds,eigvals] = cmdscale(rdm, length(rdm));

eigvals = eigvals ./ sum(eigvals);

Ds = Ds(:,1:2);

figure
for j = 1:size(Ds,1)
    if j <= 1
        scatter(Ds(j,1), Ds(j,2), 'red', 'x');
    else
        scatter(Ds(j,1), Ds(j,2), 'blue', 'o');
    end
    hold on

end
hold off

%% Find the center of dimension of marmoset vocal sounds

center_x = mean(Ds(1,1));
center_y = mean(Ds(1,2));

%% compute the distance between dimension of bird songs and the center of marmoset vocal sounds

for k = 1:size(Ds,1)
    d(k) = sqrt( (Ds(k,1)-center_x)^2 + (Ds(k,2)-center_y)^2 );
end

%% Find 15 bird songs that have minimal distances
[~,sort_idx] = sort(d);

sort_idx_bird_only = [];
for k = 1:size(Ds,1)
    if sort_idx(k) > 15
        sort_idx_bird_only = [sort_idx_bird_only, sort_idx(k)];
    end

    if length(sort_idx_bird_only) == 15
        break
    end
end

sort_idx_bird_only = sort_idx_bird_only - 1;
for k = 1:15
    [wav_bird, ~] = audioread([bird_input_dir '/'  fname_list_bird{sort_idx_bird_only(k)}]);
    if size(wav_bird, 2) == 2
        wav_bird(:,2) = [];
    end
    fname = [bird_output_dir '/'  fname_list_bird{sort_idx_bird_only(k)}];
    audiowrite(fname,wav_bird, fs);
end

% %% Iteratively find the bird songs power spectrum that has the minimal distance to each marmoset vocal
% 
% marmoset_Ds = Ds(1:15,:);
% birdsong_Ds = Ds(16:end,:);
% birdsong_Ds_temp = birdsong_Ds;
% 
% for i = 1:length(marmoset_Ds)
%     marmoset_Ds_temp = marmoset_Ds(i,:);
%     x = marmoset_Ds_temp(1);
%     y = marmoset_Ds_temp(2);
%     for j = 1:length(birdsong_Ds)
%         d(j) = sqrt( (birdsong_Ds(j,1)-x)^2 + (birdsong_Ds(j,2)-y)^2 );
%     end
%     [sort_d,sort_idx] = sort(d);
%     selected_bird_song_idx(i) = sort_idx(1);
%     birdsong_Ds(sort_idx(1),:) = NaN;
% end
% 
% figure
% for j = 1:15
%         scatter(marmoset_Ds(j,1), marmoset_Ds(j,2), 'red', 'x');
%         hold on
%         scatter(birdsong_Ds_temp(selected_bird_song_idx(j),1), birdsong_Ds_temp(selected_bird_song_idx(j),2), 'blue', 'o');
% 
% end
% hold off
% 
% 
% for k = 1:15
%     [wav_bird, ~] = audioread([bird_input_dir '/'  fname_list_bird{selected_bird_song_idx(k)}]);
%     if size(wav_bird, 2) == 2
%         wav_bird(:,2) = [];
%     end
%     fname = [bird_output_dir '/'  fname_list_bird{selected_bird_song_idx(k)}];
%     audiowrite(fname,wav_bird, fs);
% end


%% 














