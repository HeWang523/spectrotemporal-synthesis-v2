%%
clc;
clear;
%%
in_Dir = {'/MarmosetVocal_filtered', '/BirdSong_2_s'};
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
    [pxx_marmoset(i,:),f_marmoset] = compute_FFT(wav_marmoset, fs);
    % [pxx_marmoset(:,i),f_marmoset] = pspectrum(wav_marmoset, fs, 'FrequencyLimits', [f_low fs/2],'FrequencyResolution',1024);

    [ExP_marmoset(:,i), ~] = MarmosetExcitationPatternv1(f_marmoset', pxx_marmoset(i,:), f_out);
end

%% load in bird sounds and compute excitation pattern

for i = 1:bird_sound_n
    [wav_bird, ~] = audioread([bird_input_dir '/'  fname_list_bird{i}]);

    if size(wav_bird, 2) == 2
        wav_bird(:,2) = [];
    end

    wav_bird = wav_bird ./ std(wav_bird);
    [pxx_bird(i,:),f_bird] = compute_FFT(wav_bird, fs);
    % [pxx_bird(:,i),f_bird] = pspectrum(wav_bird, fs, 'FrequencyLimits', [f_low fs/2],'FrequencyResolution',1024);

    [ExP_bird(:,i), ~] = MarmosetExcitationPatternv1(f_bird', pxx_bird(i,:), f_out);
end


%% 
EXP = [ExP_marmoset ExP_bird];

% Aim: Want to get a better matching around 8 kHz
% Solution: 1. Find the idx around 8 kHz
%           2. Increase the weights for the elements of these idx
amp_factor = 5;

greater_than_8k = find(f_out > 8000);
idx_greater_than_8k = greater_than_8k(1);

less_than_8k = find(f_out < 8000);
idx_less_than_8k = less_than_8k(end);

[dim,~] = size(EXP);
temp_vector = ones(1,dim);
weight_matrix = diag(temp_vector);
weight_matrix(idx_greater_than_8k:idx_greater_than_8k+2, idx_greater_than_8k:idx_greater_than_8k+2) = amp_factor;
weight_matrix(idx_less_than_8k-2:idx_less_than_8k, idx_less_than_8k-2:idx_less_than_8k) = amp_factor;

EXP = EXP' * weight_matrix;
EXP = EXP';

% generate RDM

eu_d = pdist(EXP');
rdm= squareform(eu_d);

%% Apply MDS on RDM

[Ds,eigvals] = cmdscale(rdm, length(rdm));

eigvals = eigvals ./ sum(eigvals);

Ds = Ds(:,1:3);

figure
for j = 1:size(Ds,1)
    if j <= 15
        h1 = scatter3(Ds(j,1), Ds(j,2), Ds(j,3), 'red', 'x');
        
    else
        h2 = scatter3(Ds(j,1), Ds(j,2), Ds(j,3), 'blue', 'o');
        
    end
    hold on

end
legend([h1 h2],{"Marmoset Vocal", "Bird Songs"});
title("MDS of Marmoset Vocal and All Bird Songs");
xlabel("MDS 1");
ylabel("MDS 2");
zlabel("MDS 3");
hold off

% %% Find the center of dimension of marmoset vocal sounds
% 
% center_x = mean(Ds(1:15,1));
% center_y = mean(Ds(1:15,2));
% 
% %% compute the distance between dimension of bird songs and the center of marmoset vocal sounds
% 
% for k = 1:size(Ds,1)
%     d(k) = sqrt( (Ds(k,1)-center_x)^2 + (Ds(k,2)-center_y)^2 );
% end
% 
% %% Find 15 bird songs that have minimal distances
% [~,sort_idx] = sort(d);
% 
% sort_idx_bird_only = [];
% for k = 1:size(Ds,1)
%     if sort_idx(k) > 15
%         sort_idx_bird_only = [sort_idx_bird_only, sort_idx(k)];
%     end
% 
%     if length(sort_idx_bird_only) == 15
%         break
%     end
% end
% 
% sort_idx_bird_only = sort_idx_bird_only - 15;
% for k = 1:15
%     [wav_bird, ~] = audioread([bird_input_dir '/'  fname_list_bird{sort_idx_bird_only(k)}]);
%     if size(wav_bird, 2) == 2
%         wav_bird(:,2) = [];
%     end
%     fname = [bird_output_dir '/'  fname_list_bird{sort_idx_bird_only(k)}];
%     audiowrite(fname,wav_bird, fs);
% end
% 
%% Iteratively find the bird songs power spectrum that has the minimal distance to each marmoset vocal

marmoset_Ds = Ds(1:15,:);
birdsong_Ds = Ds(16:end,:);
birdsong_Ds_temp = birdsong_Ds;

for i = 1:length(marmoset_Ds)
    marmoset_Ds_temp = marmoset_Ds(i,:);
    x = marmoset_Ds_temp(1);
    y = marmoset_Ds_temp(2);
    z = marmoset_Ds_temp(3);
    for j = 1:length(birdsong_Ds)
        d(j) = sqrt( (birdsong_Ds(j,1)-x)^2 + (birdsong_Ds(j,2)-y)^2 + (birdsong_Ds(j,3)-z)^2 );
    end
    [sort_d,sort_idx] = sort(d);
    selected_bird_song_idx(i) = sort_idx(1);
    birdsong_Ds(sort_idx(1),:) = NaN;
end

figure
for j = 1:15
        scatter3(marmoset_Ds(j,1), marmoset_Ds(j,2), marmoset_Ds(j,3), 'red', 'x');
        hold on
        scatter3(birdsong_Ds_temp(selected_bird_song_idx(j),1), birdsong_Ds_temp(selected_bird_song_idx(j),2), birdsong_Ds_temp(selected_bird_song_idx(j),3),'blue', 'o');
        

end
title("MDS of Marmoset Vocal and Selected Bird Songs");
legend("Marmoset Vocal", "Bird Songs");
xlabel("MDS 1");
ylabel("MDS 2");
zlabel("MDS 3");
hold off


for k = 1:15
    [wav_bird, ~] = audioread([bird_input_dir '/'  fname_list_bird{selected_bird_song_idx(k)}]);
    if size(wav_bird, 2) == 2
        wav_bird(:,2) = [];
    end
    fname = [bird_output_dir '/'  fname_list_bird{selected_bird_song_idx(k)}];
    audiowrite(fname,wav_bird, fs);


    idx = find(abs(wav_bird) > 10e-4); % to remove the startup and end blank, find the indeces of wav that are greater than 0 (10e-4)

    bird_length_no_blank(k) = (idx(end) - idx(1) + 1) / fs; % store the real length of the audio
end


%% 

figure(3);
h = histogram(bird_length_no_blank,15);
h.FaceColor = [0 0.5 0.5];
h.EdgeColor = 'r';
title("Histogram of Audio of Bird Songs After Matching with Marmoset Vocal");
xlabel("Length of Audio (sec)");
ylabel("Count");












