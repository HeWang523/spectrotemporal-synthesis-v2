clc;
clear;
fs = 100000;
Dir =  '/BirdSong';
save_dir = '/BirdSong_Info/';


input_directory = [pwd Dir];

fname_list_str  = dir( fullfile(input_directory, '*.wav'));
fname_list = {fname_list_str.name};
save_directort = [pwd save_dir];


s_n = length(fname_list);

for i = 1:s_n
    [wav_orig, wav_sr] = audioread([input_directory '/' fname_list{i}]);
    fname_split = regexp(fname_list{i}, '\.', 'split');
    fname = fname_split{1};

    [pxx(i,:),f] = pspectrum(wav_orig, fs, 'FrequencyLimits', [0 25000],'FrequencyResolution',512);

    fig = figure('visible','off');
    plot(f/1000,pow2db(pxx(i,:)), 'r');
    ylim([-95, 0]);
    xlabel('Frequency (kHz)');
    ylabel('Power Spectrum (dB)');
    title("Power Spectrum of " + fname);
    set(gcf,'color','w')
    filename = save_directort + "Powerspe" + fname + ".jpg";
    saveas(fig, filename);
end


average_pxx = mean(pxx,1);
SD1=std(pow2db(pxx),1);
NSDR1=(pow2db(average_pxx)+SD1);
PSDR1=(pow2db(average_pxx)-SD1);

figure(s_n+1)
plot(f/1000,pow2db(average_pxx), 'r');
hold on
patch([f'/1000 fliplr(f'/1000)],[NSDR1 fliplr(PSDR1)], 'r', 'FaceAlpha',0.2, 'EdgeColor','none');
ylim([-95, 0]);
xlabel('Frequency (kHz)');
ylabel('Power Spectrum (dB)');
set(gcf,'color','w')
hold off





