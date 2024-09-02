% Spectral Harmonicity Preference
clear all

% S.Opt.ControlTypeIndex =	1; 
S.Opt.CallTypeIndex =       4;    
S.Level.TargetLevel =       60; % in dB SPL
S.Level.SLorSPL =           'SPL';

%% Options
S.Opt.CallType{1} =     'PH';
S.Opt.CallType{2} =     'TP';
S.Opt.CallType{3} =     'TR';
S.Opt.CallType{4} =     'TW';
% S.Opt.CallType{5} =     'MX';

S.Opt.CallTypeLengthRS{1} = 236502;
S.Opt.CallTypeLengthRS{2} = 204504;
S.Opt.CallTypeLengthRS{3} = 125601;
S.Opt.CallTypeLengthRS{4} = 196502;

S.Opt.TrialStimDur{1} = 2.4;
S.Opt.TrialStimDur{2} = 2.2;
S.Opt.TrialStimDur{3} = 1.4;
S.Opt.TrialStimDur{4} = 2.0;

% S.Opt.SeedTypeSampleLength{1} = 237500; % 2.375s  x8  = 19.00s (2.365s)
% S.Opt.SeedTypeSampleLength{2} = 205000; % 2.125s  x8  = 17.00s (2.045s)
% S.Opt.SeedTypeSampleLength{3} = 131250; % 1.3125s x16 = 21.00s (1.256s)
% S.Opt.SeedTypeSampleLength{4} = 200000; % 2.000s  x8  = 20.00s (1.956s)

S.Opt.ControlType{1} =	'Reverse';
S.Opt.ControlType{2} =	'MatchSpec';
S.Opt.ControlType{3} =  'MatchEnv_lowpass_80Hz';
S.Opt.ControlType{4} =	'MatchCochEnv';
S.Opt.ControlStr = {'Rev', 'Spe', 'Env', 'Coc'};
% S.Opt.ControlType{5} =	'MatchEnv_Hilbert';
% S.Opt.ControlType{6} =	'MatRamp';

S.Opt.AllStr = [{'Ori'}, S.Opt.ControlStr];

%% Sound Seed
S.Seed.FolderParent =   'D:\Johns Hopkins\Yueqi Guo - Vocalization\';
S.Seed.FolderTarget =   'LZ_AudFilt\';
S.Seed.FolderCtrlSeed = 'LZ_ControlSounds\';

%% System settings
S.System.SR =           100000;
S.System.NoAttSPL =     100;
S.Session.Artist =      ['Xindong @ ' datestr(now, 'yyyymmdd-HH')];

%% Read Seeds
% prepare the folders
for i = 1:length(S.Opt.AllStr)
    j = i-1;
    switch S.Opt.AllStr{i}
        case 'Ori'
            S.Seed.Path{i} = [S.Seed.FolderParent S.Seed.FolderTarget ...
                S.Opt.CallType{S.Opt.CallTypeIndex} '\'];
            S.Seed.ListFiles{i} = dir([S.Seed.Path{i} '*.wav']);
        case 'Rev'
            S.Seed.Path{i} = '';
        otherwise 
            % S.Seed.FolderCtrlSeedCur = {...
            S.Seed.Path{i} = ...
                [S.Seed.FolderParent S.Seed.FolderCtrlSeed  'LZ_',...
                S.Opt.ControlType{j} '\' S.Opt.CallType{S.Opt.CallTypeIndex} '\'];
            S.Seed.ListFiles{i} = dir([S.Seed.Path{i} '*.wav']);
    end
end
% read the files
figure;
for k = 1:36     % typically 36 seeds
    for i = 1:length(S.Opt.AllStr)
        j = i-1;
        if strcmp(S.Opt.AllStr{i}, 'Rev')
             S.Seed.Sound{i,k} = flipud(S.Seed.Sound{1,k});
             S.Seed.FS(i,k) = S.Seed.FS(1,k);
        else
            [S.Seed.Sound{i,k}, S.Seed.FS(i,k)] = audioread([...
                S.Seed.Path{i}  S.Seed.ListFiles{i}(k).name]);
        end
        S.Seed.SoundLength(i,k) = length(S.Seed.Sound{i,k});
        S.Seed.SoundRS{i,k} = resample(S.Seed.Sound{i,k},...
            S.System.SR,        S.Seed.FS(i,k));
        S.Seed.SoundRS{i,k} = S.Seed.SoundRS{i,k} / std(S.Seed.SoundRS{i,k});
        S.Seed.SoundLengthRS(i,k) = length(S.Seed.SoundRS{i,k});
    end
    subplot(6,6,k);     plot(S.Seed.SoundRS{1,k});
%     set(gca, 'Xlim', [0 max(S.Seed.SoundLengthRSTarget)]);
    title([S.Seed.ListFiles{1}(k).name '  ' ...
        num2str(S.Seed.SoundLengthRS(1,k))],...
            'interpreter',      'none');
end
% display seed lengths
figure;
    plot(S.Seed.SoundLengthRS');
    legend(S.Opt.AllStr)
    title(['mean length = ' num2str(mean(S.Seed.SoundLengthRS, 'all'))]);
% sound names
    S.Seed.NameSub = cellfun(@sprintf,...
                        repmat({'S%d,'}, 36, 1), repelem(num2cell(1:6), 6)',...
                        'UniformOutput', false);
    S.Seed.NameCall = cellfun(@sprintf,...
                        repmat({'C%d,'}, 36, 1), repmat(num2cell(1:6), 1, 6)',...
                        'UniformOutput', false);
    S.Seed.NameAll = cellfun(@sprintf,...
                        repmat({'%s%s'}, 36, 1), S.Seed.NameSub, S.Seed.NameCall,...
                        'UniformOutput', false);

%% Session
S.Trial.Duration =      5.0;
S.Trial.PreStimTime = 	1.0; 
S.Trial.StimDuration =	S.Opt.TrialStimDur{S.Opt.CallTypeIndex};
S.Trial.PostStimTime = 	S.Trial.Duration - S.Trial.PreStimTime - S.Trial.StimDuration;

S.Trial.SampleTotal =       round(S.Trial.Duration*S.System.SR);    
S.Trial.SampleBase4Sound =  round(S.Trial.PreStimTime*S.System.SR);
S.Trial.Sound =             zeros(1, S.Trial.SampleTotal);
    
S.Trial.NumberTotal = length(S.Opt.AllStr)*36;

    S.SoundTotalMat =       repmat(S.Trial.Sound, S.Trial.NumberTotal, 1);
for k = 1:36     % typically 36 seeds
    for i = 1:length(S.Opt.AllStr)
        t = (k-1)*length(S.Opt.AllStr) + i;
        S.Level.RelaAtt(t) = 20*log10(max(abs(S.Seed.SoundRS{i,k}))/sqrt(1/2));
        S.SoundTotalMat(t, ...
            S.Trial.SampleBase4Sound+(1:length(S.Seed.SoundRS{i,k}))) ...
            = S.Seed.SoundRS{i,k}/max(abs(S.Seed.SoundRS{i,k}));
        S.Trial.Atts(t) = S.System.NoAttSPL - S.Level.TargetLevel - ...
                            S.Level.RelaAtt(t);
        S.Trial.Names{t} =  [S.Seed.NameAll{k} S.Opt.AllStr{i} ' '];  
    end
end
    S.SoundTotal =          reshape(S.SoundTotalMat', 1, []);
    S.SoundTotalInt16 =     int16(S.SoundTotal*32767);
    S.Session.SampleTotal = length(S.SoundTotal);
    S.Session.Duration =    S.Session.SampleTotal / S.System.SR;

%% Write the Sound
    S.Session.TitleA = ['Vocal_' S.Opt.CallType{S.Opt.CallTypeIndex} '_'];
    S.Session.TitleMod = ...
        cellfun(@sprintf, repmat({'%s,'}, length(S.Opt.AllStr), 1), S.Opt.AllStr',...
            'UniformOutput', false);
        S.Session.TitleMod = cell2mat(S.Session.TitleMod');
        S.Session.TitleMod = [num2str(length(S.Opt.AllStr)) 'x(' S.Session.TitleMod(1:end-1) ')_'];
    S.Session.TitleSC = '6xSub_6xCall_';
    S.Session.TitleTime = ...
        sprintf('%dx(%2.1f+%2.1f+%2.1f)s_',	S.Trial.NumberTotal,...
            S.Trial.PreStimTime, S.Trial.StimDuration ,S.Trial.PostStimTime);
	S.Session.TitleLevel = ...    
        sprintf('%ddB%s_',	S.Level.TargetLevel,    S.Level.SLorSPL);
    S.Session.TitleDate =   datestr(now, 'yymmdd');
S.Session.Title = [
    S.Session.TitleA,       S.Session.TitleMod,     S.Session.TitleSC,...
    S.Session.TitleTime,    S.Session.TitleLevel,   S.Session.TitleDate];
S.SoundTotalInt16 =  int16(S.SoundTotal*32767);
S.Session.Comment = [' ',...
    'TrialNames: ',             cell2mat(S.Trial.Names),       '; ',...
    'TrialAttenuations: ',      sprintf('%3.1f ', S.Trial.Atts),'; ',...
    'TrialNumberTotal: ',       num2str(S.Trial.NumberTotal),	'; ',...
    'TrialDurTotal(sec): ',   	num2str(S.Trial.Duration),      '; ',...
    'TrialDurPreStim(sec): ',	num2str(S.Trial.PreStimTime),	'; ',...
    'TrialDurStim(sec): ',    	num2str(S.Trial.StimDuration),	'; ',...
    ''];
audiowrite([ S.Session.Title,'.wav'],...
    S.SoundTotalInt16,	S.System.SR,...
    'BitsPerSample',	16,...
    'Title',            S.Session.Title,...
    'Artist',           S.Session.Artist,...
    'Comment',          S.Session.Comment);
    