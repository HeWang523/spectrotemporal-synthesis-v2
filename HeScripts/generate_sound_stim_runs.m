clc
clear
RUN = 1;
%% Permute the order of sound stimuli 1-390

Dir =  '/Agg_sounds';
save_dir = "/saved_runs/RUN" + RUN + "/";

input_directory = [pwd Dir];

fname_list_str  = dir( fullfile(input_directory, '*.wav'));
fname_list = {fname_list_str.name};
out_directory = [pwd + save_dir];

if ~exist(out_directory, 'dir')
    mkdir(out_directory);
end

s_n = length(fname_list);

permute_idx = randperm(s_n);

fname_list_permute = fname_list(permute_idx); % sounds are permuted

% %% Re-order the permuted sound list
%
% for i = 1:s_n
%     fname_list_reorder(permute_idx(i)) = fname_list_permute(i);
% end
%% System settings
S.System.SR =           100000;
S.System.NoAttSPL =     100;
S.Session.Artist =      ['He @ ' datestr(now, 'yyyymmdd-HH')];
S.Run.permute_idx = permute_idx; % Save the permute idx
S.Run.permuteFname = fname_list_permute; % Save the permuted sound names list
%% Runs, 30 in total
S.Session.session_n = 26; % 26 sessions in total
S.Session.trial_n = 15; % 15 trials with stimuli
S.Session.slient_trial_n = 1;
S.Session.trial_n_plus_slience = S.Session.trial_n + S.Session.slient_trial_n; % 16 trials in total (with a silent trial)
S.Session.slience_index = randi([2,S.Session.trial_n]);
S.Trial.stim_dur = 2; % stimulus duration 2s
S.Trial.repeat = 5; % Each repeats 5 times
S.Trial.interval = 0.2; % 200 ms interval between each repeated sound
S.Trial.stim_plus_interval = S.Trial.stim_dur + S.Trial.interval;
S.Trial.PreStimTime = 2;
S.Trial.StimDuration = S.Trial.repeat * S.Trial.stim_plus_interval;
S.Trial.Duration = 20; % 20 s for each trial
S.Trial.PostStimTime = S.Trial.Duration - S.Trial.PreStimTime - S.Trial.StimDuration;
S.Trial.NumberTotal = S.Session.trial_n_plus_slience; 

S.Level.TargetLevel =       70; % in dB SPL
S.Level.SLorSPL =           'SPL';


%% Blocks
sound_i = 1;
for session_i = 1:S.Session.session_n
    run_mat = zeros(S.Trial.NumberTotal, S.Trial.Duration * S.System.SR);
    run_flatten = [];
    for trial_i = 1:S.Trial.NumberTotal
        if trial_i ==  S.Session.slience_index% silent block
            continue;
        else
            [wav, ~] = audioread([input_directory '/' fname_list_permute{sound_i}]);
            wav = wav ./ (abs(max(wav)));
            wav = (wav ./ rms(wav)) * ((1/10^1.5) * sqrt(1/2));  % Normalize the original sound to 70 dB
            S.Level.RelaAtt(trial_i) = 20*log10(sqrt(1/2) / rms(wav));  % Compute the relative attenuation between the 70 dB sound and 100 dB (sqrt(1/2)), supposed to be 30 dB

            S.Trial.Att(trial_i) = S.System.NoAttSPL - S.Level.TargetLevel - S.Level.RelaAtt(trial_i); % The target SPL is set to 70 dB, the attenuations should all be 0 dB

            for rep = 1:S.Trial.repeat
                l = int64(S.System.SR * S.Trial.PreStimTime + 1 + ((rep - 1) * S.Trial.stim_plus_interval * S.System.SR));
                r = int64(l + (S.Trial.stim_dur * S.System.SR) - 1);
                run_mat(trial_i, l:r) = wav;
            end
            sound_i = sound_i + 1;
        end
    end
    for trial_i = 1:S.Trial.NumberTotal
        run_flatten = [run_flatten run_mat(trial_i,:)];
    end
    S.SoundTotalInt16 = int16(run_flatten*32767);
    S.Session.SampleTotal = length(S.SoundTotalInt16);
    S.Session.Duration =    S.Session.SampleTotal / S.System.SR;
    %% Write Run
    S.Session.Sub_Title = "390_Sounds_Nat_Vocal_Bird_Ori_Model_Matched_RUN_" + RUN + "_SESSION_" + session_i + "_";
    S.Session.TitleTime = ...
        sprintf('interval_%2.1fs_%dx(%2.1f+%2.1f+%2.1f)s_',	S.Trial.interval, S.Trial.NumberTotal,...
        S.Trial.PreStimTime, S.Trial.StimDuration ,S.Trial.PostStimTime);
    S.Session.TitleLevel = ...
        sprintf('%ddB%s_',	S.Level.TargetLevel,    S.Level.SLorSPL);
    S.Session.TitleDate =   datestr(now, 'yymmdd');
    S.Session.Title = [
        S.Session.Sub_Title + S.Session.TitleTime + S.Session.TitleLevel + S.Session.TitleDate];


    S.Session.Comment = [' ',...
        'TrialNames: ',             cell2mat(fname_list_permute(sound_i - S.Session.trial_n: sound_i - 1)),       '; ',...
        'TrialAttenuations: ',      sprintf('%3.1f ', S.Trial.Att),'; ',...
        'TrialNumberTotal: ',       num2str(S.Trial.NumberTotal ),	'; ',...
        'TrialDurTotal(sec): ',   	num2str(S.Trial.Duration),      '; ',...
        'TrialDurPreStim(sec): ',	num2str(S.Trial.PreStimTime),	'; ',...
        'TrialDurStim(sec): ',    	num2str(S.Trial.StimDuration),	'; ',...
        ''];
    run_name = out_directory + S.Session.Title + ".wav";
    audiowrite(run_name,...
        S.SoundTotalInt16,	S.System.SR,...
        'BitsPerSample',	16,...
        'Title',            S.Session.Title,...
        'Artist',           S.Session.Artist,...
        'Comment',          S.Session.Comment);
    S_name = "S_session" + session_i + ".mat";
    save(out_directory + S_name, "S");
end
