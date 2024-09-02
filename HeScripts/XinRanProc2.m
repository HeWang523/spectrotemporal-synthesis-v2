function XinRanProc2(varargin)
% Xintrinsic preProcessingDATA BINNNING

global S P Tm Sys
% P:    Processed, To be saved
% Tm: 	Temporary
% Sys:  System parameters, if not in "S" yet
P = [];     Tm =[];     Sys = []; 

% P.RecCanvasHeight =     300;
% P.RecCanvasWidth =      480;
% System: Default (if not already specified in the "S")
    Tm.SysDeft.SysCamFrameRate =	80;
    Tm.SysDeft.SysCamBinNumber =	4;
    Tm.SysDeft.SysCamPixelHeight =	300;
    Tm.SysDeft.SysCamPixelWidth =	480;
        Tm.SysDeft.ProcPixelBinNum =	4;
        Tm.SysDeft.ProcFrameRate =      5;
% System: FLIR/PointGrey GS3 
    Tm.SysFlir.SysCamFrameRate =	80;
    Tm.SysFlir.SysCamBinNumber =	4;
    Tm.SysFlir.SysCamPixelHeight =	300;
    Tm.SysFlir.SysCamPixelWidth =	480;
        Tm.SysFlir.ProcPixelBinNum =	4;
        Tm.SysFlir.ProcFrameRate =      5;
%             Tm.SysFlir.ProcPixelBinNum =	1;

% System: FLIR2 /Blackfly S BFS-U3-70S7M
    Tm.SysFlir2.SysCamFrameRate =	50;
    Tm.SysFlir2.SysCamBinNumber =	4;
    Tm.SysFlir2.SysCamPixelHeight =	200;
    Tm.SysFlir2.SysCamPixelWidth =	320;
        Tm.SysFlir2.ProcPixelBinNum =	4;
        Tm.SysFlir2.ProcFrameRate =      5;
%             Tm.SysFlir2.ProcPixelBinNum =	1;
            
% System: Thorlabs sCMOS
    Tm.SysThor.SysCamFrameRate =	20;
%     Tm.SysThor.SysCamBinNumber =	4;
    Tm.SysThor.SysCamBinNumber =	6;
    Tm.SysThor.SysCamPixelHeight =	270;
    Tm.SysThor.SysCamPixelWidth =	480;
        Tm.SysThor.ProcPixelBinNum =	3;
        Tm.SysThor.ProcFrameRate =      5;

%% Get preprocessed ('*.rec') file
[~, Tm.pcname] = system('hostname');
if strcmp(Tm.pcname(1:end-1), 'FANTASIA-425')	% recording computer 
        Tm.folder = 'D:\=XINTRINSIC=\';    
else                                        % NOT recording computer
        Tm.folder = 'X:\';       
end
if nargin ==0           % Calling from direct running of the function
    Tm.RunningSource =   'D';
    [Tm.FileName, Tm.PathName, Tm.FilterIndex] = uigetfile(...
        [Tm.folder '*.rec'], 'Select raw recording files to process',...
        'MultiSelect',              'On');
    if Tm.FilterIndex == 0            
        return;                         % nothing selected
    end
    if iscell(Tm.FileName) == 0          % single file selected
        Tm.FileName = {Tm.FileName};
    end
else                    % Calling from another script
    Tm.RunningSource =   'S';
    [Tm.PathName, Tm.FileName, FileExt] = fileparts(varargin{1});
    Tm.PathName =        [Tm.PathName, '\'];
    Tm.FileName =        {[Tm.FileName, FileExt]};
end
disp(['Xintrinsic Processing Stage 1 (spatiotemporal binning) is about to start on ' ...
    num2str(length(Tm.FileName)) ' files']);

%% DATA Preprocessing for each file (Binning)
Tm.hWaitbar =	waitbar(0, 'processing');
for i = 1: length(Tm.FileName)
    % Load 'S'
    Tm.filename = [Tm.PathName, Tm.FileName{i}];
    S = load([Tm.filename(1:end-3) 'mat']);  
    try S = S.S;    end
    disp([  'Processing: "', Tm.FileName{i}, ...
            '" with the sound: "', S.SesSoundFile, '"']);
	% Default Paremeters (for files older than 2020/11/27)
                                                    Sys = Tm.SysDeft;
    if ~isfield(S, 'SysCamMain');                   S.SysCamMain = 'PointGrey';
                                                    S.SysCamDeviceName = 'Grasshopper3 GS3-U3-23S6M'; end    
	switch [S.SysCamMain '_' S.SysCamDeviceName]
       case 'PointGrey_Grasshopper3 GS3-U3-23S6M';  Sys = Tm.SysFlir;
       case 'Thorlabs_CS2100M-USB';                 Sys = Tm.SysThor;
       case 'FLIR_FLIR Blackfly S BFS-U3-70S7M';    Sys = Tm.SysFlir2;
       otherwise;                                   disp('unrecognizable camera')
	end
    if isfield(S, 'SysCamFrameRate')    % overwrite if available (for files after 2020/11/27)
        Sys.SysCamFrameRate =       S.SysCamFrameRate;
        Sys.SysCamBinNumber =       S.SysCamBinNumber;
        Sys.SysCamPixelHeight =     S.SysCamResolution(1)/S.SysCamBinNumber;
        Sys.SysCamPixelWidth =      S.SysCamResolution(2)/S.SysCamBinNumber;
    end
        Sys.SysCamFramePerTrial =	S.TrlDurTotal * Sys.SysCamFrameRate;
        Tm.SesTrlNumTotal =         length(S.SesTrlOrderVec);
    % Proc Parameters initialization for Spatial & Temporal Binning  
    P.ProcFrameRate =       Sys.ProcFrameRate;
    P.ProcFrameBinNum =     Sys.SysCamFrameRate/P.ProcFrameRate;  
    P.ProcFramePerTrial =	S.TrlDurTotal * P.ProcFrameRate;   
    P.ProcFrameNumTotal =	S.SesFrameTotal / P.ProcFrameBinNum;   
    
    P.ProcPixelBinNum =     Sys.ProcPixelBinNum;
    P.ProcCamPixelHeight =	Sys.SysCamPixelHeight/P.ProcPixelBinNum;
    P.ProcCamPixelWidth =	Sys.SysCamPixelWidth /P.ProcPixelBinNum;
    % Proc data are foced to maintain a 16:10 W/H ratio 
    %   for a consistent later visualization
    if P.ProcCamPixelWidth>(P.ProcCamPixelHeight*16/10) % original too wide
        P.ProcPixelHeight =     P.ProcCamPixelHeight;
        P.ProcPixelWidth =      round(P.ProcCamPixelHeight*16/10);
    else                                                % original too high
        P.ProcPixelHeight =     round(P.ProcCamPixelWidth*10/16);
        P.ProcPixelWidth =      P.ProcCamPixelWidth;
    end
%     P.ProcPixelHeight =     P.RecCanvasHeight/P.ProcPixelBinNum;
%     P.ProcPixelWidth =      P.RecCanvasWidth /P.ProcPixelBinNum;
    
	P.RawMeanPixel =	zeros(1, S.SesFrameTotal);
    P.RawMeanPower =	zeros(1, S.SesFrameTotal);                                
    P.ProcMeanPixel =	zeros(1, P.ProcFrameNumTotal);
    P.ProcMeanPower =	zeros(1, P.ProcFrameNumTotal);    
    P.ProcDataMat =     uint32(ones(...    
                                    S.SesCycleNumTotal,...
                                    S.TrlNumTotal,...
                                    P.ProcPixelHeight,...
                                    P.ProcPixelWidth,...
                                    P.ProcFramePerTrial...
                                    ));   
	% "ones", not "zeros": 0s at the edge would mess up normalization later
	P.ProcCamHeightIndex = (1:P.ProcPixelHeight)+round((P.ProcCamPixelHeight- P.ProcPixelHeight)/2);
	P.ProcCamWidthIndex =  (1:P.ProcPixelWidth) +round((P.ProcCamPixelWidth - P.ProcPixelWidth )/2);
% 	P.ProcDataMatHeightIndex =  (1:P.ProCamPixelHeight) + ...
%                                 round((P.ProcPixelHeight-P.ProcCamPixelHeight)/2);
% 	P.ProcDataMatWidthIndex =   (1:P.ProcCamPixelWidth) + ...
%                                 round((P.ProcPixelWidth -P.ProcCamPixelWidth )/2);
	% Downsample to 5fps and patch the dropped frames (for Thorlabs sCMOS)
        Tm.esc = PatchFrame;    
    % Read data
    if ~Tm.esc
%     Tm.fid =            fopen(Tm.filename); % Raw
            Tm.fid =            fopen([Tm.filename(1:end-4) '_' num2str(P.ProcFrameRate) 'fps.rec']); % decrease fps

    for j = 1:S.SesCycleNumTotal
        for k = 1:S.TrlNumTotal
            m = (j-1)*S.TrlNumTotal + k;
            %% Update GUI
            waitbar(m/Tm.SesTrlNumTotal, Tm.hWaitbar,...
                ['finishing ',...
                sprintf('%d out of %d total trials in the session',...
                    m, Tm.SesTrlNumTotal)] );       
            %% Read Data Batch  
%             Tm.DataRaw = 	fread(Tm.fid, [...
%                 Sys.SysCamPixelHeight * Sys.SysCamPixelWidth, ...
%                 Sys.SysCamFramePerTrial],...
%                 'uint16');
            Tm.DataProc= 	fread(Tm.fid, [...
                Sys.SysCamPixelHeight * Sys.SysCamPixelWidth, ...
                P.ProcFramePerTrial],...
                'uint32');
            %% Frame #, Trial order # location        
%             Tm.RecFramesCurrent =    ((m-1)*	Sys.SysCamFramePerTrial +1):...
%                                     (m*     Sys.SysCamFramePerTrial);
            Tm.ProcFramesCurrent =   ((m-1)*	P.ProcFramePerTrial +1):...
                                    (m*     P.ProcFramePerTrial);        
            Tm.TrialOrder =          S.SesTrlOrderVec(m);            
            %% Image Processing
            Tm.PixelMean =        mean(Tm.DataProc, 1);
            Tm.PixelMeanBinned =     Tm.PixelMean;
            P.ProcMeanPixel(Tm.ProcFramesCurrent) =  Tm.PixelMeanBinned;  
            
            switch [S.SysCamMain '_' S.SysCamDeviceName]
               case 'PointGrey_Grasshopper3 GS3-U3-23S6M'
                    Tm.ImageS0 =         reshape(Tm.DataRaw,...  
                        P.ProcPixelBinNum,     P.ProcCamPixelHeight, ...
                        P.ProcPixelBinNum,     P.ProcCamPixelWidth, ...
                        P.ProcFrameBinNum,     P.ProcFramePerTrial);
                    Tm.ImageS1 =         Tm.ImageS0;
                case 'Thorlabs_CS2100M-USB'
                    Tm.ImageS0 =         reshape(Tm.DataRaw,...  
                        P.ProcPixelBinNum,     P.ProcCamPixelWidth, ...
                        P.ProcPixelBinNum,     P.ProcCamPixelHeight, ...
                        P.ProcFrameBinNum,     P.ProcFramePerTrial);
                    Tm.ImageS1 =         permute(Tm.ImageS0, [3 4 1 2 5 6]); 
               case 'FLIR_FLIR Blackfly S BFS-U3-70S7M'
                    Tm.ImageS0 =         reshape(Tm.DataProc,...  
                        P.ProcPixelBinNum,     P.ProcCamPixelHeight, ...
                        P.ProcPixelBinNum,     P.ProcCamPixelWidth, ...
                        1,     P.ProcFramePerTrial);
                    Tm.ImageS1 =         Tm.ImageS0;
            end 
            Tm.ImageS2 =         sum(Tm.ImageS1, 1);  
            Tm.ImageS3 =         sum(Tm.ImageS2, 3); 
            Tm.ImageS4 =         sum(Tm.ImageS3, 5);
            Tm.ImageS5 =         squeeze(Tm.ImageS4);
            
            P.ProcDataMat(j, Tm.TrialOrder, :, :, :) =...
            	uint32(Tm.ImageS5(	P.ProcCamHeightIndex, ...
                                 	P.ProcCamWidthIndex, :));
            % P.ProcDataMat Dimension: 
            %   1:Cycle;    2:Trial;    3:Height;   4:Width;    5:Frame;                         
        end    
    end
    % Power Processing
        %         P.RawMeanPower =    mean(S.SesPowerMeter, 2)';
        %         P.ProcMeanPower =   mean(reshape(P.RawMeanPower,...
        %                                     P.ProcFrameBinNum,...
        %                                     P.ProcFrameNumTotal), 1 );
    % Show Figure
        Tm.timeraw =	(1:S.SesFrameTotal)/Sys.SysCamFrameRate;
        Tm.timebin =	(1:P.ProcFrameNumTotal)/P.ProcFrameRate;
        figure(     'Name',                 Tm.FileName{i},...
                    'Color',                0.8*[1 1 1]);
        Tm.hAx(1) = subplot(3,1,1);
            Tm.ColorOrder = get(gca, 'ColorOrder');
            Tm.ColorOrder = Tm.ColorOrder(2,:);
            set(Tm.hAx(1),...
                    'ColorOrder',           Tm.ColorOrder,...
                    'Toolbar',              [],...
                    'XLim',                 [0 max(Tm.timebin)],...
                    'YLim',                 8*[-1 1],...
                    'NextPlot',             'add');
%             Tm.hLineRaw =	plot( ...
%                 Tm.timeraw, 100*(P.RawMeanPixel/mean(P.RawMeanPixel)-1),...
%                     'LineWidth',            0.25);
%                 Tm.hLineRaw.Color(4) =      0.3;   
            Tm.hLineBin =	plot( ...
                Tm.timebin, 100*(P.ProcMeanPixel/mean(P.ProcMeanPixel)-1),...
                    'LineWidth',            0.75);
            Tm.hLineRaw.Color(4) =          0.8;  
        box on;
        title(Tm.hAx(1),    Tm.FileName{i},...
                    'Interpreter',          'none');
        xlabel(Tm.hAx(1),                   'Time (sec)');
        ylabel(Tm.hAx(1),               {   'Pixel Mean (%)',...
                                            sprintf('Baseline = %5.1f',...
                                                mean(P.ProcMeanPixel)) },...
                    'ButtonDownFcn',        [...
                                            'h=gcbo;  hax=h.Parent; ',...
                                            'ylim=hax.YLim;   ylim=ylim*2; ',...
                                            'if ylim(2)>20; ylim=0.125*[-1 1]; end; '...
                                            'hax.YLim=ylim; '
                                                ]);
%         [Tm.hAx,~,~] =	plotyy(	Tm.timeraw,      P.RawMeanPower, ...
%                                 Tm.timeraw,      P.RawMeanPixel);
%         xlabel(Tm.hAx(1),        'Time (sec)');
%         ylabel(Tm.hAx(1),        'Power Mean (volt)');
%         ylabel(Tm.hAx(2),        'Pixel Mean (ADU)');
%         subplot(2,1,2);
%         [Tm.hAx, Tm.hP1, Tm.hP2] = ...
%                         plotyy(	Tm.timebinned,   P.ProcMeanPower, ...
%                                 Tm.timebinned,   P.ProcMeanPixel);    
        % Tm.hP2.LineWidth =       2;                     
%         xlabel(Tm.hAx(1),        'Time (sec)');
%         ylabel(Tm.hAx(1),        'Power Mean (volt)');
%         ylabel(Tm.hAx(2),        'Pixel Mean (ADU)');
        Tm.ProcDataHWF =     squeeze(mean(mean(P.ProcDataMat, 1), 2));
        Tm.ProcDataHWpre =   squeeze(mean(Tm.ProcDataHWF(...
            :,:,1:round(S.TrlDurPreStim*P.ProcFrameRate)),3));
        Tm.ProcDataHWFnorm = Tm.ProcDataHWF./...
            repmat(Tm.ProcDataHWpre,1,1,size(Tm.ProcDataHWF,3));
        Tm.ProcDataPFprec =  (reshape(Tm.ProcDataHWFnorm, ...
            size(P.ProcDataMat,3)*size(P.ProcDataMat,4), [])-1)*100;

        Tm.hAx(2) = subplot(3,1,2);
        plot( (1:size(Tm.ProcDataPFprec,2))/P.ProcFrameRate, Tm.ProcDataPFprec' )
        set(gca,    'XTick',               [0 S.TrlDurPreStim S.TrlDurPreStim+S.TrlDurStim S.TrlDurTotal]);
        set(gca,    'XGrid',                'on');
        xlabel(Tm.hAx(2),                   'Trial time (sec)');
        ylabel(Tm.hAx(2),                   '\DeltaR/R (%)');
        set(gca,    'XLim',                [0 S.TrlDurTotal]);

        Tm.hAx(3) = subplot(3,1,3);
        plot(diff(datenum(S.SesTimestampsOri))*86400000);
        hold on
        xlim([0 length(S.SesTimestampsOri)]);
        try scatter(S.SesFrameDroppedIndices, ones(length(S.SesFrameDroppedIndices), 1)*2*1000/S.SysCamFrameRateOri + 20,'r.');
            scatter(S.SesFrameJitterIndices, ones(length(S.SesFrameJitterIndices), 1)*2*1000/S.SysCamFrameRateOri + 15,'g.');
            text(S.SesFrameQuestionIndices, ones(length(S.SesFrameQuestionIndices), 1)*2*1000/S.SysCamFrameRateOri + 10, '?', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
            yLimits = ylim;
            xLimits = xlim;
            patch([S.SesFrameOriEndIndex+1 S.SesFrameOriEndIndex+1 xLimits(2) xLimits(2)], [yLimits(1) yLimits(2) yLimits(2) yLimits(1)], [0.8 0.8 0.8], 'FaceAlpha', 0.5, 'EdgeColor', 'none');
        catch
        end
        set(gca,    'YLim',                [0 2*1000/S.SysCamFrameRateOri + 25]);
        xlabel('frame #');
        ylabel('intervel (ms)');
        title(Tm.hAx(3), ['\color{red} \bullet\color{black}(1 frm drop): ', num2str(length(S.SesFrameDroppedIndices)), ...
            ', % of total: ', sprintf('%.2f%%;' , length(S.SesFrameDroppedIndices)/length(S.SesTimestamps)*100),  ...
            '\color{green} \bullet\color{black}(jitter): ', num2str(length(S.SesFrameJitterIndices)), ...
            '; ?(others): ', num2str(length(S.SesFrameQuestionIndices))]); 
    
        figure;
        histogram(S.SesFrameIntervalMeanOf3,2.525:0.05:4.425);
        xlabel('Sum of 3 continuous intervals / Ts');
        ylabel('Counts');
        title(sprintf('Threshold of candidate patched frames = 1.25 * Ts')); 
        fprintf(['Final timestamp should be ' ...
                datestr((datenum(S.SesTimestampsOri(1,:))*86400000 + (S.SesFrameTotalOri-1) * 1000 / S.SysCamFrameRate)/86400000,'yy-mm-dd HH:MM:SS.FFF')]);
        fprintf('\n');
        fprintf(['Chosen original timestamp ' S.SesTimestampsOri(S.SesFrameOriEndIndex,:)]);
        fprintf('\n');
        % show S.SesTimestampsOri(S.SesFrameOriEndIndex+1,:)
    % Save "P"
    save([Tm.filename(1:end-4),...
        sprintf('_%dx%d@%dfps', P.ProcPixelHeight, P.ProcPixelWidth, P.ProcFrameRate),...
        '_P1.mat'], '-STRUCT', 'P', '-v7.3');     
    fclose(Tm.fid);
    end
end
%% Clean Up
close(Tm.hWaitbar);
disp('All files are processed');
clear P Tm Sys
return;

%% Downsample Original .rec and Patch Thorlabs Scientific Camera's Dropped Frames
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function esc = PatchFrame
global S Tm Sys P

%% find missing frame, modify .mat
PatchRule = 'FrameTimeStamp';
% PatchRule = 'FrameNum';
switch PatchRule(6)
    case 'N';   indexFmiss = find(S.SesFrameNum == 0);
    case 'T';   intervalF = diff(datenum(S.SesTimestamps))*86400000; % unit: ms; number: time steps - 1
                if any(intervalF > 2.5*1000/S.SysCamFrameRate) % if frame interval > 2.5 fold Ts, cannot use this program, why use 2.5? I want to see statistics
                    fprintf(['Cannot use this program to patch dropped frames!\n' ...
                        'Dropped frames > 2\n']);
                    S.SesFrameQuestionIndices = find(intervalF > 2.5*1000/S.SysCamFrameRate);
                    esc = 1;    return;
                end
                indexFmiss = find(intervalF > 1.25*1000/S.SysCamFrameRate); % candidate index, whu use 1.25? new workstation is ok, how about old workstation?

                % continuous 3 frames
                indexFmisslast = indexFmiss;      % copy indexFmiss
                indexFmisslast(indexFmiss == 1) = 3;  % if indexFmiss == 1, then indexFmisslast = 3
                indexFmisslast(indexFmiss ~= 1) = indexFmiss(indexFmiss ~= 1) - 1;  % other elements - 1
                indexFmissnext = indexFmiss;      % copy indexFmiss
                indexFmissnext(indexFmiss == length(intervalF)) = length(intervalF)-2;  % if indexFmiss == 1, then indexFmisslast = length(intervalF)-2
                indexFmissnext(indexFmiss ~= length(intervalF)) = indexFmiss(indexFmiss ~= length(intervalF)) + 1;  % other elements - 1
        
                S.SesFrameIntervalMeanOf3 = sum(intervalF([indexFmisslast, indexFmiss, indexFmissnext]), 2)/(1000/S.SysCamFrameRate);
             
                S.SesFrameJitterIndices  = indexFmiss(abs(S.SesFrameIntervalMeanOf3-3) <= 0.11);
                S.SesFrameDroppedIndices = indexFmiss(abs(S.SesFrameIntervalMeanOf3-4) <= 0.11); 
                S.SesFrameQuestionIndices  = setdiff(indexFmiss, union(S.SesFrameJitterIndices,S.SesFrameDroppedIndices));

                indexFmiss = S.SesFrameDroppedIndices;
                indexFmiss = ((1:length(indexFmiss))+indexFmiss')'; % final missing frame index to patch
                indexFmiss = indexFmiss(indexFmiss <= S.SesFrameTotal); % No considering index number more than Session frame num
                patchedNum = length(indexFmiss);
                S.SesFrameOriEndIndex = S.SesFrameTotal - patchedNum;
                S.SesTimestampsLast1Diff = 1 + (datenum(S.SesTimestamps(S.SesFrameOriEndIndex,:)) - datenum(S.SesTimestamps(1,:))) * 86400000 / (1000/S.SysCamFrameRate) - S.SesFrameTotal; % unit: frame
                S.SesFrameNumOri = S.SesFrameNum;
                S.SesFrameNum = (1:length(S.SesFrameNumOri))'; 
                S.SesFrameNum(indexFmiss) = 0;
end
if isempty(indexFmiss)   % no missing frame
        ProcFrameBinNum =     S.SysCamFrameRate/Sys.ProcFrameRate;
        unitperframe = Sys.SysCamPixelHeight * Sys.SysCamPixelWidth;
        fid_ori = fopen([Tm.filename(1:end-4) '.rec']);
        fid_ptd = fopen( [Tm.filename(1:end-4) '_' num2str(P.ProcFrameRate) 'fps.rec'], 'w');
            for i = 1:ProcFrameBinNum:length(S.SesFrameNum)
                aggregated_frame = zeros(unitperframe, 1, 'uint32'); 
            
                for j = 0:(ProcFrameBinNum-1)
                    current_index = i + j;
                    if current_index > length(S.SesFrameNum)
                        break;
                    end
            
                    % read frame
                    if S.SysCamFrameRate ~= Sys.ProcFrameRate
                        frame =         uint32(fread(fid_ori, [unitperframe, 1], 'uint16'));
                    else
                        frame =         uint32(fread(fid_ori, [unitperframe, 1], 'uint32'));
                    end

                    aggregated_frame = aggregated_frame + frame;
                end
                % write frame
                fwrite(fid_ptd, aggregated_frame, 'uint32');
                % S.SesFrameTimestamps() = S.SesFrameTimestamps;
            end    
    
            S.SesTimestampsOri = S.SesTimestamps;
            S.SesTimestamps = S.SesTimestamps(1:ProcFrameBinNum:length(S.SesFrameNum),:);
            S.SysCamFrameRateOri = S.SysCamFrameRate;
            S.SysCamFrameRate = Sys.ProcFrameRate;
            S.SesFrameTotalOri = S.SesFrameTotal;
        
        
            S.SesFrameTotal = S.SesFrameTotalOri / ProcFrameBinNum;
            S.SesFrameNum = (1:S.SesFrameTotal)'; 
            fclose(fid_ori);
            fclose(fid_ptd);
            save([Tm.filename(1:end-4) '_' num2str(P.ProcFrameRate) 'fps.mat'], 'S', '-v7.3');
            fprintf('\n');
            esc = 0;    return;
else
    esc = 0;
end                         % forward is for patch the missing frames
% available_frame = S.SesFrameNum;
indexFlast = find(([1 S.SesFrameNum(1:end-1)'] - S.SesFrameNum')>0)-1;
indexFnext = find(([S.SesFrameNum(2:end)' length(S.SesFrameNum)+1] - S.SesFrameNum')>1)+1;

disp(['missing frames#: ' sprintf('%d ', indexFmiss)]);
% fprintf('\n');
if abs(S.SesTimestampsLast1Diff) < 0.5
    namestr = sprintf('_%d_missingframes', length(indexFmiss));
    fprintf('S.SesTimestampsLast1Diff = %.2f',S.SesTimestampsLast1Diff);
    fprintf('\n');
else
    namestr = sprintf('_%d_missingframes_QM', length(indexFmiss));
    fprintf('QM! S.SesTimestampsLast1Diff = %.2f',S.SesTimestampsLast1Diff);
    fprintf('\n');
end
%% Save Original files & ...
movefile(    Tm.filename, ...
            [Tm.filename(1:end-4) namestr '.rec']);
movefile(   [Tm.filename(1:end-4) '.mat'], ...
            [Tm.filename(1:end-4) namestr '.mat']);
fid_ori = fopen([Tm.filename(1:end-4) namestr '.rec']);
fid_ptd = fopen( [Tm.filename(1:end-4) '_' num2str(P.ProcFrameRate) 'fps.rec'], 'w');
%% Patch the ends if needed
ProcFrameBinNum =     S.SysCamFrameRate/Sys.ProcFrameRate;  
unitperframe = Sys.SysCamPixelHeight * Sys.SysCamPixelWidth;
% unitperframe = S.SysCamResolution(1) / S.SysCamBinNumber * S.SysCamResolution(2) / S.SysCamBinNumber;

if indexFmiss(1) == 1
    % read the current "next" frame but also as the "0" frame
    indexCl =   indexFlast(1);
    frameCl =   uint32(fread(fid_ori, [unitperframe, 1], 'uint16'));
    % the current "next" frame;
    indexCn =   indexFnext(1);
    frameCn =   frameCl;
end


S.SesTimestampsOri = S.SesTimestamps;
for i = 1:ProcFrameBinNum:length(S.SesFrameNum)
    aggregated_frame = zeros(unitperframe, 1, 'uint32'); 

    for j = 0:(ProcFrameBinNum-1)
        current_index = i + j;
        if current_index > length(S.SesFrameNum)
            break;
        end

        % read frame
        if      ismember(current_index, indexFlast) 
            if ~ismember(current_index, indexFnext)	% last frame before gap, but already read
                indexCl =   current_index;
                frameCl =	uint32(fread(fid_ori, [unitperframe, 1], 'uint16'));
            else                        % last frame before gap, but already read    
                indexCl =   current_index;
                frameCl =	frameCn;    % as the "next" frame for the last gap                                   
            end
                indexCn =   indexFnext( find(indexFnext>current_index,1) );
            if indexCn > length(S.SesFrameNum)
                frameCn =	frameCl;    % ending frame is missing
            else                        % read the "next" frame after the gap
                frameCn =	uint32(fread(fid_ori, [unitperframe, 1], 'uint16'));
            end
            frame =	frameCl;
        elseif  ismember(current_index, indexFnext) % pure "next" frame after a gap
            frame = frameCn;
        elseif  ismember(current_index, indexFmiss) % missed frame
            frame = frameCl*(current_index-indexCl)/(indexCn-indexCl) + ...
                    frameCn*(indexCn-current_index)/(indexCn-indexCl); % patch the frame
                        S.SesFrameNum(current_index) = current_index;
            fprintf('frame#%d is patched from the "last" frame #%d and the "next" frame #%d\n',...
                    current_index, indexCl, indexCn);
            beforeInsert = S.SesTimestamps(1:current_index-1,:);
            afterInsert = S.SesTimestamps(current_index:end,:);
            S.SesTimestamps = [beforeInsert; 
                datestr((datenum(S.SesTimestamps(current_index-1,:),'yy-mm-dd HH:MM:SS.FFF')+datenum(S.SesTimestamps(current_index,:),'yy-mm-dd HH:MM:SS.FFF'))/2, ...
                'yy-mm-dd HH:MM:SS.FFF'); 
                afterInsert];

        else                            % ordinary frame
            frame =         uint32(fread(fid_ori, [unitperframe, 1], 'uint16'));
        end
        aggregated_frame = aggregated_frame + frame;
    end
    % write frame
    fwrite(fid_ptd, aggregated_frame, 'uint32');
    % S.SesFrameTimestamps() = S.SesFrameTimestamps;
end
S.SesTimestamps = S.SesTimestamps(1:ProcFrameBinNum:length(S.SesFrameNum),:);
S.SysCamFrameRateOri = S.SysCamFrameRate;
S.SysCamFrameRate = Sys.ProcFrameRate;
S.SesFrameTotalOri = S.SesFrameTotal;
S.SesFrameTotal = S.SesFrameTotalOri / ProcFrameBinNum;
S.SesFrameNum = (1:S.SesFrameTotal)'; 


%% Finish writing and saving
fclose(fid_ori);
fclose(fid_ptd);
save([Tm.filename(1:end-4) '_' num2str(P.ProcFrameRate) 'fps.mat'], 'S', '-v7.3');
fprintf('\n');
