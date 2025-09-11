function [data] = BranchContentContextTMS(subNum, visit)
% Line 24, 25, 224, 232, 243 + get rid of extra cases in menu
if nargin < 1

    subNum = input('Participant number:');
    visit = input('Visit number:'); 
    
end

%% Experiment variables
RUN_DURATION = 574; %length of a run in seconds
TRIALS_PER_RUN = 144; %number of trials per run, divisible by 3
BLOCKS_PER_RUN = 16; %number of blocks per run, divisible by 2
NUM_RUNS = 4; %number of runs in the session
BASE_RT = 1; %RT to beat to get bonus
LEFT_RESPONSE = 'f'; %left response key
RIGHT_RESPONSE = 'j'; %right response key
SEQUENCE_LENGTH = 5; %number of stimuli in a sequence
STIM_DUR = 0.5; %Stimulus duration
FB_SLIDE = 0.3; %Feedback time (correct/incorrect)
%STARTING_FIX = 0.5; %Fixation time to start each trial for timing purposes

menu_keys = [KbName('p') KbName('e') KbName('q') KbName('i')];
select_run = [KbName('y') KbName('1') KbName('1!') KbName('2') KbName('2@') KbName('3') KbName('3#') KbName('4') KbName('4$')];

%color-to-content mappings
 VerbalColor = 'Blue';
 SpatialColor = 'Yellow';

%shape-to-task mappings
 DelayShape = 'Circle';
 SwitchShape = 'Diamond';
 BranchShape = 'Cross';

%decision-response mappings
  YesKey = LEFT_RESPONSE;
  NoKey = RIGHT_RESPONSE;


%% Lists Jittered ITI and IBI
TrialJitters = zeros(TRIALS_PER_RUN,1);
TrialJitters(1:(TRIALS_PER_RUN/3)) = 2.6;
TrialJitters((TRIALS_PER_RUN/3)+1:(2*TRIALS_PER_RUN/3)) = 3;
TrialJitters((2*TRIALS_PER_RUN/3)+1:TRIALS_PER_RUN) = 3.4;

BlockJitters = zeros(BLOCKS_PER_RUN,1);
BlockJitters(1:(BLOCKS_PER_RUN/2)) = 2.6;
BlockJitters((BLOCKS_PER_RUN/2)+1:BLOCKS_PER_RUN) = 3.4;

[MatchList{1:(TRIALS_PER_RUN/2)}] = deal('Match');
[MatchList{(TRIALS_PER_RUN/2 + 1):TRIALS_PER_RUN}] = deal('NonMatch');

%% Block Parameters
%% Block Parameters
ListNum = mod(subNum,50) + 1;
% Use fullfile to ensure the path is correct across different operating systems.
% Replace the path with the correct one for your computer if different.
filename = ['/Users/aniladmello/Desktop/MATLAB/TMS/StimulusList_ContextTMS' num2str(ListNum) '.txt'];

% Attempt to open the file and capture a potential error message.
[fid, errmsg] = fopen(filename, 'r'); % 'r' for read permission

% Check if the file was successfully opened. A valid fid is >= 3.
if fid < 0
    % If fid is -1, the file could not be opened. Display a helpful error.
    error('File "%s" could not be opened. MATLAB error: %s', filename, errmsg);
end

% The file is open, so it's safe to use textscan.
d = textscan(fid, '%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\r', 'HeaderLines', 1);

% Close the file after you are done reading from it.
fclose(fid);

%ListNum = mod(subNum,50) + 1;
%filename = ['\Users\aniladmello\Desktop\MATLAB\TMS\StimulusList_ContextTMS' num2str(ListNum) '.txt'];
%fid = fopen(filename);
%d = textscan(fid,'%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\r','HeaderLines',1);
%fclose(fid);

ContentList = d{2};
WMList = d{3};
DTList = d{4};
BlockTypeList = d{5};
numTrialsList = [d{6} d{7} d{8}];


%data structure
trialStruct = struct('RESP',-1,'RT',-1,'ACC',-1,'Onset',-1,'OnsetError',-1, ...
           'Duration',-1,'DurationError',-1, ...
           'TrialType','NULL','MatchType','NULL','Letter','NULL', ...
           'Location',[-1 -1 -1 -1],'Content','NULL','WM','NULL', ...
           'DT','NULL','RunNum',-1,'BlockNum',-1,'TrialNum',-1);

data = repmat(trialStruct,NUM_RUNS*TRIALS_PER_RUN,1);

%% Initialization
Screen('Preference','VisualDebugLevel',1);
Screen('Preference','SkipSyncTests',1);
Screen('Preference', 'TextAntiAliasing', 1);

rng shuffle % creates a new random seed for this session

%Going to start KbQueue, will need a way to stop
KbName('UnifyKeyNames');
escapeKey = KbName('ESCAPE');

KbQueueCreate(-1); %create the queue
releaseQueue = onCleanup(@() KbQueueRelease);

KbQueueStart; %start the queue

ListenChar(2);
restoreKeyboard = onCleanup(@() ListenChar(0));

HideCursor;
restoreCursor = onCleanup(@() ShowCursor);

% Scanner Projector should be 1920 x 1080
%Res = Screen('Resolution',0);
screens = Screen('Screens');
screenNumber = screens(1);
% Res = Screen('Resolution', screenNumber, 1920, 1080, 60);
window = Screen('OpenWindow',screenNumber,[0 0 0]); 
Screen('TextSize', window,max(24, Res.width/80));

closeWindow = onCleanup(@() Screen('CloseAll'));

Priority(MaxPriority(0));

restorePriority = onCleanup(@() Priority(0));

FlipDur = Screen('GetFlipInterval', window);

%code for changing resolution and centering
%oldRes = Screen('Resolution',0);
%newRes.width = 1920; newRes.height = 1080;
%windowRect = CenterRect([0 0 newRes.width newRes.height],[0 0 oldRes.width oldRes.height]);
%window = Screen('OpenWindow',0,[0 0 0],windowRect);

%% Load Stimuli as Textures

FrameColors = {'Blue' 'Yellow'};
FrameShapes = {'Circle' 'Cross' 'Diamond' 'Square'};
Letters = {'T' 'A' 'B' 'L' 'E'};
StimDir = '/Users/aniladmello/Desktop/MATLAB/stimulus';
TextureMat = zeros(length(FrameColors),length(FrameShapes),length(Letters));

for colori = 1:length(FrameColors)
    for shapei = 1:length(FrameShapes)
        for letteri = 1:length(Letters)
            StimName = [FrameColors{colori} FrameShapes{shapei} Letters{letteri}];
            fileStr = sprintf('%s%s.png',StimDir,StimName);
            imagedata = imread(fileStr);
            TextureMat(colori,shapei,letteri) = Screen('MakeTexture',window,imagedata);
        end
    end
end

%Location centers in 1024 x 768 space
%LocationCenters = {[512 274] [577 473] [407 350] [617 350] [447 473]};

%Location centers in general
LocationCenters = {[Res.width/2 Res.height/2 - 110] ...
                   [Res.width/2 + 65 Res.height/2 + 89] ...
                   [Res.width/2 - 105 Res.height/2 - 34] ...
                   [Res.width/2 + 105 Res.height/2 - 34] ...
                   [Res.width/2 - 65 Res.height/2 + 89]};

imSize = 100; %each image is 100 x 100 pixels

Locations = cell(5,1);
for centeri = 1:length(LocationCenters)
    Locations{centeri} = [LocationCenters{centeri}-(imSize/2) LocationCenters{centeri}+(imSize/2)];
end

%% Fixation cross size and location
Screen('TextSize',window,40);
Screen('TextStyle',window,1); %Bold 
Screen('TextFont',window,'Courier New');
[bounds] = Screen('TextBounds',window,'+');
FixXAdjust = bounds(3)/2;
FixYAdjust = bounds(4)/2;


%Stimuli Check
%{
Screen('DrawTexture',window,BlueSquareT,[],Locations{1});
Screen('DrawTexture',window,BlueSquareA,[],Locations{2});
Screen('DrawTexture',window,BlueSquareB,[],Locations{3});
Screen('DrawTexture',window,BlueSquareL,[],Locations{4});
Screen('DrawTexture',window,BlueSquareE,[],Locations{5});
Screen('DrawText',window,'+',Res.width/2 - FixXAdjust,Res.height/2 - FixYAdjust,[255 255 255]);
Screen('Flip',window);
KbWait;
%}

%% Begin Paradigm/ Main Menu
QuitExp = 0;
RunSuggestion = 1;
scoreData.Score = 0;
scoreData.CurrACC = 0;
scoreData.CurrRT = 0;
scoreData.LastACC = 0;
scoreData.LastRT = 0;
scoreData.TotalACC = 0;
scoreData.TotalRT = 0;
scoreData.TotalTrials = 0;

%Create OnsetFile
%onsetfilename = sprintf('OnsetTimesfMRI_%d_%d',subNum, visit);
%OnsetFile = fopen(onsetfilename,'a');
OnsetFile = fopen(['NewCogFlex' num2str(subNum) '.txt'],'w');
closeOnsetFile = onCleanup(@() fclose(OnsetFile));

while ~QuitExp
    MenuStr = sprintf('Main Menu\n\n\n(I)nstructions\n\n(P)ractice\n\n(E)xperiment\n\n(Q)uit');
    DrawFormattedText(window,MenuStr,'center','center',[255 255 255]);

    Screen('Flip',window);
    
    RestrictKeysForKbCheck(menu_keys);
    [~,keycode] = KbWait([],2);
    RestrictKeysForKbCheck([]);
    
    switch KbName(keycode)
        case 'p'
            PracticeProc(scoreData);
        case 'e'
            tic;
            RunStr = sprintf('Hit (a) to continue to run %d or enter run number',RunSuggestion);
            [bounds] = Screen('TextBounds',window,RunStr);
            Screen('DrawText',window,RunStr,1920/2 - bounds(3)/2,1080/2 - bounds(4)/2,[255 255 255]);
            Screen('Flip',window);
            
            RestrictKeysForKbCheck(select_run);
            [~,keycode] = KbWait([],2);
            RestrictKeysForKbCheck([]);

            keyname = KbName(keycode);
            if keyname=='a'
                [data,scoreData] = RunProc(RunSuggestion,data,scoreData);
            else
                RunNum = str2num(keyname(1));
                if ismember(RunNum,1:NUM_RUNS)
                    [data,scoreData] = RunProc(RunNum,data,scoreData);
                else
                    while ~any(ismember(RunNum,1:NUM_RUNS)) && ~isequal(keyname,'y')
                        [~,keycode] = KbWait([],2);
                        keyname = KbName(keycode);
                        RunNum = str2num(keyname(1));
                    end
                    if KbName(keycode)=='y'
                        [data,scoreData] = RunProc(RunSuggestion,data,scoreData);
                    else
                        [data,scoreData] = RunProc(RunNum,data,scoreData);
                    end
                end
            end
            
            RestProc(scoreData);
            elapsedExTime = toc;
            disp(['Elapsed experiment time: ', num2str(elapsedExTime), ' seconds']);
            
            RunSuggestion = RunSuggestion + 1;
            if RunSuggestion > NUM_RUNS
                break;
            end
        case 'q'
            break;
        case 'i'
            Instructions();
        otherwise
            %do nothing, just loop
    end
end

%% Run Procedure
    function [data,expData] = RunProc(runi,data,expData)
        expData.EpochDur = 0;
        expData.EpochNum = 0;
        expData.EpochCorr = 0;
        expData.EpochType = '';
        expData.LastITI = 0;
        
        BlockIdx = (runi-1)*BLOCKS_PER_RUN + 1;
        TrialIdx = ((runi-1)*TRIALS_PER_RUN + 1):(runi*TRIALS_PER_RUN);
        runData = data(TrialIdx);

        BlockInfo.TrialNum = 1;
        RunInfo.TrialJitters = Shuffle(TrialJitters);
        RunInfo.BlockJitters = Shuffle(BlockJitters);
        RunInfo.MatchList = Shuffle(MatchList);

        %Wait for scanner ready
               Instr = sprintf(['Welcome! We will now begin the behavioral task practice. Remember the rules:\n\n'...
                        '1) The letter task is indicated by %s colored frames and the\n'...
                        'location task is indicated by %s colored frames\n\n'...
                        '2) The Basic task is indicated by Square frames\n\n'...
                        '3) The Delay task is indicated by %s frames\n\n'...
                        '4) The Restart task is indicated by %s frames\n\n'...
                        '5) The Dual task is indicated by %s frames\n\n'...
                        'To respond Yes press %s and to respond No press %s\n\n'...
                        '(Press any key to begin experiment)'],VerbalColor, SpatialColor,DelayShape,SwitchShape,BranchShape,YesKey,NoKey);
        DrawFormattedText(window,Instr,'center','center',[255 255 255]);
        Screen('Flip',window);
        KbWait([],2);

        [bounds] = Screen('TextBounds',window,'Get Ready!!!');
        Screen('DrawText',window,'Get Ready!!!',Res.width/2 - bounds(3)/2,Res.height/2 - bounds(4)/2,[255 255 255]);
        Screen('DrawingFinished', window);
       

        %Wait for trigger
        RunInfo.RunStart = KbWait;
        RunInfo.MyClock = RunInfo.RunStart;
        Screen('Flip',window);
     
        RunInfo.MyClock = RunInfo.MyClock + 3; %3 second ready cue

        Screen('DrawText',window,'+',Res.width/2 - FixXAdjust,Res.height/2 - FixYAdjust,[255 255 255]);
        Screen('Flip',window,RunInfo.MyClock-(FlipDur/2));   
        RunInfo.MyClock = RunInfo.MyClock + 1; %1 second post-cue fixation

        for blocki = 1:BLOCKS_PER_RUN
            %Choose Block Variables
            BlockInfo.BlockType = BlockTypeList{BlockIdx};
            BlockInfo.numTrials = numTrialsList(BlockIdx,:);
            BlockInfo.Content = ContentList{BlockIdx};
            BlockInfo.WM = WMList{BlockIdx};
            BlockInfo.DT = DTList{BlockIdx};
            BlockInfo.IBI = RunInfo.BlockJitters(blocki);
            BlockInfo.RunNum = runi;
            BlockInfo.BlockNum = blocki;

            [runData,BlockInfo,RunInfo,expData] = BlockProcedure(runData,BlockInfo,RunInfo,expData,0);

            BlockIdx = BlockIdx + 1;
        end

        %Wait about 10 seconds for HRF recovery
        RunEnd = WaitSecs('UntilTime',RunInfo.MyClock+10);

        disp(['Run duration Expected: ' num2str(RUN_DURATION)]); 
        disp(['Run duration Actual: ' num2str(RunEnd-RunInfo.RunStart)]);
        data(TrialIdx) = runData;

        %datafilename = sprintf('datafMRI_%d_%d',subNum, visit);
        %save(datafilename,'data'); %save at end of each run
        save(['dataBehav_' num2str(subNum)],'data'); %save at end of each run

                
        expData.LastACC = expData.CurrACC;
        expData.LastRT = expData.CurrRT;
        expData.CurrACC = sum([runData.ACC]);
        expData.CurrRT = sum([runData.RT]);
        expData.TotalACC = expData.TotalACC + expData.CurrACC;
        expData.TotalRT = expData.TotalRT + expData.CurrRT;
        expData.TotalTrials = expData.TotalTrials + TRIALS_PER_RUN;
    end

%% GetNewStimulus
    function [NewLetter,NewLocation] = GetNewStimulus(MatchType,Content,PhaseCase,LastLetter,LastLocation,DelayLetter,DelayLocation)
        switch PhaseCase
            case 0 %continue a sequence
                if strcmp(Content,'verbal')
                    NextLetter = LastLetter + 1;
                    if NextLetter > SEQUENCE_LENGTH
                        NextLetter = 1;
                    end
                    NextLocation = 0;
                else
                    NextLocation = LastLocation + 1;
                    if NextLocation > SEQUENCE_LENGTH
                        NextLocation = 1;
                    end
                    NextLetter = 0;
                end
                
                if strcmp(MatchType,'Match')
                    if strcmp(Content,'verbal')
                        %match verbal trial, take next letter in sequence
                        NewLetter = NextLetter;
                        
                        %get random new location that does not match last
                        %location
                        cand = 1:SEQUENCE_LENGTH;
                        cand = cand(cand ~= LastLocation);
                        NewLocation = cand(randi(length(cand),1));
                    else
                        %match spatial trial, take next locatoin in
                        %sequence
                        NewLocation = NextLocation;
                        
                        %get random new letter that does not match last
                        %letter
                        cand = 1:SEQUENCE_LENGTH;
                        cand = cand(cand ~= LastLetter);
                        NewLetter = cand(randi(length(cand),1));
                    end
                else
                    %NonMatch trial, get random letters and locations that
                    %do not match last trial or next in sequence
                    if strcmp(Content,'verbal')
                        cand = 1:SEQUENCE_LENGTH;
                        cand = cand(cand ~= LastLetter & cand ~= NextLetter);
                        NewLetter = cand(randi(length(cand),1));
                        
                        cand = 1:SEQUENCE_LENGTH;
                        cand = cand(cand ~= LastLocation);
                        NewLocation = cand(randi(length(cand),1));
                    else
                        cand = 1:SEQUENCE_LENGTH;
                        cand = cand(cand ~= LastLetter);
                        NewLetter = cand(randi(length(cand),1));
                        
                        cand = 1:SEQUENCE_LENGTH;
                        cand = cand(cand ~= LastLocation & cand ~= NextLocation);
                        NewLocation = cand(randi(length(cand),1));
                    end
                end
            case 1 %start a new sequence
                if strcmp(MatchType,'Match')
                    if strcmp(Content,'verbal')
                        %match verbal trial, take start of sequence
                        NewLetter = 1;
                        
                        %get random new location that does not match last
                        %location
                        cand = 1:SEQUENCE_LENGTH;
                        cand = cand(cand ~= LastLocation);
                        NewLocation = cand(randi(length(cand),1));
                    else
                        %match spatial trial, take start of sequence
                        NewLocation = 1;
                        
                        %get random new letter that does not match last
                        %letter
                        cand = 1:SEQUENCE_LENGTH;
                        cand = cand(cand ~= LastLetter);
                        NewLetter = cand(randi(length(cand),1));
                    end
                else
                    %NonMatch trial, get random letters and locations that
                    %do not match last trial or first position
                    if strcmp(Content,'verbal')
                        cand = 1:SEQUENCE_LENGTH;
                        cand = cand(cand ~= LastLetter & cand ~= 1);
                        NewLetter = cand(randi(length(cand),1));
                        
                        cand = 1:SEQUENCE_LENGTH;
                        cand = cand(cand ~= LastLocation);
                        NewLocation = cand(randi(length(cand),1));
                    else
                        cand = 1:SEQUENCE_LENGTH;
                        cand = cand(cand ~= LastLetter);
                        NewLetter = cand(randi(length(cand),1));
                        
                        cand = 1:SEQUENCE_LENGTH;
                        cand = cand(cand ~= LastLocation & cand ~= 1);
                        NewLocation = cand(randi(length(cand),1));
                    end
                end
            case 2
                %return to a sequence
                if strcmp(Content,'verbal')
                    NextLetter = DelayLetter + 1;
                    if NextLetter > SEQUENCE_LENGTH
                        NextLetter = 1;
                    end
                    NextLocation = 0;
                else
                    NextLocation = DelayLocation + 1;
                    if NextLocation > SEQUENCE_LENGTH
                        NextLocation = 1;
                    end
                    NextLetter = 0;
                end
                
                if strcmp(MatchType,'Match')
                    if strcmp(Content,'verbal')
                        %match verbal trial, take next letter in sequence
                        NewLetter = NextLetter;
                        
                        %get random new location that does not match last
                        %location
                        cand = 1:SEQUENCE_LENGTH;
                        cand = cand(cand ~= LastLocation);
                        NewLocation = cand(randi(length(cand),1));
                    else
                        %match spatial trial, take next locatoin in
                        %sequence
                        NewLocation = NextLocation;
                        
                        %get random new letter that does not match last
                        %letter
                        cand = 1:SEQUENCE_LENGTH;
                        cand = cand(cand ~= LastLetter);
                        NewLetter = cand(randi(length(cand),1));
                    end
                else
                    %NonMatch trial, get random letters and locations that
                    %do not match last trial or pre-delay
                    if strcmp(Content,'verbal')
                        cand = 1:SEQUENCE_LENGTH;
                        cand = cand(cand ~= LastLetter & cand ~= NextLetter & cand ~= DelayLetter);
                        NewLetter = cand(randi(length(cand),1));
                        
                        cand = 1:SEQUENCE_LENGTH;
                        cand = cand(cand ~= LastLocation);
                        NewLocation = cand(randi(length(cand),1));
                    else
                        cand = 1:SEQUENCE_LENGTH;
                        cand = cand(cand ~= LastLetter);
                        NewLetter = cand(randi(length(cand),1));
                        
                        cand = 1:SEQUENCE_LENGTH;
                        cand = cand(cand ~= LastLocation & cand ~= NextLocation & cand ~= DelayLocation);
                        NewLocation = cand(randi(length(cand),1));
                    end
                end
        end
    end

%% Trial Procedure
    function [trialData,TrialInfo,RunInfo,expData] = TrialProc(TrialInfo,BlockInfo,RunInfo,expData,IsPractice)
        %start with some fixation
        %this will give time to figure out what stimulus to present
                
        %Screen('DrawText',window,'+',Res.width/2 - FixXAdjust,Res.height/2 - FixYAdjust,[255 255 255]);
        %Screen('DrawingFinished',window);
        
        %% Determine Stimulus
        MatchType = TrialInfo.MatchType;
        Content = TrialInfo.Content;
        WM = TrialInfo.WM;
        DT = TrialInfo.DT;
        Phase = TrialInfo.Phase;
        PhaseSwitch = TrialInfo.PhaseSwitch;
        LastLetter = TrialInfo.LastLetter;
        LastLocation = TrialInfo.LastLocation;
        DelayLetter = TrialInfo.DelayLetter;
        DelayLocation = TrialInfo.DelayLocation;
        ITI = TrialInfo.ITI;
        
        if strcmp(MatchType,'Match')
            CorrectResponse = YesKey;
            IncorrectResponse = NoKey;
        else
            CorrectResponse = NoKey;
            IncorrectResponse = YesKey;
        end
        
        if strcmp(Content,'verbal')
            StimColor = VerbalColor;
        else
            StimColor = SpatialColor;
        end
        
        if Phase == 2
            switch [WM DT]
                case 'yesyes'
                    %branch
                    StimShape = BranchShape;
                case 'yesno'
                    %delay
                    StimShape = DelayShape;
                    CorrectResponse = NoKey;
                case 'noyes'
                    %restart
                    StimShape = SwitchShape;
                case 'nono'
                    %control
                    StimShape = 'Square';
            end
        else
            StimShape = 'Square';
        end
        
        switch num2str([Phase PhaseSwitch])
            case '1  1'
                %Init Block
                [LastLetter,LastLocation] = GetNewStimulus(MatchType,Content,1,LastLetter,LastLocation,0,0);
                TrialType = [Content 'Init'];
            case '1  0'
                %Baseline 1
                [LastLetter,LastLocation] = GetNewStimulus(MatchType,Content,0,LastLetter,LastLocation,0,0);
                TrialType = [Content 'Baseline1'];
            case '2  1'
                %Init Sub-task
                DelayLetter = LastLetter;
                DelayLocation = LastLocation;
                if isequal([WM DT],'nono')
                    %Control block, continue as normal
                    [LastLetter,LastLocation] = GetNewStimulus(MatchType,Content,0,LastLetter,LastLocation,0,0);
                else
                    [LastLetter,LastLocation] = GetNewStimulus(MatchType,Content,1,LastLetter,LastLocation,0,0);
                end
                TrialType = ['Init' BlockInfo.BlockType];
            case '2  0'
                %Sub-task
                [LastLetter,LastLocation] = GetNewStimulus(MatchType,Content,0,LastLetter,LastLocation,0,0);
                TrialType = BlockInfo.BlockType;
            case '3  1'
                %Return
                switch [WM DT]
                    case 'nono'
                        %Control, continue as normal
                        [LastLetter,LastLocation] = GetNewStimulus(MatchType,Content,0,LastLetter,LastLocation,0,0);
                    case 'noyes'
                        %Restart, restart sequence
                        [LastLetter,LastLocation] = GetNewStimulus(MatchType,Content,1,LastLetter,LastLocation,0,0);
                    case {'yesno','yesyes'}
                        %Return trial with delay
                        [LastLetter,LastLocation] = GetNewStimulus(MatchType,Content,2,LastLetter,LastLocation,DelayLetter,DelayLocation);
                end
                TrialType = ['Return' BlockInfo.BlockType];
            case '3  0'
                %Baseline 2
                [LastLetter,LastLocation] = GetNewStimulus(MatchType,Content,0,LastLetter,LastLocation,0,0);
                TrialType = [Content 'Baseline2'];
        end

        Letter = Letters{LastLetter};
        Location = Locations{LastLocation};
        
        TrialInfo.LastLetter = LastLetter;
        TrialInfo.LastLocation = LastLocation;
        TrialInfo.DelayLetter = DelayLetter;
        TrialInfo.DelayLocation = DelayLocation;
                
        c = find(strcmp(FrameColors,StimColor));
        s = find(strcmp(FrameShapes,StimShape));
        l = find(strcmp(Letters,Letter));
        
        Stimulus = TextureMat(c,s,l);
        
%% Present Stimuli
        %presenting starting fixation now
        %TrialOnset = Screen('Flip',window,RunInfo.MyClock-FlipDur/2);
        %disp(['Trial Onset Error: ' num2str(TrialOnset-RunInfo.MyClock)]);
        %RunInfo.MyClock = RunInfo.MyClock + STARTING_FIX;
        
        %Stimulus
        Screen('DrawTexture',window,Stimulus,[],Location);
        Screen('DrawText',window,'+',Res.width/2 - FixXAdjust,Res.height/2 - FixYAdjust,[255 255 255]);      
        StimOnset = Screen('Flip',window,RunInfo.MyClock-FlipDur/2);
        StimOnsetError = StimOnset - RunInfo.MyClock;
        RunInfo.MyClock = RunInfo.MyClock + STIM_DUR; 
        KbQueueFlush; %Flush buffer so only responses after stim onset are recorded     

        %ITI
        Screen('DrawText',window,'+',Res.width/2 - FixXAdjust,Res.height/2 - FixYAdjust,[255 255 255]);
        ITIOnset = Screen('Flip',window,RunInfo.MyClock-FlipDur/2);
        RunInfo.MyClock = RunInfo.MyClock + ITI;
        
        %Wait the ITI for a response, giving some buffer to perform other
        %code
        WaitSecs('UntilTime',RunInfo.MyClock - 0.3);

        %Get Response Data
        [pressed,firstPress] = KbQueueCheck; %check for a keypress
        response = '';
        thisRT = 0;
        if pressed
            if firstPress(escapeKey)
                %escape the program, saving what data we can

                %datafilename = sprintf('datafMRI_%d_%d_R%d',subNum, visit,num2str(BlockInfo.RunNum));
                %save(datafilename,'data');
                save(['dataBehav_' num2str(subNum) '_' num2str(BlockInfo.RunNum)],'data');

                KbQueueRelease;
                ListenChar(0);
                ShowCursor;
                error('Script aborted by escape command');
            end
                
            firstPress(firstPress==0)=NaN; %get rid of zeros
            [secs,idx] = min(firstPress); %get RT of first keypress
            disp('Getting response... ');
            response = KbName(idx);
            disp('Response: ');
            disp(response);
            thisRT = secs-StimOnset;
        end
        
        if IsPractice
            disp("inside IsPractice")
            if isequal(response,CorrectResponse)
                [bounds] = Screen('TextBounds',window,'Correct!!');
                Screen('DrawText',window,'Correct!!',Res.width/2 - bounds(3)/2,Res.height/2 - bounds(4)/2,[0 255 0]);
                Screen('Flip',window);
                RunInfo.MyClock = RunInfo.MyClock + FB_SLIDE;
        
                
            elseif isequal(response,IncorrectResponse)
                [bounds] = Screen('TextBounds',window,'Incorrect');
                Screen('DrawText',window,'Incorrect',Res.width/2 - bounds(3)/2,Res.height/2 - bounds(4)/2,[255 0 0]);
                Screen('Flip',window);
                RunInfo.MyClock = RunInfo.MyClock + FB_SLIDE;

            else
                [bounds] = Screen('TextBounds',window,'No input detected');
                Screen('DrawText',window,'No input detected',Res.width/2 - bounds(3)/2,Res.height/2 - bounds(4)/2,[255 0 0]);
                Screen('Flip',window);
                RunInfo.MyClock = RunInfo.MyClock + FB_SLIDE;
            end

            WaitSecs('UntilTime',RunInfo.MyClock - 0.3);
            Screen('DrawText',window,'+',Res.width/2 - FixXAdjust,Res.height/2 - FixYAdjust,[255 255 255]);
            Screen('Flip',window,RunInfo.MyClock-FlipDur/2);
            RunInfo.MyClock = RunInfo.MyClock + ITI/2;
            
        end
        
        if isempty(response)
            response= 'NA';
            thisACC= 0;
            thisRT= 0;
        elseif sum(response==CorrectResponse) >= 1
                thisACC = 1;
        else
                thisACC = 0;
        end
        
        %thisACC = (response==CorrectResponse);

        trialData.RESP = response;
        trialData.RT = thisRT;
        trialData.ACC = thisACC;
        trialData.Onset = StimOnset - RunInfo.RunStart;
        trialData.OnsetError = StimOnsetError;
        trialData.Duration = ITIOnset-StimOnset;
        trialData.DurationError = trialData.Duration - STIM_DUR;
        trialData.TrialType = TrialType;
        trialData.MatchType = TrialInfo.MatchType;
        trialData.Letter = Letter;
        trialData.Location = Location;
        trialData.Content = BlockInfo.Content;
        trialData.WM = BlockInfo.WM;
        trialData.DT = BlockInfo.DT;
        trialData.RunNum = BlockInfo.RunNum;
        trialData.BlockNum = BlockInfo.BlockNum;
        trialData.TrialNum = BlockInfo.TrialNum;
        
        if ~IsPractice
            %Update Score
            if trialData.ACC
                if trialData.RT < BASE_RT
                    expData.Score = expData.Score + 1;
                    if PhaseSwitch && Phase == 3
                        expData.Score = expData.Score + 2;
                    end
                end
            else
                expData.Score = expData.Score - 1;
            end

            %Update Onset File
            if isequal(trialData.RESP,LEFT_RESPONSE)
                fprintf(OnsetFile,'1%d\t%s\t %3.2f\t%2.2f\t%2.2f\t %2.2f\n', ...
                            trialData.RunNum,'LeftResponse',trialData.Onset,0, ...
                            1,trialData.RT);
            elseif isequal(trialData.RESP,RIGHT_RESPONSE)
                fprintf(OnsetFile,'%d\t%s\t%3.2f\t%2.2f\t%2.2f\t%2.2f\n', ...
                            trialData.RunNum,'RightResponse',trialData.Onset,0, ...
                            1,trialData.RT);
            end   

            switch TrialType
                case {'verbalInit' 'spatialInit'}
                    %write init transient
                    fprintf(OnsetFile,'%d\t%s\t%3.2f\t%2.2f\t%2.2f\t%2.2f\n', ...
                            trialData.RunNum,TrialType,trialData.Onset,0, ...
                            trialData.ACC,trialData.RT);

                case {'InitVerbalControl','InitVerbalDelay','InitVerbalSwitch','InitVerbalBranch', ...
                      'InitSpatialControl','InitSpatialDelay','InitSpatialSwitch','InitSpatialBranch'}
                    %write out baseline epoch
                    %write init transient
                    if (expData.EpochDur-expData.LastITI) < 1
                        fprintf(OnsetFile,'%d\t%s\t%3.2f\t%2.2f\t%2.2f\t%2.2f\n', ...
                            trialData.RunNum,expData.EpochType,expData.EpochOns,0, ...
                            expData.EpochCorr/expData.EpochNum,1);
                    else
                        fprintf(OnsetFile,'%d\t%s\t%3.2f\t%2.2f\t%2.2f\t%2.2f\n', ...
                            trialData.RunNum,expData.EpochType,expData.EpochOns,expData.EpochDur-expData.LastITI, ...
                            expData.EpochCorr/expData.EpochNum,1);
                    end
                    fprintf(OnsetFile,'%d\t%s\t%3.2f\t%2.2f\t%2.2f\t%2.2f\n', ...
                            trialData.RunNum,TrialType,trialData.Onset,0, ...
                            trialData.ACC,trialData.RT);

                    expData.EpochOns = trialData.Onset;
                    expData.EpochType = TrialType;
                    expData.EpochDur = 0;
                    expData.EpochNum = 0;
                    expData.EpochCorr = 0;

                case {'ReturnVerbalControl','ReturnVerbalDelay','ReturnVerbalSwitch','ReturnVerbalBranch', ...
                      'ReturnSpatialControl','ReturnSpatialDelay','ReturnSpatialSwitch','ReturnSpatialBranch'}
                    %write out sub-task epoch
                    %write return transient
                    fprintf(OnsetFile,'%d\t%s\t%3.2f\t%2.2f\t%2.2f\t%2.2f\n', ...
                            trialData.RunNum,expData.EpochType,expData.EpochOns,expData.EpochDur-expData.LastITI, ...
                            expData.EpochCorr/expData.EpochNum,1);
                    fprintf(OnsetFile,'%d\t%s\t%3.2f\t%2.2f\t%2.2f\t%2.2f\n', ...
                            trialData.RunNum,TrialType,trialData.Onset,0, ...
                            trialData.ACC,trialData.RT);

                    expData.EpochDur = 0;
                    expData.EpochNum = 0;
                    expData.EpochCorr = 0;
                    expData.EpochType = '';

                case {'VerbalControl','VerbalDelay','VerbalSwitch','VerbalBranch','verbalBaseline1','verbalBaseline2', ...
                      'SpatialControl','SpatialDelay','SpatialSwitch','SpatialBranch','spatialBaseline1','spatialBaseline2'}
                    if expData.EpochDur == 0
                        %starting a new epoch
                        expData.EpochOns = trialData.Onset;
                        expData.EpochType = TrialType;
                    end

                    expData.EpochDur = expData.EpochDur + STIM_DUR + ITI;
                    expData.EpochNum = expData.EpochNum + 1;
                    expData.EpochCorr = expData.EpochCorr + trialData.ACC;

                    %record errors for separate co-variate
                    if trialData.ACC == 0
                        fprintf(OnsetFile,'%d\t%s\t%3.2f\t%2.2f\t%2.2f\t%2.2f\n', ...
                            trialData.RunNum,'Error',trialData.Onset,0, ...
                            trialData.ACC,trialData.RT);
                    end
            end

            expData.LastITI = ITI;
        end
                        
    end %Trial

%% Block
    function [runData,BlockInfo,RunInfo,expData] = BlockProcedure(runData,BlockInfo,RunInfo,expData,IsPractice)
        trialStart = BlockInfo.TrialNum;
        TrialNum = BlockInfo.TrialNum;
        numTrials = BlockInfo.numTrials;
        
        TrialInfo.BlockType = BlockInfo.BlockType;
        TrialInfo.Content = BlockInfo.Content;
        TrialInfo.WM = BlockInfo.WM;
        TrialInfo.DT = BlockInfo.DT;
                
        TrialInfo.LastLetter = 0;
        TrialInfo.LastLocation = 0;
        TrialInfo.DelayLetter = 0;
        TrialInfo.DelayLocation = 0;
              
        %phase 1 - task
        TrialInfo.Phase = 1;
        TrialInfo.PhaseSwitch = 1;
        
        for trials1 = 1:numTrials(1)
            TrialInfo.TrialNum = TrialNum;
            TrialInfo.ITI = RunInfo.TrialJitters(TrialNum);
            TrialInfo.MatchType = RunInfo.MatchList{TrialNum};
            [trialData,TrialInfo,RunInfo,expData] = TrialProc(TrialInfo,BlockInfo,RunInfo,expData,IsPractice);
            runData(TrialNum) = trialData;
            TrialNum = TrialNum + 1;
            TrialInfo.PhaseSwitch = 0;
        end
        
        %phase 2 - sub-task
        TrialInfo.Phase = 2;
        TrialInfo.PhaseSwitch = 1;
        
        for trials2 = 1:numTrials(2)
            TrialInfo.TrialNum = TrialNum;
            TrialInfo.ITI = RunInfo.TrialJitters(TrialNum);
            TrialInfo.MatchType = RunInfo.MatchList{TrialNum};
            [trialData,TrialInfo,RunInfo,expData] = TrialProc(TrialInfo,BlockInfo,RunInfo,expData,IsPractice);
            runData(TrialNum) = trialData;
            TrialNum = TrialNum + 1;
            TrialInfo.PhaseSwitch = 0;
        end
        
        %phase 3 - task
        TrialInfo.Phase = 3;
        TrialInfo.PhaseSwitch = 1;
        
        for trials3 = 1:numTrials(3)
            TrialInfo.TrialNum = TrialNum;
            TrialInfo.ITI = RunInfo.TrialJitters(TrialNum);
            TrialInfo.MatchType = RunInfo.MatchList{TrialNum};
            [trialData,TrialInfo,RunInfo,expData] = TrialProc(TrialInfo,BlockInfo,RunInfo,expData,IsPractice);
            runData(TrialNum) = trialData;
            TrialNum = TrialNum + 1;
            TrialInfo.PhaseSwitch = 0;
        end
        trialEnd = TrialNum-1;
        
        %Finish the ITI
        %Screen('DrawText',window,'+',Res.width/2 - FixXAdjust,Res.height/2 - FixYAdjust,[255 255 255]);
        %Screen('Flip',window,RunInfo.MyClock-FlipDur/2);
        %RunInfo.MyClock = RunInfo.MyClock + STARTING_FIX;
        
        %Feedback
        numCorrect = sum([runData(trialStart:trialEnd).ACC]);
        numTrials = trialEnd-trialStart+1;
        
        if ~IsPractice
            %Modify Onset File
            %write out baseline epoch
            %write out feedback transient
            if (expData.EpochDur-expData.LastITI) < 1
                fprintf(OnsetFile,'%d\t%s\t%3.2f\t%2.2f\t%2.2f\t%2.2f\n', ...
                    trialData.RunNum,expData.EpochType,expData.EpochOns,0, ...
                    expData.EpochCorr/expData.EpochNum,1);
            else
                fprintf(OnsetFile,'%d\t%s\t%3.2f\t%2.2f\t%2.2f\t%2.2f\n', ...
                    trialData.RunNum,expData.EpochType,expData.EpochOns,expData.EpochDur-expData.LastITI, ...
                    expData.EpochCorr/expData.EpochNum,1);
            end
            fprintf(OnsetFile,'%d\t%s\t%3.2f\t%2.2f\t%2.2f\t%2.2f\n', ...
                    trialData.RunNum,'Feedback',RunInfo.MyClock-RunInfo.RunStart,0,1,1);

            expData.EpochOns = 0;
            expData.EpochType = '';
            expData.EpochDur = 0;
            expData.EpochNum = 0;
            expData.EpochCorr = 0;
        end
        
        str = sprintf('%d/%d Correct',numCorrect,numTrials);
        [bounds] = Screen('TextBounds',window,str);
        Screen('DrawText',window,str,Res.width/2 - bounds(3)/2,Res.height/2 - bounds(4)/2,[255 255 255]);
        Screen('Flip',window,RunInfo.MyClock-FlipDur/2);
        RunInfo.MyClock = RunInfo.MyClock + STIM_DUR;
        
        %IBI
        Screen('DrawText',window,'+',Res.width/2 - FixXAdjust,Res.height/2 - FixYAdjust,[255 255 255]);
        Screen('Flip',window,RunInfo.MyClock-FlipDur/2);
        RunInfo.MyClock = RunInfo.MyClock+BlockInfo.IBI;
        
        BlockInfo.TrialNum = TrialNum;
        
    end %Block

%% Practice
    function PracticeProc(scoreData)
        tic;
        Screen('TextSize',window,34);
        disp('inside PracticeProc');
        disp(YesKey);
        disp(NoKey);
     
        
        Instr = sprintf(['Welcome back! You will now receive practice on all of the possible tasks while\n', ...
                        'we take a high-resolution image of your brain. It is VERY important to keep\n', ...
                        'your head still during this time.\n\n Remember the rules:\n\n'...
                        '1) The letter task is indicated by %s colored frames and the\n'...
                        'location task is indicated by %s colored frames\n\n'...
                        '2) The Basic task is indicated by Square frames\n\n'...
                        '3) The Delay task is indicated by %s frames\n\n'...
                        '4) The Restart task is indicated by %s frames\n\n'...
                        '5) The Dual task is indicated by %s frames\n\n'...
                        'To respond Yes press the %s button and to respond No press the %s button\n\n'...
                        '(Press B to begin practice)'],VerbalColor, SpatialColor,DelayShape,SwitchShape,BranchShape,YesKey,NoKey);

        DrawFormattedText(window,Instr,'center','center',[255 255 255]);
        Screen('Flip',window);
        KbWait([],2);
   
        
        Screen('TextSize',window,40);
        
        BlockIdx = 1;
        TrialIdx = 1:(TRIALS_PER_RUN/2);
        NUM_PRACTICE_BLOCKS = 8;
        
        practiceData = data(TrialIdx);

        BlockInfo.TrialNum = 1;
        RunInfo.TrialJitters = Shuffle(TrialJitters);
        RunInfo.BlockJitters = Shuffle(BlockJitters);
        RunInfo.MatchList = Shuffle(MatchList);
        
        %Create Lists for Practice
        BlockTypePracticeList = {'VerbalControl' 'VerbalDelay' 'VerbalSwitch' 'VerbalBranch', ...
                                 'SpatialControl' 'SpatialDelay' 'SpatialSwitch' 'SpatialBranch'};
        ContentPracticeList = {'verbal' 'verbal' 'verbal' 'verbal' 'spatial' 'spatial' 'spatial' 'spatial'};
        WMPracticeList = {'no' 'yes' 'no' 'yes' 'no' 'yes' 'no' 'yes'};
        DTPracticeList = {'no' 'no' 'yes' 'yes' 'no' 'no' 'yes' 'yes'};
        
        randIdx = randperm(8);
        BlockTypePracticeList = BlockTypePracticeList(randIdx);
        ContentPracticeList = ContentPracticeList(randIdx);
        WMPracticeList = WMPracticeList(randIdx);
        DTPracticeList = DTPracticeList(randIdx);
        
        numTrialsPracticeList = [randi(2,NUM_PRACTICE_BLOCKS,1)+1 randi(3,NUM_PRACTICE_BLOCKS,1)+2 randi(2,NUM_PRACTICE_BLOCKS,1)+1];
        

        %Wait for scanner ready
        Instr = 'Please wait and remain still while the scan is prepared';
        [bounds] = Screen('TextBounds',window,Instr);
        Screen('DrawText',window,Instr,Res.width/2 - bounds(3)/2,Res.height/2 - bounds(4)/2,[255 255 255]);
        Screen('Flip',window);
        KbWait([],2);

        [bounds] = Screen('TextBounds',window,'Get Ready!!!');
        Screen('DrawText',window,'Get Ready!!!',Res.width/2 - bounds(3)/2,Res.height/2 - bounds(4)/2,[255 255 255]);
        Screen('DrawingFinished', window);

        %Wait for trigger
        RunInfo.RunStart = KbWait;
        RunInfo.MyClock = RunInfo.RunStart;
        Screen('Flip',window);
        RunInfo.MyClock = RunInfo.MyClock + 3; %3 second ready cue

        Screen('DrawText',window,'+',Res.width/2 - FixXAdjust,Res.height/2 - FixYAdjust,[255 255 255]);
        Screen('Flip',window,RunInfo.MyClock-(FlipDur/2));   
        RunInfo.MyClock = RunInfo.MyClock + 1; %1 second post-cue fixation

        for blocki = 1:NUM_PRACTICE_BLOCKS
            %Choose Block Variables
            BlockInfo.BlockType = BlockTypePracticeList{BlockIdx};
            BlockInfo.numTrials = numTrialsPracticeList(BlockIdx,:);
            BlockInfo.Content = ContentPracticeList{BlockIdx};
            BlockInfo.WM = WMPracticeList{BlockIdx};
            BlockInfo.DT = DTPracticeList{BlockIdx};
            BlockInfo.IBI = RunInfo.BlockJitters(blocki);
            BlockInfo.RunNum = 1;
            BlockInfo.BlockNum = blocki;

            [practiceData,BlockInfo,RunInfo,scoreData] = BlockProcedure(practiceData,BlockInfo,RunInfo,scoreData,1);

            BlockIdx = BlockIdx + 1;
        end
        
        save(['practiceData_' num2str(subNum) '_' num2str(visit)],'practiceData');
        
        Instr = sprintf(['Practice has ended. The experiment will begin soon. Please remember to maintain\n', ...
                        'your gaze at the center fixation at all times and perform your best.\n', ...
                        'Remember to keep your head still!!!\n\n']);
        DrawFormattedText(window,Instr,'center','center',[255 255 255]);
        Screen('Flip',window);
        KbWait([],2); 
        
        %clear out score data
        scoreData.Score = 0;
        scoreData.CurrACC = 0;
        scoreData.CurrRT = 0;
        scoreData.LastACC = 0;
        scoreData.LastRT = 0;
        scoreData.TotalACC = 0;
        scoreData.TotalRT = 0;
        scoreData.TotalTrials = 0;
%         KbQueueStop;
%         KbQueueRelease;
        elapsedTime = toc;
        disp(['Elapsed practice time: ', num2str(elapsedTime), ' seconds']);
    end

%% Instructions
    function Instructions
        Screen('TextSize',window,32);
        Instr = sprintf('Welcome to the behavioral experiment!\n\n(Press any key to continue)');
        DrawFormattedText(window,Instr,'center','center',[255 255 255]);
        Screen('Flip',window);
        KbWait([],2);
        
        Instr = sprintf(['We will begin by reviewing the tasks you will be performing for this experiment.\n', ...
                        'You will then receive a short practice session to get you familiar with the tasks.\n', ...
                        'Afterwards, the experiment will begin.\n\n', ...
                        '(Press any key to continue)']);
        DrawFormattedText(window,Instr,'center','center',[255 255 255]);
        Screen('Flip',window);
        KbWait([],2);
        
        Instr = sprintf(['In this experiment, you will be asked to make decisions\n', ...
                        'based on letters and spatial locations.\n\n', ...
                        '(Press any key to continue)']);
        DrawFormattedText(window,Instr,'center','center',[255 255 255]);
        Screen('Flip',window);
        KbWait([],2);
        
        Instr = sprintf(['On each trial, you will be presented with a letter at a particular location on\n', ...
                        'the screen, surrounded by a colored frame.\n\n', ...
                        'On some trials, your decision will be based upon the letter.\n', ...
                        'On other trials, your decision will be based upon the location.\n\n', ...
                        'The color of the frame will indicate which aspect of the stimulus\n', ...
                        '(either letter or location) is important for the trial.\n\n', ...
                        '(Press any key to continue)']);
        DrawFormattedText(window,Instr,'center','center',[255 255 255]);
        Screen('Flip',window);
        KbWait([],2);
        
        Screen('DrawTexture',window,TextureMat(1,4,1),[],Locations{1});
        Screen('DrawText',window,'+',Res.width/2 - FixXAdjust,Res.height/2 - FixYAdjust,[255 255 255]);
        DrawFormattedText(window,'This is an example stimulus','center',200,[255 255 255]);
        DrawFormattedText(window,'(Press any key to continue)','center',Res.height-200,[255 255 255]);
        Screen('Flip',window);
        KbWait([],2);
        
        Instr = sprintf(['If the color of the frame is %s, you should pay attention to the letter\n', ...
                         'If the color of the frame is %s, you should pay attention to the location\n\n', ...
                         '(Press any key to continue)'],VerbalColor,SpatialColor);
        DrawFormattedText(window,Instr,'center','center',[255 255 255]);
        Screen('Flip',window);
        KbWait([],2);
        
        Instr = sprintf(['Please use the index fingers of each hand to make the responses.\n',...
                         'To make yes responses, press the %s button\n', ...
                         'To make no responses, press the %s button\n\n', ...
                         '(Press any key to continue)'],YesKey,NoKey);
        DrawFormattedText(window,Instr,'center','center',[255 255 255]);
        Screen('Flip',window);
        KbWait([],2);
        
%Check Button Input
        Instr = sprintf(['Let''s test the buttons.\n\n', ...
                         'During the experiment, please keep one index finger on the %s\n\n', ...
                         '(Press %s to continue)'],YesKey,YesKey);
        DrawFormattedText(window,Instr,'center','center',[255 255 255]);
        Screen('Flip',window);
        RestrictKeysForKbCheck(KbName(YesKey));
        KbWait([],2);
        
        Instr = sprintf(['Good! Now, please keep your other index finger on the %s\n\n', ...
                         '(Press %s button to continue)'],NoKey,NoKey);
        DrawFormattedText(window,Instr,'center','center',[255 255 255]);
        Screen('Flip',window);
        RestrictKeysForKbCheck(KbName(NoKey));
        KbWait([],2);
        RestrictKeysForKbCheck([]);
        
%Basic Task Instructions        
        Instr = sprintf(['Your basic task is to determine whether the current stimulus\n', ...
                         'follows the previous stimulus in a given sequence. For the\n', ...
                         'very first stimulus, you will determine whether\n',...
                         'it is the first of the sequence.\n\n',...
                         '(Press any key to continue)']);
        DrawFormattedText(window,Instr,'center','center',[255 255 255]);
        Screen('Flip',window);
        KbWait([],2);
        
%Letter Task Instructions
        Instr = sprintf(['For the letter task, the sequence follows the word "TABLET".\n', ...
                        'For each letter, you will determine whether it follows the previous\n', ...
                        'letter in the word "TABLET".\n\n', ...
                        'For example, if the stimulus is "B" and the previous stimulus was "A",\n', ...
                        'you would make a "yes" response.\n\n', ...
                        'At the start of every block, you will determine whether the first\n',...
                        'letter is "T"(i.e. the first letter of "TABLET"). Afterwards,\n',...
                        'you will perform the task as described above.\n\n', ...
                        '(Press any key to continue)']);
        DrawFormattedText(window,Instr,'center','center',[255 255 255]);
        Screen('Flip',window);
        KbWait([],2);

%Letter Task: Example
        instructColor = find(strcmp(FrameColors,VerbalColor));
        Screen('DrawTexture',window,TextureMat(instructColor,4,2),[],Locations{1});
        Screen('DrawText',window,'+',Res.width/2 - FixXAdjust,Res.height/2 - FixYAdjust,[255 255 255]);
        Instr = sprintf(['Let''s do some example trials of the letter task:\n',...
                        'What would your response be if this was the very first stimulus of a block?']);
        DrawFormattedText(window,Instr,'center',200,[255 255 255]);
        DrawFormattedText(window,'(respond to stimulus to continue)','center',Res.height-200,[255 255 255]);
        Screen('Flip',window);

        [~,keyCode,~] = KbWait([],2);
            if KbName(keyCode)== NoKey
                Instr = sprintf(['Correct! Because this is the first stimulus of the block\n',...
                                'you are looking for a "T". Let''s try another.\n\n',...
                                '(Press any key to continue)']);
                DrawFormattedText(window,Instr,'center','center',[255 255 255]);
                Screen('Flip',window);
            else 
                Instr = sprintf(['Incorrect. Because this is the first stimulus of the block you are\n',...
                                'looking for a "T" (the first letter of "TABLET"). Let''s try another.\n\n',...
                                '(Press any key to continue)']);
                DrawFormattedText(window,Instr,'center','center',[255 255 255]);
                Screen('Flip',window);
            end
        KbWait([],2);

        Screen('DrawTexture',window,TextureMat(instructColor,4,1),[],Locations{5});
        Screen('DrawText',window,'+',Res.width/2 - FixXAdjust,Res.height/2 - FixYAdjust,[255 255 255])
        Instr = sprintf('What would your response be if this was the very first stimulus of the block?');
        DrawFormattedText(window,Instr,'center',200,[255 255 255]);
        DrawFormattedText(window,'(respond to stimulus to continue)','center',Res.height-200,[255 255 255]);
        Screen('Flip',window);

        [~,keyCode,~] = KbWait([],2);
            if KbName(keyCode)== YesKey
                Instr = sprintf(['Correct! Because this is the first stimulus of the block\n',...
                                'you are looking for a "T".\n\n',...
                                '(Press any key to continue)']);
                DrawFormattedText(window,Instr,'center','center',[255 255 255]);
                Screen('Flip',window);
            else 
                Instr = sprintf(['Incorrect. Because this is the first stimulus of the block\n',...
                                'you are looking for a "T".\n\n',...
                                '(Press any key to continue)']);
                DrawFormattedText(window,Instr,'center','center',[255 255 255]);
                Screen('Flip',window);
            end
        KbWait([],2);

%Letter Task Practice Start Screen
        Instr = sprintf(['Now let''s try a whole block (multiple trials).\n\n',...
                        'Remember, after every stimulus you must give a yes or no response.\n',...
                        'After the first stimulus, the question you should be asking is:\n',...
                        'Does this letter follow the previous letter, in the word "TABLET"?\n\n',...
                        '(Press any key to begin practice)']);
        DrawFormattedText(window,Instr,'center', 'center',[255 255 255]);
        Screen('Flip',window);
        KbWait([],2);

%Practice Basic Task: Letter  
        keyCode = 'r';
        while keyCode == 'r'
                instructPractice({'VerbalControl'},{'verbal'},{'no'},{'no'},1,1)

                Screen('TextSize',window,32);

                Instr = sprintf(['How did you do? Don''t worry if you did poorly.\n',...
                                'If you would like to try again hit the ''r'' key, otherwise\n',...
                                'press any other key to continue with the instructions.']);
                DrawFormattedText(window,Instr,'center', 'center',[255 255 255]);
                Screen('Flip',window);
                KbWait([],2);

                keyDown = 0;
                while ~keyDown
                    [keyDown,~,keyPress] = KbCheck;
                    keyCode = KbName(keyPress); 
                end
        end 

%Location Instructions
        Instr = sprintf(['Now let''s go over the location task.\n\n',...
                        '(Press any key to continue)']);
        DrawFormattedText(window,Instr,'center', 'center',[255 255 255]);
        Screen('Flip',window);
        KbWait([],2);
        
        KbQueueFlush;
        LCount = 0; 
        keyDown = 0;
        while ~keyDown
            ii = mod(LCount,5) + 1; 
           
            Instr = sprintf('For the location task, the sequence of locations are as follows:');
            
            DrawFormattedText(window,Instr,'center',200,[255 255 255]);
            instructColor = find(strcmp(FrameColors,SpatialColor));
            Screen('DrawTexture',window,TextureMat(instructColor,4,1),[],Locations{ii});
            Screen('DrawText',window,'+',Res.width/2 - FixXAdjust,Res.height/2 - FixYAdjust,[255 255 255]);
           
            Instr = sprintf(['You can think of the locations as points on a star.\n', ...
                             'Please take a moment to memorize this sequence.\n', ...
                             'At the start of the series, you will determine whether the location\n', ...
                             'is the top point. Afterwards, you will perform the location task\n', ...
                             'just as described for the letter task.\n\n',...
                             '(Hold down on any key to continue)']);
            DrawFormattedText(window,Instr,'center',Res.height-300,[255 255 255]);
            Screen('Flip',window);

            LCount = LCount + 1;
            WaitSecs(1); 
            %[keyDown] = KbCheck;

            % Check for a key press
            [keyDown, ~, keyPress] = KbCheck; % Get the key state

            % Optional: Debugging output to see which key was pressed
            if keyDown
                pressedKey = KbName(find(keyPress)); % Get the name of the pressed key
                disp(['Key pressed: ', pressedKey]); % Display the key pressed
            end
                    
        end
        
%Location Task Practice Start Screen
        Instr = sprintf(['Again, just like the letter task, after the first stimulus,\n'...
                        'the question you should be asking is:\n',...
                        'Does this position follow the previous position in the location sequence?\n\n',...
                        'Let''s jump right in and try a whole block of the location task.\n\n',...
                        '(Press any key to begin practice)']);
        DrawFormattedText(window,Instr,'center','center',[255 255 255]);
        Screen('Flip',window);
        KbWait([],2);
        
%Practice Basic Task: Location
        keyCode = 'r';
        while keyCode == 'r'
                instructPractice({'SpatialControl'},{'spatial'},{'no'},{'no'},1,1)

                Screen('TextSize',window,32);

                Instr = sprintf(['How did you do? Don''t worry if you did poorly.\n',...
                                'If you would like to try again hit the ''r'' key, otherwise\n',...
                                'press any other key to continue with the instructions.']);
                DrawFormattedText(window,Instr,'center', 'center',[255 255 255]);
                Screen('Flip',window);
                KbWait([],2);

                keyDown = 0;
                while ~keyDown
                    [keyDown,~,keyPress] = KbCheck;
                    keyCode = KbName(keyPress); 
                end
        end
        
%Feedback and Fixation Rules
        Instr = sprintf(['After each block, you will receive feedback\n', ...
                         'regarding your performance (as you''ve seen).\n\n', ...
                         'Feedback will be in the form of: # correct/ # trials in block\n\n', ...
                         'After feedback, a new block will automatically start.\n\n', ...
                         '(Press any key to continue)']);
        DrawFormattedText(window,Instr,'center','center',[255 255 255]);
        Screen('Flip',window);
        KbWait([],2);
        
        Instr = sprintf(['At all times, we ask that you maintain your gaze at the center fixation\n', ...
                         'cross ("+"). So, please do NOT move your eyes to look at the stimuli.\n', ...
                         'Do your best to process the stimuli using your peripheral vision.\n', ...
                         '(Press any key to continue)']);
        DrawFormattedText(window,Instr,'center','center',[255 255 255]);
        Screen('Flip',window);
        KbWait([],2);
        
%Changing Task Instructions
        Instr = sprintf(['Each block of trials will begin with either the basic letter or location task.\n', ...
                         'After a few trials, you may be required to perform one of three additional tasks.\n', ...
                         'These additional tasks will be cued by a change in the shape of the frame\n', ...
                         'surrounding the letter. The frame may change from a square to either:\n', ...
                         'a circle, diamond, or cross.\n\n', ...
                         '(Press any key to continue)']);
        DrawFormattedText(window,Instr,'center','center',[255 255 255]);
        Screen('Flip',window);
        KbWait([],2);
                     
        Instr = sprintf(['For example, the square frames indicate the basic task for a few trials.\n', ...
                         'Then the frames may change shape (e.g. circle) indicating a new task for a few\n', ...
                         'trials. Finally, they will return to squares. All blocks will begin and end\n', ...
                         'with square frames(the basic task).\n\n', ...
                         '(Press any key to continue)']);
        DrawFormattedText(window,Instr,'center','center',[255 255 255]);
        Screen('Flip',window);
        KbWait([],2);

%Delay Task Instructions
        Instr = sprintf(['One possible additional task is the "Delay" task. When the frame changes to a %s\n', ...
                         'you must:\n\n', ...
                         '1) Remember where you were in the sequence before the frame changed to a %s.\n\n', ...
                         '2) While the frame remains a %s, respond "no" to all stimuli regardless of\n', ...
                         'their sequence.\n\n', ...
                         '3) When the frame returns to the square shape, resume the sequence from where you\n', ...
                         'left off in step 1.\n\n', ...
                         '(Press any key to continue)'],DelayShape,DelayShape,DelayShape);
        DrawFormattedText(window,Instr,'center','center',[255 255 255]);
        Screen('Flip',window);
        KbWait([],2);

%Delay Practice Start Screen
        Instr = sprintf(['Let''s practice doing a "Delay" and put everything you''ve learned together.\n\n',...
                         'For this round, you will be shown TWO blocks. One block will be the\n',...
                         'location task (%s frames) and the other will be the letter task (%s frames).\n\n',...
                         '(Press any key to begin practice)'],SpatialColor,VerbalColor);
        DrawFormattedText(window,Instr,'center','center',[255 255 255]);
        Screen('Flip',window);
        KbWait([],2);
         
%Delay Practice: 2 Blocks
        keyCode = 'r';
        while keyCode == 'r'
                instructPractice({'SpatialDelay' 'VerbalDelay'},{'spatial' 'verbal'},{'yes' 'yes'},{'no' 'no'},1,2)

                Screen('TextSize',window,32);

                Instr = sprintf(['How did you do? Don''t worry if you did poorly.\n',...
                                'If you would like to try again hit the ''r'' key, otherwise\n',...
                                'press any other key to continue with the instructions.']);
                DrawFormattedText(window,Instr,'center', 'center',[255 255 255]);
                Screen('Flip',window);
                KbWait([],2);

                keyDown = 0;
                while ~keyDown
                    [keyDown,~,keyPress] = KbCheck;
                    keyCode = KbName(keyPress); 
                end
        end 
          
%Restart Task Instructions
        Instr = sprintf(['Another possible additional task is the "Restart" task. When the frame changes to a %s\n', ...
                         'you must:\n\n', ...
                         '1)Begin the task again. If you were on the letter task and the frame\n',...
                         'changed shape, determine if the letter inside the circle is a "T".\n'...
                         'For the location task, determine if the location is the top point.\n\n'...
                         '2) While the frame remains the same, continue the task as normal.\n\n', ...
                         '3) If the frame changes back to squares, "RESTART" again (back to step 1).\n\n', ...
                         '(Press any key to continue)'],SwitchShape);
        DrawFormattedText(window,Instr,'center','center',[255 255 255]);
        Screen('Flip',window);
        KbWait([],2);
        
%Restart Practice Start Screen
        Instr = sprintf(['Let''s practice doing a "Restart" and put everything you''ve learned together.\n\n',...
                         'For this round, you will be shown two blocks. One block will be the\n',...
                         'location task (%s frames) and the other will be the letter task (%s frames).\n\n',...
                         '(Press any key to begin practice)'],SpatialColor,VerbalColor);
        DrawFormattedText(window,Instr,'center','center',[255 255 255]);
        Screen('Flip',window);
        KbWait([],2);
         
%Restart Practice: 2 Blocks
        keyCode = 'r';
        while keyCode == 'r'
                instructPractice({'SpatialSwitch' 'VerbalSwitch'},{'spatial' 'verbal'},{'no' 'no'},{'yes' 'yes'},1,2)

                Screen('TextSize',window,32);

                Instr = sprintf(['How did you do? Don''t worry if you did poorly.\n',...
                                'If you would like to try again hit the ''r'' key, otherwise\n',...
                                'press any other key to continue with the instructions.']);
                DrawFormattedText(window,Instr,'center', 'center',[255 255 255]);
                Screen('Flip',window);
                KbWait([],2);

                keyDown = 0;
                while ~keyDown
                    [keyDown,~,keyPress] = KbCheck;
                    keyCode = KbName(keyPress); 
                end
        end 

%Dual Task Instructions
        Instr = sprintf(['The last possible additional task is the "Dual" task. When the frame changes to a %s\n', ...
                         'you must:\n\n', ...
                         '1) Remember where you were in the sequence before the frame changed to a %s\n', ...
                         'AND restart the basic task.\n\n', ...
                         '2) While the frame remains a %s, continue the basic task while continuing to\n', ...
                         'remember where you were before the frame changed.\n\n', ...
                         '3) When the frame returns to the square shape, resume the sequence from where you\n', ...
                         'left off in step 1.\n\n', ...
                         '(Press any key to continue)'],BranchShape,BranchShape,BranchShape);
        DrawFormattedText(window,Instr,'center','center',[255 255 255]);
        Screen('Flip',window);
        KbWait([],2);

%Dual Practice Start Screen
        Instr = sprintf(['Let''s practice doing a "Dual" and put everything you''ve learned together.\n\n',...
                         'For this round, you will be shown two blocks. One block will be the\n',...
                         'location task (%s frames) and the other will be the letter task (%s frames).\n\n',...
                         '(Press any key to begin practice)'],SpatialColor,VerbalColor);
        DrawFormattedText(window,Instr,'center','center',[255 255 255]);
        Screen('Flip',window);
        KbWait([],2);
         
%Dual Practice: 2 Blocks
        keyCode = 'r';
        while keyCode == 'r'
                instructPractice({'SpatialBranch' 'VerbalBranch'},{'spatial' 'verbal'},{'yes' 'yes'},{'yes' 'yes'},1,2)

                Screen('TextSize',window,32);

                Instr = sprintf(['How did you do? Don''t worry if you did poorly.\n',...
                                'If you would like to try again hit the ''r'' key, otherwise\n',...
                                'press any other key to continue with the instructions.']);
                DrawFormattedText(window,Instr,'center', 'center',[255 255 255]);
                Screen('Flip',window);
                KbWait([],2);

                keyDown = 0;
                while ~keyDown
                    [keyDown,~,keyPress] = KbCheck;
                    keyCode = KbName(keyPress); 
                end
        end 
        
%End Instructions
        Instr = sprintf(['Note that if the frames do not change shape on a given block,\n',...
                        'just continue the task as normal for the entire block.\n\n', ...
                        '(Press any key to continue)']);
        DrawFormattedText(window,Instr,'center','center',[255 255 255]);
        Screen('Flip',window);
        KbWait([],2);

        Instr = sprintf(['This concludes the instructions.\n', ...
                         'Please ask the experimenter if anything is unclear.\n\n', ...
                         '(Press any key to continue)']);
        DrawFormattedText(window,Instr,'center','center',[255 255 255]);
        Screen('Flip',window);
        KbWait([],2);
        
        Screen('TextSize',window,40);
    end
%% Rest Proc
    function RestProc(scoreData)
        Screen('TextSize',window,32);
        
        Screen('FillRect',window,[0 76 151]);
        Instr = sprintf(['Please take a moment to rest.\n\n', ...
                        'Current Score: %d\n\n', ...
                        'This Run Accuracy: %2.2f%%\n', ...
                        'This Run Average Speed: %d ms\n\n', ...
                        'Last Run Accuracy: %2.2f%%\n', ...
                        'Last Run Average Speed: %d ms\n\n', ...
                        'Total Accuracy: %2.2f%%\n', ...
                        'Total Speed: %d ms\n\n'], ...
                        scoreData.Score, ...
                        (scoreData.CurrACC/TRIALS_PER_RUN)*100, ...
                        round((scoreData.CurrRT/TRIALS_PER_RUN)*1000), ...
                        (scoreData.LastACC/TRIALS_PER_RUN)*100, ...
                        round((scoreData.LastRT/TRIALS_PER_RUN)*1000), ...
                        (scoreData.TotalACC/scoreData.TotalTrials)*100, ...
                        round((scoreData.TotalRT/scoreData.TotalTrials)*1000));
                    
        DrawFormattedText(window,Instr,'center','center',[255 255 255]);
        Screen('Flip',window);
        WaitSecs(10); %enforce 10 second rest
        
        Instr = sprintf(['Please take a moment to rest.\n\n', ...
                        'Current Score: %d\n\n', ...
                        'This Run Accuracy: %2.2f%%\n', ...
                        'This Run Average Speed: %d ms\n\n', ...
                        'Last Run Accuracy: %2.2f%%\n', ...
                        'Last Run Average Speed: %d ms\n\n', ...
                        'Total Accuracy: %2.2f%%\n', ...
                        'Total Speed: %d ms\n\n', ...
                        '(Press any key to continue)'], ...
                        scoreData.Score, ...
                        (scoreData.CurrACC/TRIALS_PER_RUN)*100, ...
                        round((scoreData.CurrRT/TRIALS_PER_RUN)*1000), ...
                        (scoreData.LastACC/TRIALS_PER_RUN)*100, ...
                        round((scoreData.LastRT/TRIALS_PER_RUN)*1000), ...
                        (scoreData.TotalACC/scoreData.TotalTrials)*100, ...
                        round((scoreData.TotalRT/scoreData.TotalTrials)*1000));
                    
        DrawFormattedText(window,Instr,'center','center',[255 255 255]);
        Screen('Flip',window);
        
        KbWait([],2);
        Screen('FillRect',window,[0 0 0]);
        Screen('TextSize',window,40);
    end         

%% Instructional Practice
function instructPractice(BlockTypePracticeList,ContentPracticeList,WMPracticeList,DTPracticeList,BlockIdx,NUM_PRACTICE_BLOCKS)
        %Input Variables: BlockTypePracticeList,ContentPracticeList...
        %WMPracticeList,DTPracticeList,PreSwitchPracticeList...
        %BlockIdx, and NUM_PRACTICE_BLOCKS

        Screen('TextSize',window,40);
        
        BlockInfo.TrialNum = 1;
        TrialIdx = 1:(TRIALS_PER_RUN/2);
        practiceData = data(TrialIdx);
        numTrialsPracticeList = [randi(2,NUM_PRACTICE_BLOCKS,1)+1 randi(3,NUM_PRACTICE_BLOCKS,1)+2 randi(2,NUM_PRACTICE_BLOCKS,1)+1];
              
        RunInfo.MatchList = Shuffle(MatchList);
        RunInfo.TrialJitters = Shuffle(TrialJitters);
        RunInfo.BlockJitters = Shuffle(BlockJitters);

        %Wait for scanner ready
        [bounds] = Screen('TextBounds',window,'Get Ready!!!');
        Screen('DrawText',window,'Get Ready!!!',Res.width/2 - bounds(3)/2,Res.height/2 - bounds(4)/2,[255 255 255]);
        Screen('DrawingFinished', window);

        %Wait for trigger
        RunInfo.RunStart = KbWait;
        RunInfo.MyClock = RunInfo.RunStart;
        Screen('Flip',window);
        RunInfo.MyClock = RunInfo.MyClock + 3; %3 second ready cue

        Screen('DrawText',window,'+',Res.width/2 - FixXAdjust,Res.height/2 - FixYAdjust,[255 255 255]);
        Screen('Flip',window,RunInfo.MyClock-(FlipDur/2));   
        RunInfo.MyClock = RunInfo.MyClock + 1; %1 second post-cue fixation

        for blocki = 1:NUM_PRACTICE_BLOCKS
        %Choose Block Variables
            BlockInfo.BlockType = BlockTypePracticeList{BlockIdx};
            BlockInfo.numTrials = numTrialsPracticeList(BlockIdx,:);
            BlockInfo.Content = ContentPracticeList{BlockIdx};
            BlockInfo.WM = WMPracticeList{BlockIdx};
            BlockInfo.DT = DTPracticeList{BlockIdx};
            BlockInfo.IBI = RunInfo.BlockJitters(blocki);
            BlockInfo.ITI = RunInfo.TrialJitters(blocki);
            BlockInfo.RunNum = 1;
            BlockInfo.BlockNum = blocki;

            [practiceData,BlockInfo,RunInfo,scoreData] = BlockProcedure(practiceData,BlockInfo,RunInfo,scoreData,1);

            BlockIdx = BlockIdx + 1;
        end
    %clear out score data
        scoreData.Score = 0;
        scoreData.CurrACC = 0;
        scoreData.CurrRT = 0;
        scoreData.LastACC = 0;
        scoreData.LastRT = 0;
        scoreData.TotalACC = 0;
        scoreData.TotalRT = 0;
        scoreData.TotalTrials = 0;
end         

%% Termination
%These are now handled by onCleanup
%ListenChar(0);
%ShowCursor;
%sca
end