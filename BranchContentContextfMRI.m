function [data] = BranchContentContextfMRI(subNum, visit)

if nargin < 1
    subNum = input('Subject number:');
    visit = input('Visit number:');
end

%% Experiment variables
%RUN_DURATION = 574; %length of a run in seconds
RUN_DURATION = 583; %length of a run in seconds, pre-task time was added for dis acqs, DEN 1/17/25
TRIALS_PER_RUN = 144; %number of trials per run, divisible by 3
BLOCKS_PER_RUN = 16; %number of blocks per run, divisible by 2
NUM_RUNS = 6; %number of runs in the session
BASE_RT = 1; %RT to beat to get bonus
LEFT_RESPONSE = '4$'; %left response key
RIGHT_RESPONSE = '1!'; %right response key
SEQUENCE_LENGTH = 5; %number of stimuli in a sequence
STIM_DUR = 0.5; %Stimulus duration
FB_SLIDE = 0.3; %Feedback time (correct/incorrect)
%STARTING_FIX = 0.5; %Fixation time to start each trial for timing purposes

%variables for tool problem testing
KbName('UnifyKeyNames');
any_key_cont= [KbName('1') KbName('1!') KbName('2') KbName('2@') KbName('3') KbName('3#') KbName('4') KbName('4$') KbName('ESCAPE')];
select_run = [KbName('y') KbName('1') KbName('1!') KbName('2') KbName('2@') KbName('3') KbName('3#') KbName('4') KbName('4$') KbName('5') KbName('5%') KbName('6') KbName('6^')];
menu_keys = [KbName('i') KbName('e') KbName('q') KbName('p')];
ttl_inputs= [KbName('5') KbName('5%') KbName('6') KbName('6^')];
button_keys= zeros(1,256);
button_keys(any_key_cont)= 1;

%% Counter-balancing
%color-to-content mappings
switch mod(subNum,4)
    case {1,2}
        VerbalColor = 'Blue';
        SpatialColor = 'Yellow';
    case {3,0}
        VerbalColor = 'Yellow';
        SpatialColor = 'Blue';
end

%shape-to-task mappings
switch mod(subNum,24)
    case {1,2,3,4}
        DelayShape = 'Circle';
        SwitchShape = 'Diamond';
        BranchShape = 'Cross';
    case {5,6,7,8}
        DelayShape = 'Cross';
        SwitchShape = 'Circle';
        BranchShape = 'Diamond';
    case {9,10,11,12}
        DelayShape = 'Diamond';
        SwitchShape = 'Cross';
        BranchShape = 'Circle';
    case {13,14,15,16}
        DelayShape = 'Circle';
        SwitchShape = 'Cross';
        BranchShape = 'Diamond';
    case {17,18,19,20}
        DelayShape = 'Cross';
        SwitchShape = 'Diamond';
        BranchShape = 'Circle';
    case {21,22,23,0}
        DelayShape = 'Diamond';
        SwitchShape = 'Circle';
        BranchShape = 'Cross';
end

%decision-response mappings
if mod(subNum,2)
    YesKey = LEFT_RESPONSE;
    NoKey = RIGHT_RESPONSE;
else
    YesKey = RIGHT_RESPONSE;  
    NoKey = LEFT_RESPONSE;
end


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
ListNum = mod(subNum,50) + 1;
filename = ['\Users\aniladmello\Downloads\untitled folder\MATLAB\Lists\StimulusList_ContextTMS' num2str(ListNum) '.txt'];
%filename = ['.\Lists\ListsfMRI1Session\StimulusListRandOneSession' num2str(ListNum) '.txt'];

fid = fopen(filename);
d = textscan(fid,'%d\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\r','HeaderLines',1);
fclose(fid);

ContentList = d{3};
WMList = d{4};
DTList = d{5};
BlockTypeList = d{6};
numTrialsList = [d{7} d{8} d{9}];

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
Res = Screen('Resolution',0);
window = Screen('OpenWindow',0,[0 0 0]);

closeWindow = onCleanup(@() Screen('CloseAll'));

Priority(MaxPriority(0));

restorePriority = onCleanup(@() Priority(0));

FlipDur = Screen('GetFlipInterval', window);

%code for changing resolution and centering
%oldRes = Screen('Resolution',0);
%newRes.width = 1024; newRes.height = 1024;
%windowRect = CenterRect([0 0 newRes.width newRes.height],[0 0 oldRes.width oldRes.height]);
%window = Screen('OpenWindow',0,[0 0 0],windowRect);

%% Load Stimuli as Textures

FrameColors = {'Blue' 'Yellow'};
FrameShapes = {'Circle' 'Cross' 'Diamond' 'Square'};
Letters = {'T' 'A' 'B' 'L' 'E'};
StimDir = 'C:\Users\Help Desk\Documents\My Experiments\Alexa\ContextfMRI\Stimuli\';
%StimDir = '.\Stimuli\';
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

%Create Onset File
onsetfilename = sprintf('OnsetTimesfMRI_%d_%d.txt',subNum,visit);
OnsetFile = fopen(onset ,'a');
closeOnsetFile = onCleanup(@() fclose(OnsetFile));

while ~QuitExp
    MenuStr = sprintf('Main Menu\n\n(I)nstructions\n(P)ractice\n(E)xperiment\n(Q)uit');
    %MenuStr = sprintf('Main Menu\n\n(P)ractice\n(E)xperiment\n(Q)uit');
    DrawFormattedText(window,MenuStr,'center','center',[255 255 255]);
    Screen('Flip',window);
    
    RestrictKeysForKbCheck(menu_keys);
    [~,keycode] = KbWait([],2);
    RestrictKeysForKbCheck([]);

    switch KbName(keycode)
        case 'i'
            Instructions;
        case 'p'
            PracticeProc(scoreData);
        case 'e'
            RunStr = sprintf('Hit (y) to continue to run %d or enter run number',RunSuggestion);
            [bounds] = Screen('TextBounds',window,RunStr);
            Screen('DrawText',window,RunStr,Res.width/2 - bounds(3)/2,Res.height/2 - bounds(4)/2,[255 255 255]);
            Screen('Flip',window);
            
            RestrictKeysForKbCheck(select_run);
            [~,keycode] = KbWait([],2);
            RestrictKeysForKbCheck([]);
            
            keyname = KbName(keycode);
            if keyname=='y'
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
                        RunSuggestion = RunNum; %update run to what user input, DEN, 1/17/25
                    end
                end
            end
            
            RestProc(scoreData);
            
            RunSuggestion = RunSuggestion + 1;
            if RunSuggestion > NUM_RUNS
                break;
            end
        case 'q'
            break;
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
        Instr = 'Please wait and remain still while the next scan is prepared';
        [bounds] = Screen('TextBounds',window,Instr);
        Screen('DrawText',window,Instr,Res.width/2 - bounds(3)/2,Res.height/2 - bounds(4)/2,[255 255 255]);
        Screen('Flip',window);

        %KbTriggerWait (Wait for trigger)
        KbQueueStop;
        KbQueueRelease;
        %RunInfo.RunStart= KbTriggerWait(ttl_inputs);
        KbQueueCreate(-1);
        KbQueueStart;

        %changing from KbTriggerWait because TTLs are being missed
        %occasionally, DEN, 1/17/25

        pressed = false;
        while ~pressed
            [pressed,firstPress] = KbQueueCheck; %check for a keypress  
            firstPress(firstPress==0)=NaN; %get rid of zeros
            [secs,idx] = min(firstPress); %get time of first keypress
            %response = KbName(idx);
            if ~any(ismember(idx,ttl_inputs))
                pressed = false;
            end
        end
        RunInfo.RunStart = secs;

        %start a new queue that only accepts the button box keys to avoid
        %TTLs being recorded as responses
        KbQueueStop;
        KbQueueRelease;
        KbQueueCreate(-1,button_keys);
        KbQueueStart;

        RunInfo.MyClock = RunInfo.RunStart;
        Screen('Flip',window);
        RunInfo.MyClock = RunInfo.MyClock + 4; %First 4 second fixation (for dummy scan)

        Screen('DrawText',window,'+',Res.width/2 - FixXAdjust,Res.height/2 - FixYAdjust,[255 255 255]);
        Screen('Flip',window,RunInfo.MyClock-(FlipDur/2));   
        RunInfo.MyClock = RunInfo.MyClock + 4; %Second 4 second fixation (for dummy scan)
                
        [bounds] = Screen('TextBounds',window,'Get Ready!!!');
        Screen('DrawText',window,'Get Ready!!!',Res.width/2 - bounds(3)/2,Res.height/2 - bounds(4)/2,[255 255 255]);
        Screen('Flip',window,RunInfo.MyClock-(FlipDur/2));
        RunInfo.MyClock = RunInfo.MyClock + 2; %2 second pause

        Screen('DrawText',window,'+',Res.width/2 - FixXAdjust,Res.height/2 - FixYAdjust,[255 255 255]);
        Screen('Flip',window,RunInfo.MyClock-(FlipDur/2));   
        RunInfo.MyClock = RunInfo.MyClock + 3; %3 second post-cue fixation

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
        save(['datafMRI_' num2str(subNum) '_' num2str(visit)],'data'); %save at end of each run
        
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
                save(['datafMRI_' num2str(subNum) '_' num2str(visit) '_' num2str(BlockInfo.RunNum)],'data');
                KbQueueRelease;
                ListenChar(0);
                ShowCursor;
                error('Script aborted by escape command');
            end
                
            firstPress(firstPress==0)=NaN; %get rid of zeros
            [secs,idx] = min(firstPress); %get RT of first keypress
            response = KbName(idx);
            thisRT = secs-StimOnset;
        end
        
        if IsPractice
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
                fprintf(OnsetFile,'%d\t%s\t%3.3f\t%2.3f\t%2.3f\t%2.3f\n', ...
                            trialData.RunNum,'LeftResponse',trialData.Onset,0, ...
                            1,trialData.RT);
            elseif isequal(trialData.RESP,RIGHT_RESPONSE)
                fprintf(OnsetFile,'%d\t%s\t%3.3f\t%2.3f\t%2.3f\t%2.3f\n', ...
                            trialData.RunNum,'RightResponse',trialData.Onset,0, ...
                            1,trialData.RT);
            end   

            switch TrialType
                case {'verbalInit' 'spatialInit'}
                    %write init transient
                    fprintf(OnsetFile,'%d\t%s\t%3.3f\t%2.3f\t%2.3f\t%2.3f\n', ...
                            trialData.RunNum,TrialType,trialData.Onset,0, ...
                            trialData.ACC,trialData.RT);

                case {'InitVerbalControl','InitVerbalDelay','InitVerbalSwitch','InitVerbalBranch', ...
                      'InitSpatialControl','InitSpatialDelay','InitSpatialSwitch','InitSpatialBranch'}
                    %write out baseline epoch
                    %write init transient
                    if (expData.EpochDur-expData.LastITI) < 1
                        fprintf(OnsetFile,'%d\t%s\t%3.3f\t%2.3f\t%2.3f\t%2.3f\n', ...
                            trialData.RunNum,expData.EpochType,expData.EpochOns,0, ...
                            expData.EpochCorr/expData.EpochNum,1);
                    else
                        fprintf(OnsetFile,'%d\t%s\t%3.3f\t%2.3f\t%2.3f\t%2.3f\n', ...
                            trialData.RunNum,expData.EpochType,expData.EpochOns,expData.EpochDur-expData.LastITI, ...
                            expData.EpochCorr/expData.EpochNum,1);
                    end
                    fprintf(OnsetFile,'%d\t%s\t%3.3f\t%2.3f\t%2.3f\t%2.3f\n', ...
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
                    fprintf(OnsetFile,'%d\t%s\t%3.3f\t%2.3f\t%2.3f\t%2.3f\n', ...
                            trialData.RunNum,expData.EpochType,expData.EpochOns,expData.EpochDur-expData.LastITI, ...
                            expData.EpochCorr/expData.EpochNum,1);
                    fprintf(OnsetFile,'%d\t%s\t%3.3f\t%2.3f\t%2.3f\t%2.3f\n', ...
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
                        fprintf(OnsetFile,'%d\t%s\t%3.3f\t%2.3f\t%2.3f\t%2.3f\n', ...
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
                fprintf(OnsetFile,'%d\t%s\t%3.3f\t%2.3f\t%2.3f\t%2.3f\n', ...
                    trialData.RunNum,expData.EpochType,expData.EpochOns,0, ...
                    expData.EpochCorr/expData.EpochNum,1);
            else
                fprintf(OnsetFile,'%d\t%s\t%3.3f\t%2.3f\t%2.3f\t%2.3f\n', ...
                    trialData.RunNum,expData.EpochType,expData.EpochOns,expData.EpochDur-expData.LastITI, ...
                    expData.EpochCorr/expData.EpochNum,1);
            end
            fprintf(OnsetFile,'%d\t%s\t%3.3f\t%2.3f\t%2.3f\t%2.3f\n', ...
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
        Screen('TextSize',window,32);
        
        if isequal(YesKey,LEFT_RESPONSE)
            YesInstruct = 'Left-Most';
            NoInstruct = 'Right-Most';
        else
            YesInstruct = 'Right-Most';
            NoInstruct = 'Left-Most';
        end
        
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
                        '(Press any key to begin practice)'],VerbalColor, SpatialColor,DelayShape,SwitchShape,BranchShape,YesInstruct,NoInstruct);
        DrawFormattedText(window,Instr,'center','center',[255 255 255]);
        Screen('Flip',window);
        %KbWait([],2);
        %KbQueueStop;
        %KbQueueRelease;
        %RunInfo.RunStart= KbTriggerWait(any_key_cont);
        %KbQueueCreate(-1);
        %KbQueueCreate(-1,button_keys);
        %KbQueueStart;
        KbQueueWait([],2); %simplifying this, DEN, 1/17/25
        
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
        %KbWait([],2);

        [bounds] = Screen('TextBounds',window,'Get Ready!!!');
        Screen('DrawText',window,'Get Ready!!!',Res.width/2 - bounds(3)/2,Res.height/2 - bounds(4)/2,[255 255 255]);
        Screen('DrawingFinished', window);

        %Wait for trigger
        %KbQueueStop;
        %KbQueueRelease;
        %RunInfo.RunStart= KbTriggerWait(ttl_inputs);
        %KbQueueCreate(-1,button_keys);
        %KbQueueStart;
        
        %changing from KbTriggerWait because TTLs are being missed
        %occasionally, DEN, 1/17/25

        pressed = false;
        while ~pressed
            [pressed,firstPress] = KbQueueCheck; %check for a keypress  
            firstPress(firstPress==0)=NaN; %get rid of zeros
            [secs,idx] = min(firstPress); %get time of first keypress
            %response = KbName(idx);
            if ~any(ismember(idx,ttl_inputs))
                pressed = false;
            end
        end
        RunInfo.RunStart = secs;

        %start a new queue that only accepts the button box keys to avoid
        %TTLs being recorded as responses
        KbQueueStop;
        KbQueueRelease;
        KbQueueCreate(-1,button_keys);
        KbQueueStart;

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
        
        %Not clear why this is here, removing, DEN, 1/17/25
        %KbQueueStop;
        %KbQueueRelease;
        %RunInfo.RunStart= KbTriggerWait(ttl_inputs);
        %KbQueueCreate(-1);
        %KbQueueCreate(-1,button_keys);
        %KbQueueStart; 
        
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

    end

%% Instructions
    function Instructions
        Screen('TextSize',window,32);
        Instr = sprintf('Welcome to the experiment!\n\n(Press any key to continue)');
        DrawFormattedText(window,Instr,'center','center',[255 255 255]);
        Screen('Flip',window);
        KbQueueWait([],2);
        
        Instr = sprintf(['We will begin by reviewing the instructions you received earlier\n', ...
                        'in order to refresh you with the task.\n\n', ...
                        'You will then receive a short practice session to get you familiar\n', ...
                        'with performing the task in the scanner environment.\n\n', ...
                        'Afterwards, the task will begin.\n\n', ...
                        '(Press any key to continue)']);
        DrawFormattedText(window,Instr,'center','center',[255 255 255]);
        Screen('Flip',window);
        KbQueueWait([],2);
        
        Instr = sprintf(['In this experiment, you will need to memorize and make decisions\n', ...
                        'upon letters and spatial locations.\n\n', ...
                        '(Press any key to continue)']);
        DrawFormattedText(window,Instr,'center','center',[255 255 255]);
        Screen('Flip',window);
        KbQueueWait([],2);
        
        Instr = sprintf(['On each trial, you will be presented with a letter at a particular location on\n', ...
                        'the screen. The letter will be surrounded by a colored frame.\n\n', ...
                        'On some trials, your decision will be based upon the letter.\n', ...
                        'On some trials, your decision will be based upon the location.\n\n', ...
                        'The color of the frame will indicate which aspect of the stimulus\n', ...
                        '(either letter or location) is important for the trial.\n\n', ...
                        '(Press any key to continue)']);
        DrawFormattedText(window,Instr,'center','center',[255 255 255]);
        Screen('Flip',window);
        KbQueueWait([],2);
        
        Screen('DrawTexture',window,TextureMat(1,4,1),[],Locations{1});
        Screen('DrawText',window,'+',Res.width/2 - FixXAdjust,Res.height/2 - FixYAdjust,[255 255 255]);
        DrawFormattedText(window,'This is an example stimulus','center',200,[255 255 255]);
        DrawFormattedText(window,'(Press any key to continue)','center',Res.height-200,[255 255 255]);
        Screen('Flip',window);
        KbQueueWait([],2);
        
        Instr = sprintf(['If the color of the frame is %s, you should pay attention to the letter\n', ...
                         'If the color of the frame is %s, you should pay attention to the location\n\n', ...
                         '(Press any key to continue)'],VerbalColor,SpatialColor);
        DrawFormattedText(window,Instr,'center','center',[255 255 255]);
        Screen('Flip',window);
        KbQueueWait([],2);
        
        Instr = sprintf(['Your basic task will be to determine whether the presented stimulus\n', ...
                         'follows the previous stimulus in a sequence.\n\n', ...
                        'For the letter task, the sequence consists of the word "TABLET".\n', ...
                        'For each letter, you will determine whether it follows the previous\n', ...
                        'letter in the word "TABLET".\n\n', ...
                        'For example, if the stimulus is "B" and the previous stimulus was "A",\n', ...
                        'you would make a "yes" response.\n\n', ...
                        'At the start of the sequence, you will determine whether the first letter is "T"\n', ...
                        '(i.e. the first letter of "TABLET"). Afterwards, you will perform the sequence\n', ...
                        'letter-matching task described above.\n\n', ...
                        '(Press any key to continue)']);
        DrawFormattedText(window,Instr,'center','center',[255 255 255]);
        Screen('Flip',window);
        KbQueueWait([],2);
        
        Instr = sprintf('For the location task, the sequence is:');
        DrawFormattedText(window,Instr,'center',200,[255 255 255]);
        SpatialSequence = Screen('MakeTexture',window,imread('SpatialSequence.png'));
        %Image is 400 x 400
        Screen('DrawTexture',window,SpatialSequence,[],[Res.width/2-200 Res.height/2-200 Res.width/2+200 Res.height/2+200]);
        Instr = sprintf(['You can think of the locations as points on a star.\n', ...
                         'Please take a moment to make sure you remember this sequence.\n', ...
                         'At the start of the sequence, you will determine whether the location\n', ...
                         'is the top point. Afterwards, you will perform the sequence location-matching\n', ...
                         'task described above.\n\n', ...
                         '(Press any key to continue)']);
        DrawFormattedText(window,Instr,'center',Res.height-300,[255 255 255]);
        Screen('Flip',window);
        KbQueueWait([],2);
        
        Instr = sprintf(['Sequences will be grouped into a block of trials. After each block, you\n', ...
                         'will receive feedback regarding your performance. This feedback will be\n', ...
                         'in the form of:\n\n', ...
                         '# correct/# trials in block\n\n', ...
                         'After feedback, a new sequence will start.\n\n', ...
                         '(Press any key to continue)']);
        DrawFormattedText(window,Instr,'center','center',[255 255 255]);
        Screen('Flip',window);
        KbQueueWait([],2);
        
        if isequal(YesKey,LEFT_RESPONSE)
            YesInstruct = 'Left-Most';
            NoInstruct = 'Right-Most';
        else
            YesInstruct = 'Right-Most';
            NoInstruct = 'Left-Most';
        end
        
        Instr = sprintf(['To make a yes response, press the %s button\n', ...
                         'To make a no response, press the %s button\n', ...
                         'Please use the index fingers of each hand to make the responses.\n\n', ...
                         '(Press any key to continue)'],YesInstruct,NoInstruct);
        DrawFormattedText(window,Instr,'center','center',[255 255 255]);
        Screen('Flip',window);
        KbQueueWait([],2);
        
        Instr = sprintf(['Let''s test the buttons.\n', ...
                         'Please respond with your left index finger on the left-most button\n\n', ...
                         '(Press left-most button to continue)']);
        DrawFormattedText(window,Instr,'center','center',[255 255 255]);
        Screen('Flip',window);
        %RestrictKeysForKbCheck(KbName(LEFT_RESPONSE));
        %KbQueueWait([],2);
        pressed = false;
        while ~pressed
            [pressed,firstPress] = KbQueueCheck; %check for a keypress  
            firstPress(firstPress==0)=NaN; %get rid of zeros
            [secs,idx] = min(firstPress); %get RT of first keypress
            response = KbName(idx);
            if ~isequal(response,LEFT_RESPONSE)
                pressed = false;
            end
        end


        Instr = sprintf(['Good! Now, please respond with your right index finger on the right-most button.\n\n', ...
                         '(Press right-most button to continue)']);
        DrawFormattedText(window,Instr,'center','center',[255 255 255]);
        Screen('Flip',window);
        %RestrictKeysForKbCheck(KbName(RIGHT_RESPONSE));
        %KbQueueWait([],2);
        %RestrictKeysForKbCheck([]);
        pressed = false;
        while ~pressed
            [pressed,firstPress] = KbQueueCheck; %check for a keypress  
            firstPress(firstPress==0)=NaN; %get rid of zeros
            [secs,idx] = min(firstPress); %get RT of first keypress
            response = KbName(idx);
            if ~isequal(response,RIGHT_RESPONSE)
                pressed = false;
            end
        end
        
        
        Instr = sprintf(['At all times, we ask that you maintain your gaze at the center fixation\n', ...
                         'cross ("+"). So, please do NOT move your eyes to look at the stimuli.\n', ...
                         'Do your best to process the stimuli using your peripheral attention.\n', ...
                         'We will be monitoring your eye movements for adherence to these instructions.\n\n', ...
                         '(Press any key to continue)']);
        DrawFormattedText(window,Instr,'center','center',[255 255 255]);
        Screen('Flip',window);
        KbQueueWait([],2);
        
        Instr = sprintf(['Each block of trials will begin with either the basic letter or location task.\n', ...
                         'After a few trials, you may be required to perform one of three additional tasks.\n', ...
                         'These additional tasks will be cued by a change in the shape of the frame\n', ...
                         'surrounding the letter. The frame may change from a square to either\n', ...
                         'a circle, diamond, or cross.\n\n', ...
                         'For example, the frames will be squares indicating the basic task for a few\n', ...
                         'trials. Then they may change shape (e.g. circle) indicating a new task for a few\n', ...
                         'trials. Finally, they will return to squares. All blocks will begin and end\n', ...
                         'with the basic task (square frames).\n\n', ...
                         '(Press any key to continue)']);
        DrawFormattedText(window,Instr,'center','center',[255 255 255]);
        Screen('Flip',window);
        KbQueueWait([],2);
        
        Instr = sprintf(['One possible task will be the "Delay" task. When the frame changes to a %s\n', ...
                         'you must:\n\n', ...
                         '1) Remember where you were in the sequence before the frame changed to a %s.\n\n', ...
                         '2) While the frame remains a %s, respond "no" to all stimuli regardless of\n', ...
                         'their sequence.\n\n', ...
                         '3) When the frame returns to the square shape, resume the sequence from where you\n', ...
                         'left off in 1).\n\n', ...
                         '(Press any key to continue)'],DelayShape,DelayShape,DelayShape);
        DrawFormattedText(window,Instr,'center','center',[255 255 255]);
        Screen('Flip',window);
        KbQueueWait([],2);
        
        Instr = sprintf(['Another possible task will be the "Restart" task. When the frame changes to a %s\n', ...
                         'you must:\n\n', ...
                         '1) Restart the basic task. If you were on the letter task, start by determining\n',...
                         'whether the stimulus is a "T". If you were on the location task, start by\n',...
                         'determing whether the location is the top point.\n\n', ...
                         '2) While the frame remains a %s, continue the basic task.\n\n', ...
                         '3) When the frame returns to the square shape, restart the basic task again.\n\n', ...
                         '(Press any key to continue)'],SwitchShape,SwitchShape);
        DrawFormattedText(window,Instr,'center','center',[255 255 255]);
        Screen('Flip',window);
        KbQueueWait([],2);
        
        Instr = sprintf(['The last possible task will be the "Dual" task. When the frame changes to a %s\n', ...
                         'you must:\n\n', ...
                         '1) Remember where you were in the sequence before the frame changed to a %s\n', ...
                         'AND restart the basic task.\n\n', ...
                         '2) While the frame remains a %s, continue the basic task while continuing to\n', ...
                         'remember where you were before the frame changed.\n\n', ...
                         '3) When the frame returns to the square shape, resume the sequence from where you\n', ...
                         'left off in 1).\n\n', ...
                         '(Press any key to continue)'],BranchShape,BranchShape,BranchShape);
        DrawFormattedText(window,Instr,'center','center',[255 255 255]);
        Screen('Flip',window);
        KbQueueWait([],2);
        
        Instr = sprintf(['Note that if the frames do not change shape on a given block, just continue the\n', ...
                         'basic task the entire block.\n\n', ...
                         '(Press any key to continue)']);
        DrawFormattedText(window,Instr,'center','center',[255 255 255]);
        Screen('Flip',window);
        KbQueueWait([],2);
        
        Instr = sprintf(['It is important that you perform at your best at all times. To facilitate this,\n', ...
                         'you will be rewarded for fast and accurate responses. Incorrect and slow\n', ...
                         'responses will incur a penalty. Since we realize that juggling multiple tasks is\n', ...
                         'especially hard, rewards will be increased for the trial following a switch back\n', ...
                         'to the basic task (e.g. the trial when frames change to squares after having been\n', ...
                         'some other shape).\n\n', ...
                         'During rest breaks, you will receive feedback regarding your speed, accuracy,\n', ...
                         'and current score. The score will be translated into a bonus payment at the end of\n', ...
                         'the experiment. So, perform your best to receive the greatest reward.\n\n', ...
                         '(Press any key to continue)']);
        DrawFormattedText(window,Instr,'center','center',[255 255 255]);
        Screen('Flip',window);
        KbQueueWait([],2);
        
        Instr = sprintf(['We will be collecting high-resolution images of your brain while you perform the\n', ...
                         'task. It is EXTREMELEY important that you keep your head as still as possible\n', ...
                         'during this time. Even 1 mm of movement can distort the images. Note that when\n', ...
                         'you move your arms and legs, your head moves too! So, please keep as still as\n', ...
                         'possible.\n\n', ...
                         '(Press any key to continue)']);
        DrawFormattedText(window,Instr,'center','center',[255 255 255]);
        Screen('Flip',window);
        KbQueueWait([],2);
        
        Instr = sprintf(['This concludes the instructions.\n', ...
                         'Please ask the experimenter if anything is unclear.\n\n', ...
                         '(Press any key to continue)']);
        DrawFormattedText(window,Instr,'center','center',[255 255 255]);
        Screen('Flip',window);
        KbQueueWait([],2);
        
        Screen('TextSize',window,40);
    end
            
%% Rest Proc
    function RestProc(scoreData)
        Screen('TextSize',window,32);
        
        %make the feedback screen FSU colors
        Screen('FillRect',window,[120 47 64]);
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
                    
        DrawFormattedText(window,Instr,'center','center',[206 184 136]);
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
                    
        DrawFormattedText(window,Instr,'center','center',[206 184 136]);
        Screen('Flip',window);
        
        KbQueueWait([],2);
        Screen('FillRect',window,[0 0 0]);
        Screen('TextSize',window,40);
    end

%% Instructional Practice
% function instructPractice(BlockTypePracticeList,ContentPracticeList,WMPracticeList,DTPracticeList,BlockIdx,NUM_PRACTICE_BLOCKS)
%         %Input Variables: BlockTypePracticeList,ContentPracticeList...
%         %WMPracticeList,DTPracticeList,PreSwitchPracticeList...
%         %BlockIdx, and NUM_PRACTICE_BLOCKS
% 
%         Screen('TextSize',window,40);
%         
%         BlockInfo.TrialNum = 1;
%         TrialIdx = 1:(TRIALS_PER_RUN/2);
%         practiceData = data(TrialIdx);
%         numTrialsPracticeList = [randi(2,NUM_PRACTICE_BLOCKS,1)+1 randi(3,NUM_PRACTICE_BLOCKS,1)+2 randi(2,NUM_PRACTICE_BLOCKS,1)+1];
%         
% %         trialCount = sum(numTrialsPracticeList,'all');
% %         matchFirstHalf = floor(trialCount/2);
% %         [MatchList{1 : matchFirstHalf}] = deal('Match');
% %         [MatchList{matchFirstHalf + 1 : trialCount}] = deal('NonMatch');
% %         RunInfo.MatchList = repmat(Shuffle(MatchList),1,NUM_PRACTICE_BLOCKS+1);
% %        
%         RunInfo.MatchList = Shuffle(MatchList);
%         RunInfo.TrialJitters = Shuffle(TrialJitters);
%         RunInfo.BlockJitters = Shuffle(BlockJitters);
% 
%         %Wait for scanner ready
%         [bounds] = Screen('TextBounds',window,'Get Ready!!!');
%         Screen('DrawText',window,'Get Ready!!!',Res.width/2 - bounds(3)/2,Res.height/2 - bounds(4)/2,[255 255 255]);
%         Screen('DrawingFinished', window);
% 
%         %Wait for trigger
%         RunInfo.RunStart = KbWait;
%         RunInfo.MyClock = RunInfo.RunStart;
%         Screen('Flip',window);
%         RunInfo.MyClock = RunInfo.MyClock + 3; %3 second ready cue
% 
%         Screen('DrawText',window,'+',Res.width/2 - FixXAdjust,Res.height/2 - FixYAdjust,[255 255 255]);
%         Screen('Flip',window,RunInfo.MyClock-(FlipDur/2));   
%         RunInfo.MyClock = RunInfo.MyClock + 1; %1 second post-cue fixation
% 
%         for blocki = 1:NUM_PRACTICE_BLOCKS
%         %Choose Block Variables
%             BlockInfo.BlockType = BlockTypePracticeList{BlockIdx};
%             BlockInfo.numTrials = numTrialsPracticeList(BlockIdx,:);
%             BlockInfo.Content = ContentPracticeList{BlockIdx};
%             BlockInfo.WM = WMPracticeList{BlockIdx};
%             BlockInfo.DT = DTPracticeList{BlockIdx};
%             %BlockInfo.preSwitches = PreSwitchPracticeList(BlockIdx);
%             BlockInfo.IBI = RunInfo.BlockJitters(blocki);
%             BlockInfo.ITI = RunInfo.TrialJitters(blocki);
%             BlockInfo.RunNum = 1;
%             BlockInfo.BlockNum = blocki;
% 
%             [practiceData,BlockInfo,RunInfo,scoreData] = BlockProcedure(practiceData,BlockInfo,RunInfo,scoreData,1);
% 
%             BlockIdx = BlockIdx + 1;
%         end
%     %clear out score data
%         scoreData.Score = 0;
%         scoreData.CurrACC = 0;
%         scoreData.CurrRT = 0;
%         scoreData.LastACC = 0;
%         scoreData.LastRT = 0;
%         scoreData.TotalACC = 0;
%         scoreData.TotalRT = 0;
%         scoreData.TotalTrials = 0;
% end         
%% Termination
%These are now handled by onCleanup
%ListenChar(0);
%ShowCursor;
%sca
end