%% Beginning
close all; clear all;clc;
% data_directory= 'C:\Users\Yawen.Wang\Documents\Matalb Scripts\NBM_fMRI_Behavior\data';
subjno=input('Please input date:', 's');
subjnum = '11';
% subjnum = input('Please input subject number:', 's');
subjname=input('Please input subject name:', 's');
filename = sprintf('%s_%s_%s',subjnum,subjname,subjno);
ET = 1;
if ET    
edffile = [upper(subjname) '.edf'];
data_directory=['D:\Yawen\VBehavior\data\NBM_fMRI_Behavior_data\Staircase']; % set a new directory

% eyedata_directory=['D:\Yawen\VBehavior\data\NBM_fMRI_Behavior_data\Eyedata'];
else
data_directory=['D:\Yawen\VBehavior\data\NBM_fMRI_Behavior_data\Staircase']; % set a new directory

end
% data_directory=['D:\Yawen\VBehavior\data\NBM_fMRI_Behavior'];
% data_directory=['D:\Yawen\VBehavior\data\NBM_fMRI_Behavior_data\Vonly']; % set a new directory
if ~exist(data_directory)
    mkdir(data_directory);
end

run_experiment=1;

KbName('UnifyKeyNames');
%% settings for the experiment
BlockCondition = {'Attend to Low Frequency',
    'Attend to High Frequency',
    'Attend to Left',
    'Attend to Right'};
cue = {'H','L','<+<','>+>'};
SDur  = 1.0; FDur  = 0.3; CDur = 1.0;respwindow = 3;
%% settings for visual stimulus
ngabors = 2;
% Phase of underlying sine grating in degrees:
phase = 0;
% Spatial constant of the exponential "hull"
sc = 10.0;
% Frequency of sine grating:
sf = .05;
% Contrast of grating:
contrast = 20;
% Aspect ratio width vs. height:
aspectratio = 1.0;


%% settings for the staircase
MaxTrials = 120;
MaxReverals = 14;
IgnoreReversals = 4;

% Define ranges of stimulus
InitialContra = 12;
MaxStimContr = 20;
ConstContr = 20;

% StepSize(1) = 5;
% for block = 2:15
%     StepSize(block) = StepSize(block-1)/1.2;
% end
% Define the type of staircase
UpNum = 1;    % number of incorrect answer to go one step up
DownNum = 4;  % number of correct answers to go one step down
% Ratio of Up and Down stepsize
ratio = 0.8415;
% load('1_test_20180613Vstaircase.mat');
ContrChg = [];TID=[];
%% create log file 
% colheaders = {'subId','subName','blockNum','blockorder','cue','trialNum','catchTrial','soundfile','ISI','VctPostion','pressedkey','rt','rttoexpstart'}

if run_experiment
    
    %%
    AssertOpenGL;
    
    %% Initiate psychtoolbox and basic setting up
    try
        HideCursor;
        Screen('Preference','SkipSyncTests',1);
        Priority(1);
        PsychImaging('PrepareConfiguration');
        screennumber = max(Screen('Screens'));
        PsychImaging('AddTask', 'General', 'FloatingPoint32BitIfPossible');
%           PsychImaging('AddTask','AllViews','RestrictProcessing',CenterRect([0 0 512 512],Screen('Rect',0)));
        [w,screenrect]=PsychImaging('OpenWindow',screennumber, 128,[]);
        Screen('BlendFunction',w,GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
        Screen('TextSize',w,20);
%         Screen('Preference', 'DefaultFontSize', 35);
%         Screen('Preference', 'DefaultFontStyle',1);
        Screen('Flip', w);
        ifi = Screen('GetFlipInterval',w);
        Screen('BlendFunction', w, GL_ONE, GL_ONE);
        % get the center position of the screen
        [cx,cy] = RectCenter(screenrect);
        
        eyetoscreen =62; %cm 7T nova coil 99
        ImageWidth = 40.7;  %cm 7T nova coil  30
        ImageHeight = 30; % cm 7T nova coil 18
        xdistance = 4; % deg
        pixperdeg= floor(.5*screenrect(4)/(atan2(0.5*ImageHeight,eyetoscreen)*(180/pi))); % size of the stimuli
        xdisPix = floor(xdistance*(.5*screenrect(3)/(atan2(0.5*ImageWidth,eyetoscreen)*(180/pi))));
        gaborw = 2*pixperdeg+1; gaborh = 2*pixperdeg + 1;
        textRect = [0; 0; 2*gaborw; 2*gaborh;];
        dstRects(:,1) = CenterRectOnPoint(textRect',cx-xdisPix,cy);% location of the stimuli
        dstRects(:,2) = CenterRectOnPoint(textRect',cx+xdisPix,cy);        % for play auditory stimulus
        %%      set up for eyelink
        if ET
            el=EyelinkInitDefaults(w);
            
            % Initialization of the connection with the Eyelink Gazetracker.
            % exit program if this fails.
            dummymode=0;
            if ~EyelinkInit(dummymode)
                fprintf('Eyelink Init aborted.\n');
                cleanup;  % cleanup function
                return;
            end
            % make sure that we get gaze data from the Eyelink
            Eyelink('command', 'link_sample_data = LEFT,RIGHT,GAZE,AREA');
            
            % open file to record data to
            
            res = Eyelink('openfile', edffile);%20180611
            if res~=0
                fprintf('Cannot create EDF file ''%s'' ', edffilename);
                cleanup;
                %         Eyelink( 'Shutdown');
                return;
            end
                        
            % make sure we're still connected.
            if Eyelink('IsConnected')~=1 && ~dummymode
                fprintf('Not connected. exiting');
                cleanup;
                return;
            end
            % STEP 4
            % Calibrate the eye tracker
            EyelinkDoTrackerSetup(el);
            
            % do a final check of calibration using driftcorrection
            EyelinkDoDriftCorrection(el);
            WaitSecs(0.1);
            Eyelink('StartRecording');
            eye_used = Eyelink('EyeAvailable'); % get eye that's tracked
            if eye_used == el.BINOCULAR; % if both eyes are tracked
                eye_used = el.LEFT_EYE; % use left eye
            end

        end
%% key of interest
        keysOfInterest = zeros(1,256);
        keysOfInterest(KbName({'5%','ESCAPE','space','q'})) = 1;
        KbQueueCreate(-1,keysOfInterest);
        %% fixation cross before start
        DrawFormattedText(w, 'Press Space Key To Begin', 'center', 'center', [255,255,255]);
        expstart=Screen('Flip', w);
        KbStrokeWait;    
        KbQueueFlush(-1);
        %% experimental loops
        SDurFrames = round(SDur/ifi);CDurFrames = round(CDur/ifi); FDurFrames = round(FDur/ifi);respFrames = round(2/ifi);
        gabortex = CreateProceduralGabor(w,gaborw,gaborh,1);
        trialcounter = 0;
        % frame to change the contrast
        loc = round(4400/44100/ifi); % time to change, correspond to the center of individual beep
        if ET
            Eyelink('Message','EVENTID %s', 'Begin');
            eye_data = [];
        end

        for block = 3:4
             if ET
                 Eyelink('Message','BLOCKID %d',block);
            end
            blocktext = cue{block};
            DrawFormattedText(w, blocktext, 'center', 'center', 255);
            blockstartvbl = Screen('Flip',w);
            WaitSecs(5.2);
            TrialIdx = ones(1,11);
            TrialIdx(1:11) = 1;
            revrs = []; chgcontrast = [];corrans = [];count=0;UpFlag=0; DownFlag=0;
            revrs_trial = [];
            trial=0;forder = zeros(1,11);
            forder(1:11) = [1 1 1 0 1 0 1 0 1 0 1];
            chgcontrast(1:11) = [1,10,5,20,18,20,13,20,7,20,InitialContra];

        while length(revrs) <= MaxReverals  && trial <= MaxTrials 
            x=[];    y=[]; t=[];t_ix=[];
            trialcounter = trialcounter + 1;
            trial = trial + 1;
            if ET
                Eyelink('Message', 'TRIALID %d', trial);
                % This supplies the title at the bottom of the eyetracker display
                Eyelink('command', 'record_status_message "TRIAL %d/%d"',trial, 120);
                
            end

            rotAngles = [rand * 135 rand*45];
            phase = randi(180);
            whichbeep = randi([2,4],1);
            whichframe = whichbeep*loc;
            
            tstarttime= GetSecs;
            KbQueueStart(-1);
            DrawFormattedText(w, '+', 'center','center',255);
            vbl = Screen('Flip', w,tstarttime +0.5*ifi);
            if ET
                Eyelink('Message', 'EVENTID %s', 'F');
            end
            if ET
                wait_fixation = 1;
                while wait_fixation
                    if Eyelink( 'NewFloatSampleAvailable') > 0
                        % get the sample in the form of an event structure
                        evt = Eyelink( 'NewestFloatSample');
                        if eye_used ~= -1 % do we know which eye to use yet?
                            % if we do, get current gaze position from sample
                            if (abs(evt.gx(eye_used+1)-cx)/pixperdeg)<5  && (abs(evt.gy(eye_used+1)-cy)/pixperdeg)<5
                                wait_fixation=0;
                            end
                        end
                    end
                    DrawFormattedText(w, '+', 'center','center',128);
                    vbl = Screen('Flip', w);
                end
                eyelink_times(1) =Eyelink('TrackerTime');
                break_fixation=0;tic
            end
            
            %                 break_fixation=0;%tic %% comment it
            
            for frame = 1:FDurFrames - 1
                % Draw the fixation point
                DrawFormattedText(w, '+', 'center','center',255);
                % Flip to the screen
                vbl = Screen('Flip', w, vbl + 0.5 * ifi);
                if ET
                    if Eyelink( 'NewFloatSampleAvailable') > 0
                        % get the sample in the form of an event structure
                        evt = Eyelink( 'NewestFloatSample');
                        if eye_used ~= -1 % do we know which eye to use yet?
                            % if we do, get current gaze position from sample
                            x(end+1) = evt.gx(eye_used+1); % +1 as we're accessing MATLAB array
                            y(end+1) = evt.gy(eye_used+1);
                            t(end+1)=toc;
                            t_ix(end+1)=1;
                            if (abs(evt.gx(eye_used+1)-cx)/pixperdeg)>3  && (abs(evt.gy(eye_used+1)-cy)/pixperdeg)>3
                                break_fixation=1;
                                beep
                                break
                            end
                            
                        end
                    end
                end
                
            end
            break_fixation = 0;
                
                if ~break_fixation
                    for frame = 1 : CDurFrames
                        if ET
                            if frame == 1
                                Eyelink('Message', 'EVENTID %s','C');
                            end
                        end
                        
                        %                 cueidx = randi([3 4],1);
                        %                         DrawFormattedText(w, cue{block}, 'center','center',255);
                        if ET
                            if Eyelink('NewFloatSampleAvailable') > 0
                                evt = Eyelink('NewestFloatSample');
                                if eye_used ~= -1 % do we know which eye to use yet?
                                    % if we do, get current gaze position from sample
                                    x(end+1) = evt.gx(eye_used+1); % +1 as we're accessing MATLAB array
                                    y(end+1) = evt.gy(eye_used+1);
                                    t(end+1)=toc;
                                    t_ix(end+1)=2;
                                    if (abs(evt.gx(eye_used+1)-cx)/pixperdeg)>3  && (abs(evt.gy(eye_used+1)-cy)/pixperdeg)>3
                                        break_fixation=1;
                                        beep
                                        break
                                    end
                                end
                            end
                        end
                        
                        
                        
                        contrast = 20;
                        mypars = repmat([phase+2*frame, sf, sc, contrast, aspectratio, 0, 0, 0]', 1, ngabors);
                        if block == 3 && forder(trial)>0 && frame>= whichframe && frame <  whichframe + 2
                            mypars = [[phase+2*frame, sf, sc, chgcontrast(trial), aspectratio, 0, 0, 0]',[phase+2*frame, sf, sc, 20, aspectratio, 0, 0, 0]'];
                            
                        end
                        
                        if block == 4 && forder(trial) > 0 && frame>= whichframe && frame <  whichframe + 2 
                            
                            mypars = [[phase+2*frame, sf, sc, 20, aspectratio, 0, 0, 0]',[phase+2*frame, sf, sc, chgcontrast(trial), aspectratio, 0, 0, 0]'];
                            
                        end
                        
                        DrawFormattedText(w,'+', 'center','center',255);
                        
                        Screen('DrawTextures', w, gabortex, [], dstRects, rotAngles, [], [], [], [], kPsychDontDoRotation, mypars);
                        
                        
                        vbl = Screen('Flip',w,vbl+0.5*ifi);
                    end
                end
            if ET
                eye_data{trialcounter}.x=x;    eye_data{trialcounter}.y=y; eye_data{trialcounter}.t=t;eye_data{trialcounter}.t_ix=t_ix;
            end
            
            % record response
            respStart = GetSecs;
            for frame  = 1 : respFrames
                if ET
                    if frame == 1
                        Eyelink('Message', 'EVENTID %s', 'R');
                    end
                end
                
                DrawFormattedText(w, '+', 'center','center',[0,255,0]);
                vbl = Screen('Flip',w,vbl+0.5*ifi);
            end
            [pressed, firstPress, firstRelease, lastPress, lastRelease] = KbQueueCheck(-1);
             if pressed
                       
                        if length(find(strcmp(KbName(find(lastPress)),'space')))>=1
                            exp_term=0;
                            if  forder(trial)>0
                                
                                 corrans(trial) = 1;
                            else
                                 corrans(trial) = 0;
                            end
                            pressedKey = 32;
                            presstime = lastPress(pressedKey) - respStart;
                            log(trialcounter).allpressedkeyname = KbName(find(lastPress));
                            log(trialcounter).allpressedkey = find(lastPress);
                            log(trialcounter).allpressedkeyrt = firstPress(pressedKey) - respStart;
                        end

                if strcmp(KbName(find(lastPress)),'q')
                            exp_term=1;
                            pressedKey = 27;
                            presstime = lastPress(pressedKey) - respStart;

                            if exp_term
                                %
                                log(trialcounter).subjnum = str2num(subjnum);
                                log(trialcounter).currblock = block;
                                log(trialcounter).currtrial = trial;
                                log(trialcounter).blockstartime = blockstartvbl-expstart;
                                log(trialcounter).trialstartime = tstarttime - blockstartvbl;
                                log(trialcounter).fixationstart = tstarttime - blockstartvbl + 0.5*ifi;
                                log(trialcounter).cuestart = tstarttime - blockstartvbl +(FDurFrames-0.5)*ifi;
                                log(trialcounter).stimulustart = tstarttime - blockstartvbl +(FDurFrames +19 -0.5)*ifi;
                                log(trialcounter).correctans = corrans(trial);
                                log(trialcounter).key = pressedKey;
                                log(trialcounter).rt = presstime;
                                log(trialcounter).dur = (GetSecs-tstarttime);
                                log(trialcounter).forder = forder(trial);
                                log(trialcounter).frames = whichframe +19;
                                log(trialcounter).dur = (GetSecs-tstarttime);
                                log(trialcounter).chgcontrast = chgcontrast(trial);
                                log(trialcounter).allpressedkeyname = KbName(find(lastPress));
                                log(trialcounter).allpressedkey = find(lastPress);
                                log(trialcounter).allpressedkeyrt = lastPress(pressedKey) - respStart;
                                
                                tmp = [str2num(subjnum),block,trial,forder(trial),chgcontrast(trial),blockstartvbl-expstart,tstarttime - blockstartvbl, 0.5*ifi, (FDurFrames-0.5)*ifi, (FDurFrames+19-0.5)*ifi,whichframe,corrans(trial),pressedKey,presstime,(GetSecs-tstarttime),rotAngles,phase];
                                %                                     tmp = [cnd(trial,1),block,cnd(trial,2),cnd(trial,3),cnd(trial,4),cnd(trial,5),cnd(trial,6),blockstartvbl-expstart,tstarttime, 0.5*ifi, (FDurFrames-0.5)*ifi, (FDurFrames+19-0.5)*ifi,NAndG{cnd(trial,4)}.bpos,vctpos,corrans, pressedKey,presstime,(GetSecs-tstarttime)];
                                results(trialcounter) = tmp;
                                cd(data_directory)
                                
                                if ET
                                    save(sprintf('%s%s%s',filename,'Vstaircase','.mat'),'results','eye_data','log','R','Threshold','C','ContrChg','TID');
                                    Eyelink('Message', 'EVENTID %s', 'Q');
                                    Eyelink('StopRecording');
                                    Eyelink('CloseFile');
                                    status = Eyelink('ReceiveFile');
                                    if status > 0
                                        fprintf('ReceiveFile status %d\n', status);
                                    end
                                    if 2==exist(edfFile, 'file')
                                        fprintf('Data file ''%s'' can be found in ''%s''\n', edfFile, pwd );
                                    end
                                    
                                else
                                    
                                    save(sprintf('%s%s%s',filename,'Vstaircase','.mat'),'results','log','R','Threshold','C','ContrChg','TID');
                                end
                                Priority(0);
                                ShowCursor;
                                Screen('CloseAll');
                                return;
                            end
                        end

             else
                 pressedKey = 0;
                 presstime = 0;
                 break_fixation = 0;
                 log(trialcounter).allpressedkeyname = '0';
                 log(trialcounter).allpressedkey = 0;
                 log(trialcounter).allpressedkeyrt = 0;
                 
                 
                 if forder(trial)<=0
                     corrans(trial) = 1;
                 else
                     corrans(trial) = 0;
                 end
                 
                 
             end
             

           if  trial >= 11
                        trialidx = rand;
                        TrialIdx(trial+1)= trialidx;
                        if  TrialIdx(trial+1) <= 0.2
                            %contrast don't change
                            %the trial that don't change contrast won't count in the staircase 
                            forder(trial+1) = 0;
                            chgcontrast(trial+1) = chgcontrast(trial);

                        else
                            % Contrast change
                            forder(trial+1) = 1;
                            
                            if corrans(trial)
                                count = count + 1; DownFlag = 1;
                                chgcontrast(trial+1) = chgcontrast(trial);
                               
                                if count >= DownNum
                                    count = 0;
                                    chgcontrast(trial+1) = chgcontrast(trial)*1.1;
                                    
                                    %                                     chgcontrast(trial+1) = chgcontrast(trial)*1.2*ratio;
                                    if chgcontrast(trial+1) >= MaxStimContr
                                        chgcontrast(trial+1) = MaxStimContr*0.95;
                                    end
                                    if UpFlag == 1
                                        revrs = [revrs chgcontrast(trial)];
                                        revrs_trial = [revrs_trial trial];
                                        
                                    end
                                    UpFlag = 0;
                                end
                                
                            else
                                count = 0;
                                UpFlag = 1;
                                chgcontrast(trial+1) = chgcontrast(trial)/1.1;
                                if chgcontrast(trial+1) <= InitialContra
                                    chgcontrast(trial+1) = InitialContra;
                                end
                                if DownFlag == 1
                                    revrs = [revrs chgcontrast(trial)];
                                    revrs_trial = [revrs_trial trial];
                                end
                                DownFlag = 0;
                            end
                        end
           else
               chgcontrast(trial+1) = chgcontrast(trial+1);
           end
           
           
           
           log(trialcounter).subjnum = str2num(subjnum);
           log(trialcounter).currblock = block;
           log(trialcounter).currtrial = trial;
           log(trialcounter).blockstartime = blockstartvbl-expstart;
           log(trialcounter).trialstartime = tstarttime - blockstartvbl;
           log(trialcounter).fixationstart = tstarttime - blockstartvbl + 0.5*ifi;
           log(trialcounter).cuestart = tstarttime - blockstartvbl +(FDurFrames-0.5)*ifi;
           log(trialcounter).stimulustart = tstarttime - blockstartvbl +(FDurFrames +19 -0.5)*ifi;
           log(trialcounter).correctans = corrans(trial);
           log(trialcounter).key = pressedKey;
           log(trialcounter).rt = presstime;
           log(trialcounter).dur = (GetSecs-tstarttime);
           log(trialcounter).forder = forder(trial);
           log(trialcounter).frames = whichframe+19;
           log(trialcounter).dur = (GetSecs-tstarttime);
           log(trialcounter).chgcontrast = chgcontrast(trial);

           
           tmp = [str2num(subjnum),block,trial,forder(trial),chgcontrast(trial),blockstartvbl-expstart,tstarttime - blockstartvbl, 0.5*ifi, (FDurFrames-0.5)*ifi, (FDurFrames+19-0.5)*ifi,whichframe+19,corrans(trial),pressedKey,presstime,(GetSecs-tstarttime),rotAngles,phase];
           results(trialcounter,1:length(tmp)) = tmp;
           KbQueueFlush(-1);
        end
        if ET
            Eyelink('StopRecording');
        end
        
        TID(block-2).trialid = TrialIdx;
        ContrChg(block-2).conchg = chgcontrast;
        R(block-2).revrs = revrs; R(block-2).revrstrial = revrs_trial;
        Threshold(block-2).threshold = mean(revrs(end-10:end));
        C(block-2).correctans  = corrans;
        DrawFormattedText(w, 'Have a Rest.Press any key to start', 'center','center',255);
        blockendtime = Screen('Flip', w); blockdur = blockendtime - blockstartvbl;
        bend(block) = blockendtime; bdur(block) = blockdur;
        KbStrokeWait;
        if ET
            EyelinkDoTrackerSetup(el);
            KbStrokeWait;
            Eyelink('StartRecording');
        else
            KbStrokeWait;
        end
        WaitSecs(0.1);
        cd(data_directory)
        
        
        if ET
            save(sprintf('%s%s%s',filename,'Vstaircase','.mat'),'results','log','R','eye_data','Threshold','C','ContrChg','TID');
        else
            save(sprintf('%s%s%s',filename,'Vstaircase','.mat'),'results','log','R','Threshold','C','ContrChg','TID');
            
        end

        end
        
        DrawFormattedText(w, 'End of the experiments', 'center','center',255);
        endexp = Screen('Flip', w);
        KbStrokeWait;
        
        expdur = endexp - expstart;
        cd(data_directory)
                KbQueueRelease(-1);
        if ET
            Eyelink('CloseFile');
            status = Eyelink('ReceiveFile');
            if status > 0
                fprintf('ReceiveFile status %d\n', status);
            end
            if 2==exist(edfFile, 'file')
                fprintf('Data file ''%s'' can be found in ''%s''\n', edfFile, pwd );
            end
        end

        save(sprintf('%s%s%s',filename,'Vstaircase','.mat'),'results','log','R','Threshold','C','expdur','ContrChg','TID');
        Screen('CloseAll');

            
    catch
        Screen('CloseAll');
end
end


