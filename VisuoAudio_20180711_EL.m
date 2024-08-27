%% Beginning
close all; clear all;
subjnum = input('Please input subject number:', 's');
subjno=input('Please input date:', 's');
subjname=input('Please input subject name:', 's');
SettingsForVisuoAudio20180711
% load('minyecondition.mat')
SE_Frequency
KbName('UnifyKeyNames');
ET = 1;
filename = sprintf('%s_%s_%s',subjnum,subjname,subjno);
if ET    
edffile = [upper(subjname) '.edf'];
% eyedata_directory=['D:\Yawen\VBehavior\data\NBM_fMRI_Behavior_data\Eyedata'];
end
data_directory='D:\Yawen\VBehavior\data\NBM_fMRI_Behavior_data\AudioVisual'; % set a new directory
% load('DavidConditons.mat') 201080615 (start from block =5 ,from second time)
% data_directory=['D:\Yawen\VBehavior\data\NBM_fMRI_Behavior']; % set a new directory
if ~exist(data_directory)
    mkdir(data_directory);
end

%% create log file
% colheaders = {'subId','subName','blockNum','blockorder','cue','trialNum','catchTrial','soundfile','ISI','VctPostion','pressedkey','rt','rttoexpstart'}



%%
AssertOpenGL;

%% Initiate psychtoolbox and basic setting up
try
    HideCursor;
    Priority(1);
    PsychImaging('PrepareConfiguration');
    screennumber = max(Screen('Screens'));
    PsychImaging('AddTask', 'General', 'FloatingPoint32BitIfPossible');
    %   PsychImaging('AddTask','AllViews','RestrictProcessing',CenterRect([0 0 512 512],Screen('Rect',0)));
    [w,screenrect]=PsychImaging('OpenWindow',screennumber, 128,[]);
    Screen('BlendFunction',w,GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
    Screen('TextSize',w,32);
    %         Screen('Preference', 'DefaultFontSize', 35);
    %         Screen('Preference', 'DefaultFontStyle',1);
    Screen('Flip', w);
    ifi = Screen('GetFlipInterval',w);
    loc = round(4400/44100/ifi);
    % Enable alpha-blending, set it to a blend equation useable for linear
    % superposition with alpha-weighted source. This allows to linearly
    % superimpose gabor patches in the mathematically correct manner, should
    % they overlap. Alpha-weighted source means: The 'globalAlpha' parameter in
    % the 'DrawTextures' can be used to modulate the intensity of each pixel of
    % the drawn patch before it is superimposed to the framebuffer image, ie.,
    % it allows to specify a global per-patch contrast value:
    Screen('BlendFunction', w, GL_ONE, GL_ONE);
    % get the center position of the screen
    [cx,cy] = RectCenter(screenrect);
    eyetoscreen =62; %cm 7T nova coil 99
    ImageWidth = 40.7;  %cm 7T nova coil  30
    ImageHeight = 30; % cm 7T nova coil 18
    xdistance = 4; % deg
    pixperdeg= floor(.5*screenrect(4)/(atan2(0.5*ImageHeight,eyetoscreen)*(180/pi))); % size of the stimuli
    xdisPix = floor(xdistance*(.5*screenrect(3)/(atan2(0.5*ImageWidth,eyetoscreen)*(180/pi))));
    Screen('TextSize',w,20);
    gaborw = 2*pixperdeg+1; gaborh = 2*pixperdeg + 1;
    textRect = [0; 0; 2*gaborw; 2*gaborh;];
    dstRects(:,1) = CenterRectOnPoint(textRect',cx-xdisPix,cy);% location of the stimuli
    dstRects(:,2) = CenterRectOnPoint(textRect',cx+xdisPix,cy);        % for play auditory stimulus
    gabortex = CreateProceduralGabor(w, gaborw,gaborh);
    InitializePsychSound;
    pahandle = PsychPortAudio('Open', [], [], 0, Fs, 1);
    
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
    
    %%
    
    keysOfInterest = zeros(1,256);
    keysOfInterest(KbName({'5%',responsekeys})) = 1; %20180411 change
    KbQueueCreate(-1,keysOfInterest);
    KbQueueStart(-1);
    
    %         PsychHID('KbQueueCreate',-1,keysOfInterest);
    
    %% fixation cross before start
    DrawFormattedText(w, 'Press Any Key To Begin', 'center', 'center', [255,255,255]);
    expstart=Screen('Flip', w);
    if ET
        Eyelink('Message','EVENTID %s', 'Begin')
    end
    
    KbStrokeWait;
    KbQueueFlush(-1);
    
    %         PsychHID('KbQueueFlush',-1);
    
    % text = 'Please keep your head still and focus on the center cross and cue.\n\nPress "up" when you detect an increase of illuminance or gap in high frequency.\n\nPress "down" when you detect an decrease of illuminance of the grating or gap in low frequency. \n\n';
    % wt=RectWidth(Screen('TextBounds',w,text));
    % ht=RectHeight(Screen('TextBounds',w,text));
    % Screen('DrawText',w,text,cx,cy,[255 255 255]);
    % Screen('Flip',w);
    
    %% experimental loops
    FDur = 0.3;
    SDurFrames = round(SDur/ifi);CDurFrames = round(CDur/ifi); respFrames = round(1.5/ifi); % same as the scanner
    FDurFrames = round(FDur/ifi);
    
    for  block  = 1 : size(Conditions,3)
        if ET
            Eyelink('Message','BLOCKID %d',block);
        end
        cnd = Conditions(:,:,block);
        blocktext = cue{cnd(1,2)};
        DrawFormattedText(w, blocktext, 'center', 'center', 255);
        blockstartvbl = Screen('Flip',w);
        bstart(block) = blockstartvbl -expstart;
        gabortex = CreateProceduralGabor(w,gaborw,gaborh,1);
        WaitSecs(5.2);
        %       blobtex = CreateProceduralGaussBlob(w,gaborw*2,gaborh*2);
        for  trial =  1: size(Conditions,1)
            trialcounter = nTrials*(block-1)+trial;
            if ET
                Eyelink('Message', 'TRIALID %d', trial);
                % This supplies the title at the bottom of the eyetracker display
                Eyelink('command', 'record_status_message "TRIAL %d/%d/%d"',trial, size(Conditions,1),block);
                
            end
            x=[];    y=[]; t=[];t_ix=[];actframe = 0;corrans=0;
            rotAngles = cnd(trial,9:10);
            %                 rotAngles = [rand * 135 rand*45];
            
            whichframe = 0;
            if cnd(trial,5)
                tsound = NewSound(cnd(trial,5)).sound; % low change
            end
            if cnd(trial,4)
                tsound = NewSound(cnd(trial,4)).sound; % high change
            end
            if ~cnd(trial,5) && ~cnd(trial,4)
                tsound = EqSound(1).sound;
            end
            phase = cnd(trial,8);
            ps = mean(tsound);% change the sound here
            KbQueueStart(-1);
            %                 PsychHID('KbQueueStart',-1);
            tstarttime= GetSecs;
            DrawFormattedText(w, '+', 'center','center',255);
            vbl = Screen('Flip', w,tstarttime +0.5*ifi);
            fixstarttime = vbl;
            %                 VPX_SendCommand('dataFile_InsertMarker F');
            
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
                            if (abs(evt.gx(eye_used+1)-cx)/pixperdeg)<3  && (abs(evt.gy(eye_used+1)-cy)/pixperdeg)<3
                                wait_fixation=0;
                            end
                        end
                    end
                    DrawFormattedText(w, '+', 'center','center',255);
                    vbl = Screen('Flip', w);
                end
                eyelink_times(1) =Eyelink('TrackerTime');
                break_fixation=0;tic
            end
            %                 break_fixation=0;
            %                 if ET
            %                     tic
            %                 end %% comment it
            % Now we present the isi interval with fixation point minus one frame because we presented the fixation point once already when getting a time stamp
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
                            if (abs(evt.gx(eye_used+1)-cx)/pixperdeg)>10  && (abs(evt.gy(eye_used+1)-cy)/pixperdeg)>10
                                break_fixation=1;
                                beep
                                break
                            end
                            
                        end
                    end
                end
            end
            
            
            if ~break_fixation
                PsychPortAudio('FillBuffer', pahandle, ps); % loads data into buffer
                soundstarttime = PsychPortAudio('Start', pahandle,1,vbl+0.5*ifi,1);
                %                     PsychPortAudio('Start', pahandle,1,inf); %starts sound at infinity
                %                     PsychPortAudio('RescheduleStart', pahandle, vbl+0.5*ifi, 0) %reschedules startime to
                %
                for frame = 1 : CDurFrames
                    if ET
                        if frame == 1
                            Eyelink('Message', 'EVENTID %s', 'C');
                        end
                        if Eyelink('NewFloatSampleAvailable') > 0
                            evt = Eyelink('NewestFloatSample');
                            if eye_used ~= -1 % do we know which eye to use yet?
                                % if we do, get current gaze position from sample
                                x(end+1) = evt.gx(eye_used+1); % +1 as we're accessing MATLAB array
                                y(end+1) = evt.gy(eye_used+1);
                                t(end+1)=toc;
                                t_ix(end+1)=2;
                                if (abs(evt.gx(eye_used+1)-cx)/pixperdeg)>10  && (abs(evt.gy(eye_used+1)-cy)/pixperdeg)>10
                                    break_fixation=1;
                                    beep
                                    break
                                end
                            end
                        end
                    end
                    
                    
                    contrast = 20;
                    mypars = repmat([phase+2*frame, sf, sc, contrast, aspectratio, 0, 0, 0]', 1, ngabors);
                    if cnd(trial,6)
                        whichframe = cnd(trial,6)*loc;
                        if frame>= whichframe && frame <  whichframe + 2
                            contrast = 15.74;%13.29;
                            actframe = frame;
                            mypars = [[phase+2*frame, sf, sc, contrast, aspectratio, 0, 0, 0]',[phase+2*frame, sf, sc, 20, aspectratio, 0, 0, 0]'];
                        end
                    end
                    
                    if cnd(trial,7)
                        whichframe = cnd(trial,7)*loc;
                        if frame>= whichframe && frame <  whichframe + 2
                            contrast = 17.00;%15.86;
                            actframe = frame;
                            mypars = [[phase+2*frame, sf, sc, 20, aspectratio, 0, 0, 0]',[phase+2*frame, sf, sc, contrast, aspectratio, 0, 0, 0]'];
                        end
                    end
                    DrawFormattedText(w, '+', 'center','center',255);
                    Screen('DrawTextures', w, gabortex, [], dstRects, rotAngles, [], [], [], [], kPsychDontDoRotation, mypars);
                    vbl = Screen('Flip',w,vbl+0.5*ifi);
                end
            end
            if ET
                eye_data{trialcounter}.x=x;    eye_data{trialcounter}.y=y; eye_data{trialcounter}.t=t;eye_data{trialcounter}.t_ix=t_ix;
            end
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
            
            %                         [secs, keyCode, deltaSecs] = KbStrokeWait(-1);
            [pressed, firstPress, firstRelease, lastPress, lastRelease] = KbQueueCheck(-1);
            if pressed
                if  length(find(strcmp(KbName(find(lastPress)),'space')))>=1
                    %                             if ET
                    %                                 vpx_SendCommandString('dataFile_InsertMarker K');
                    %                             end
                    exp_term=0;
                    if cnd(trial,2) == 1 && cnd(trial,5)%low
                        corrans = 1;
                    elseif cnd(trial,2) == 2 && cnd(trial,4)% high
                        corrans = 1;
                    elseif cnd(trial,2) == 3 && cnd(trial,6)%  left
                        corrans = 1;
                    elseif cnd(trial,2) == 4 && cnd(trial,7)%  right
                        corrans = 1;
                    end
                    if ~exp_term
                        % 20180411
                        pressedKey = 32;
                        keyname = KbName(pressedKey);
                        presstime = lastPress(32) - respStart;
                        log{trialcounter}.allpressedkeyname = KbName(find(lastPress));
                        log{trialcounter}.allpressedkey = find(lastPress);
                        log{trialcounter}.allpressedkeyrt = lastPress(pressedKey) - respStart;
                        
                    end
                end
                if strcmp(KbName(find(firstPress)),'q')
                    exp_term=1;
                    
                    
                    trialendtime = GetSecs;
                    
                    log{trialcounter}.subjnum = cnd(trial,1);
                    log{trialcounter}.currblock = block;
                    log{trialcounter}.blockcond = cnd(trial,2);
                    log{trialcounter}.currtrial = cnd(trial,3);
                    log{trialcounter}.lowcond = cnd(trial,4);
                    log{trialcounter}.highcond = cnd(trial,5);
                    log{trialcounter}.leftcond = cnd(trial,6);
                    log{trialcounter}.rightcond = cnd(trial,7);
                    log{trialcounter}.leftori = cnd(trial,8);
                    log{trialcounter}.rightori = cnd(trial,9);
                    log{trialcounter}.phase     = cnd(trial,10);
                    log{trialcounter}.blockstartime = blockstartvbl-expstart;
                    log{trialcounter}.trialstartime = tstarttime - expstart;
                    log{trialcounter}.fixstartime = fixstarttime -expstart;
                    log{trialcounter}.stimulustart = soundstarttime -expstart;
                    log{trialcounter}.respstart = respStart-expstart;
                    log{trialcounter}.trialendtime = trialendtime - expstart;
                    log{trialcounter}.frames = whichframe;
                    log{trialcounter}.actframe = actframe;
                    log{trialcounter}.correctans = corrans;
                    log{trialcounter}.key = pressedKey;
                    log{trialcounter}.rt = presstime;
                    log{trialcounter}.dur = (GetSecs-tstarttime);
                    log{trialcounter}.allpressedkeyname = KbName(find(lastPress));
                    log{trialcounter}.allpressedkey = find(lastPress);
                    log{trialcounter}.allpressedkeyrt = lastPress(pressedKey) - respStart;
                    
                    tmp = [cnd(trial,1),block,cnd(trial,2:10),blockstartvbl-expstart,tstarttime-expstart,fixstarttime -expstart,soundstarttime -expstart, respStart-expstart, trialendtime - expstart,whichframe,actframe,corrans, pressedKey,presstime,iti(trial)];
                    results(trialcounter,1:length(tmp)) = tmp;
                    
                    cd(data_directory);
                    if ET
                        save([filename,'.mat'],'results','log','eye_data','NewSound','AllConditions');
                        save(sprintf('%s%s%s',filename,'blockinfo','.mat'),'bend','bdur','expdur','bstart');
                        Eyelink('Message', 'EVENTID %s', 'Q');
                        Eyelink('CloseFile');
                        status = Eyelink('ReceiveFile');
                        if status > 0
                            fprintf('ReceiveFile status %d\n', status);
                        end
                        if 2==exist(edfFile, 'file')
                            fprintf('Data file ''%s'' can be found in ''%s''\n', edfFile, pwd );
                        end
                        
                    else
                        save(sprintf('%s%s%s',filename,'blockinfo','.mat'),'bend','bdur','expdur','bstart');
                        save([filename,'.mat'],'results','log','NewSound','AllConditions');
                    end
                    Priority(0);
                    ShowCursor;
                    Screen('CloseAll');
                    return;
                end
                
            else
                pressedKey = 0;
                presstime = 0;
                if cnd(trial,2) == 1 && ~cnd(trial,5)% low
                    corrans = 1;
                elseif cnd(trial,2) == 2 && ~cnd(trial,4)% high
                    corrans = 1;
                elseif cnd(trial,2) == 3 && ~cnd(trial,6)%  left
                    corrans = 1;
                elseif cnd(trial,2) == 4 && ~cnd(trial,7)%  right
                    corrans = 1;
                end
                
                if length(find(cnd(trial,4:7)==0)) ==4
                    corrans = 1;
                end
                
                %break_fixation = 1;
                log{trialcounter}.allpressedkeyname = '0';
                log{trialcounter}.allpressedkey = 0;
                log{trialcounter}.allpressedkeyrt = 0;
                
                %                     WaitSecs(ISI(trial)); % gaya did
            end
            WaitSecs(iti(trial));
            trialendtime = GetSecs;
            
            log{trialcounter}.subjnum = cnd(trial,1);
            log{trialcounter}.currblock = block;
            log{trialcounter}.blockcond = cnd(trial,2);
            log{trialcounter}.currtrial = cnd(trial,3);
            log{trialcounter}.lowcond = cnd(trial,4);
            log{trialcounter}.highcond = cnd(trial,5);
            log{trialcounter}.leftcond = cnd(trial,6);
            log{trialcounter}.rightcond = cnd(trial,7);
            log{trialcounter}.leftori = cnd(trial,8);
            log{trialcounter}.rightori = cnd(trial,9);
            log{trialcounter}.phase     = cnd(trial,10);
            log{trialcounter}.blockstartime = blockstartvbl-expstart;
            log{trialcounter}.trialstartime = tstarttime - expstart;
            log{trialcounter}.fixstartime = fixstarttime -expstart;
            log{trialcounter}.stimulustart = soundstarttime -expstart;
            log{trialcounter}.respstart = respStart-expstart;
            log{trialcounter}.trialendtime = trialendtime - expstart;
            log{trialcounter}.frames = whichframe;
            log{trialcounter}.actframe = actframe;
            log{trialcounter}.correctans = corrans;
            log{trialcounter}.key = pressedKey;
            log{trialcounter}.rt = presstime;
            log{trialcounter}.dur = (GetSecs-tstarttime);
        
            tmp = [cnd(trial,1),block,cnd(trial,2:10),blockstartvbl-expstart,tstarttime-expstart,fixstarttime-expstart,soundstarttime-expstart,respStart-expstart,trialendtime-expstart,whichframe,actframe,corrans, pressedKey,presstime,iti(trial)];
            results(trialcounter,1:length(tmp)) = tmp;
            
            cd(data_directory)
            if ET
                save([filename,'.mat'],'results','log','eye_data','NewSound','AllConditions');
            else
                save([filename,'.mat'],'results','log','NewSound','AllConditions');
                
            end
            KbQueueFlush(-1);
            
            
        end
        if ET
            Eyelink('StopRecording');
        end
        DrawFormattedText(w, 'Have a Rest.Press any key to start', 'center','center',255);
        vbl = Screen('Flip', w);blockendtime = GetSecs;%blockendtime = KbStrokeWait;
        blockdur = blockendtime - blockstartvbl;
        bend(block) = blockendtime -expstart;; bdur(block) = blockdur;
        KbStrokeWait;
        if ET
            EyelinkDoTrackerSetup(el);
            KbStrokeWait;
            Eyelink('StartRecording');
            KbStrokeWait;
        else
            KbStrokeWait;
        end
        WaitSecs(0.1);
        %             KbStrokeWait;
        
    end
    DrawFormattedText(w, 'End of the experiments,Press any key to exit', 'center','center',255);
    endexp = Screen('Flip', w);
    KbStrokeWait;
    expdur = endexp - expstart;
    save(sprintf('%s%s%s',filename,'blockinfo','.mat'),'bend','bdur','expdur','bstart');
    KbQueueRelease(-1);
    
    %         PsychHID('KbQueueRelease',-1);
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
    Screen('CloseAll');
catch
    Screen('CloseAll');
end


% block = 8;
if block == 8
    i = 1;
    rtemp = results;
    highCondition    = rtemp(find(rtemp(:,3)==2),:);
    NumOfChgInHigh   = highCondition(find(highCondition(:,5)),:); % all the trials that subj should press the key
    NumOfNChgInHigh = highCondition(find(highCondition(:,5)==0),:);
    
    NumOfHitInHigh   = highCondition(find(highCondition(:,21)),:); % total hit numbers including the false hit
    NumOfCorrRespForHigh = highCondition(find(highCondition(:,20)),:); % all the correct response in high condition
    NumOfCorrHitForHigh = NumOfChgInHigh(find(NumOfChgInHigh(:,21)),:);
    NumOfCorrRejForHigh = NumOfNChgInHigh(find(NumOfNChgInHigh(:,21)==0),:);
    NumOfMissHitForHigh = NumOfChgInHigh(find(NumOfChgInHigh(:,21)==0),:);% True negetive
    
    NumOfHitLowInHigh = NumOfHitInHigh(find(NumOfHitInHigh(:,6)),:);
    NumOfHitLeftInHigh = NumOfHitInHigh(find(NumOfHitInHigh(:,7)),:);
    NumOfHitRightInHigh = NumOfHitInHigh(find(NumOfHitInHigh(:,8)),:);
    
    Hit(i).highCondition = highCondition;
    Hit(i).hithigh = NumOfHitInHigh;
    Hit(i).chginhigh = NumOfChgInHigh;
    Hit(i).corrhithigh = NumOfCorrHitForHigh;
    Hit(i).misshigh = NumOfMissHitForHigh;
    Hit(i).corresphigh = NumOfCorrRespForHigh;
    
    Hit(i).corrhitratehigh = size(NumOfCorrHitForHigh,1)/size(NumOfHitInHigh,1); % correct hit rate in the high
    Hit(i).corratehigh = size(NumOfCorrRespForHigh,1)/size(highCondition,1); % correct response in the high condition (correct hit and rejection)
    Hit(i).falsehithigh = 1 - size(NumOfCorrHitForHigh,1)/size(NumOfHitInHigh,1);
    Hit(i).missratehigh = size(NumOfMissHitForHigh,1)/size(NumOfChgInHigh,1);
    Hit(i).trueposhigh  = size(NumOfCorrHitForHigh,1)/size(NumOfChgInHigh,1);
    Hit(i).truenegtivehigh = size(NumOfCorrRejForHigh,1)/size(NumOfNChgInHigh,1);
    
    Hit(i).lowInhigh= size(NumOfHitLowInHigh,1)/size(NumOfHitInHigh,1);
    Hit(i).leftInhigh= size(NumOfHitLeftInHigh,1)/size(NumOfHitInHigh,1);
    Hit(i).rightInhigh= size(NumOfHitRightInHigh,1)/size(NumOfHitInHigh,1);
    
    % for condition 'attend to low frequency' and 'correct hit'
    lowCondition       = rtemp(find(rtemp(:,3)==1) ,:);
    % lowCondition       = rtemp(find(rtemp(:,2)==8),:); % only for david
    
    NumOfChgInLow      = lowCondition(find(lowCondition(:,6)),:);
    NumOfNChgInLow = lowCondition(find(lowCondition(:,6)==0),:);
    
    NumOfHitInLow     =  lowCondition(find(lowCondition(:,21)),:);
    NumOfCorrRespForLow = lowCondition(find(lowCondition(:,20)),:); % all the correct response in low condition
    NumOfCorrHitForLow = NumOfChgInLow(find(NumOfChgInLow(:,21)),:);% correct hit in all the hit trials
    NumOfCorrRejForLow = NumOfNChgInLow(find(NumOfNChgInLow(:,21)==0),:);
    
    NumOfMissHitForLow = NumOfChgInLow(find(NumOfChgInLow(:,21)==0),:);% miss rate
    
    NumOfHitHighInLow = NumOfHitInLow(find(NumOfHitInLow(:,5)),:);
    NumOfHitLeftInLow = NumOfHitInLow(find(NumOfHitInLow(:,7)),:);
    NumOfHitRightInLow = NumOfHitInLow(find(NumOfHitInLow(:,8)),:);
    
    Hit(i).lowCondition = lowCondition;
    Hit(i).hitlow = NumOfHitInLow;
    Hit(i).chginlow = NumOfChgInLow;
    Hit(i).corrhitlow = NumOfCorrHitForLow;
    Hit(i).misslow = NumOfMissHitForLow;
    Hit(i).corresplow =  NumOfCorrRespForLow;
    
    Hit(i).corrhitratelow  = size(NumOfCorrHitForLow,1)/size(NumOfHitInLow,1);%correct hit ratio in all the hit trials
    Hit(i).corratelow = size(NumOfCorrRespForLow,1)/size(lowCondition,1);% correct response rate in the low condition (correct hit and rejection)
    Hit(i).falsehitlow = 1 - (size(NumOfCorrHitForLow,1)/size(NumOfHitInLow,1));
    Hit(i).missratelow = size(NumOfMissHitForLow,1)/size(NumOfChgInLow,1);
    Hit(i).trueposlow  = size(NumOfCorrHitForLow,1)/size(NumOfChgInLow,1);
    Hit(i).truenegtivelow = size(NumOfCorrRejForLow,1)/size(NumOfNChgInLow,1);
    
    Hit(i).highInlow= size(NumOfHitHighInLow,1)/size(NumOfHitInLow,1);
    Hit(i).leftInlow= size(NumOfHitLeftInLow,1)/size(NumOfHitInLow,1);
    Hit(i).rightInlow= size(NumOfHitRightInLow,1)/size(NumOfHitInLow,1);
    
    
    % for condition 'attend to left' and 'correct hit'
    leftCondition = rtemp(find(rtemp(:,3)==3),:);
    NumOfCchgInLeft  = leftCondition(find(leftCondition(:,7)),:);
    NumOfNChgInLeft = leftCondition(find(leftCondition(:,7)==0),:);
    
    NumOfHitInLeft = leftCondition(find(leftCondition(:,21)),:);
    NumOfCorrRespForLeft = leftCondition(find(leftCondition(:,20)),:); % all the correct response in low condition
    NumOfCorrHitForLeft =  NumOfCchgInLeft(find(NumOfCchgInLeft(:,21)),:);% correct hit in all the hit trials
    NumOfCorrRejForLeft = NumOfNChgInLeft(find(NumOfNChgInLeft(:,21)==0),:);
    
    NumOfMissHitForLeft = NumOfCchgInLeft(find(NumOfCchgInLeft(:,21)==0),:);% Miss rate
    
    NumOfHitHighInLeft = NumOfHitInLeft(find(NumOfHitInLeft(:,5)),:);
    NumOfHitLowInLeft = NumOfHitInLeft(find(NumOfHitInLeft(:,6)),:);
    NumOfHitRightInLeft = NumOfHitInLeft(find(NumOfHitInLeft(:,8)),:);
    
    
    Hit(i).leftcondition = leftCondition;
    Hit(i).hitleft = NumOfHitInLeft;
    Hit(i).cchgleft = NumOfCchgInLeft;
    Hit(i).corrhitleft = NumOfCorrHitForLeft;
    Hit(i).missleft = NumOfMissHitForLeft;
    Hit(i).correspleft = NumOfCorrRespForLeft;
    
    Hit(i).corrhitrateleft  = size(NumOfCorrHitForLeft,1)/size(NumOfHitInLeft,1);
    Hit(i).corrateleft = size(NumOfCorrRespForLeft,1)/size(leftCondition,1);
    Hit(i).falsehitleft = 1 -(size(NumOfCorrHitForLeft,1)/size(NumOfHitInLeft,1));
    Hit(i).missrateleft = size(NumOfMissHitForLeft,1)/size(NumOfCchgInLeft,1);
    Hit(i).trueposleft = size(NumOfCorrHitForLeft,1)/size(NumOfCchgInLeft,1);
    Hit(i).truenegtiveleft = size(NumOfCorrRejForLeft,1)/size(NumOfNChgInLeft,1);
    
    Hit(i).highInleft= size(NumOfHitHighInLeft,1)/size(NumOfHitInLeft,1);
    Hit(i).lowInleft= size(NumOfHitLowInLeft,1)/size(NumOfHitInLeft,1);
    Hit(i).rightInleft= size(NumOfHitRightInLeft,1)/size(NumOfHitInLeft,1);
    
    
    
    % for condition 'attend to right' and 'correct hit'
    rightCondition   = rtemp(find(rtemp(:,3)==4),:);
    % rightCondition   = rtemp(find(rtemp(:,2)==5),:);%only for david
    
    NumOfCchgInRight = rightCondition(find(rightCondition(:,8)),:);
    NumOfNChgInRight = rightCondition(find(rightCondition(:,8)==0),:);
    
    NumOfHitInRight = rightCondition(find(rightCondition(:,21)),:);
    NumOfCorrRespForRight = rightCondition(find(rightCondition(:,20)),:); % all the correct response in low condition
    NumOfCorrHitForRight = NumOfCchgInRight(find(NumOfCchgInRight(:,21)),:);% correct hit in all the hit trials
    NumOfCorrRejForRight = NumOfNChgInRight(find(NumOfNChgInRight(:,21)==0),:);
    
    NumOfMissHitForRight = NumOfCchgInRight(find(NumOfCchgInRight(:,21)==0),:);% True negetive
    
    NumOfHitHighInRight = NumOfHitInRight(find(NumOfHitInRight(:,5)),:);
    NumOfHitLowInRight = NumOfHitInRight(find(NumOfHitInRight(:,6)),:);
    NumOfHitLeftInRight = NumOfHitInRight(find(NumOfHitInRight(:,7)),:);
    
    
    Hit(i).rightCondition = rightCondition;
    Hit(i).hitright = NumOfHitInRight;
    Hit(i).corrhitright = NumOfCorrHitForRight;
    Hit(i).missright = NumOfMissHitForRight;
    Hit(i).cchgright = NumOfCchgInRight;
    
    Hit(i).corrhitrateright  = size(NumOfCorrHitForRight,1)/size(NumOfHitInRight,1);
    Hit(i).corrateright = size(NumOfCorrRespForRight,1)/size(rightCondition,1);
    Hit(i).falsehitright= 1 -(size(NumOfCorrHitForRight,1)/size(NumOfHitInRight,1));
    Hit(i).missrateright = size(NumOfMissHitForRight,1)/size(NumOfCchgInRight,1);
    Hit(i).trueposright = size(NumOfCorrHitForRight,1)/size(NumOfCchgInRight,1);
    Hit(i).truenegtiveright = size(NumOfCorrRejForRight,1)/size(NumOfNChgInRight,1);
    
    Hit(i).highInright= size(NumOfHitHighInRight,1)/size(NumOfHitInRight,1);
    Hit(i).lowInright= size(NumOfHitLowInRight,1)/size(NumOfHitInRight,1);
    Hit(i).leftInright= size(NumOfHitLeftInRight,1)/size(NumOfHitInRight,1);
    
    
    save(sprintf('%s%s%s',filename,'HitInfo','.mat'),'Hit');
    
    % %% plot figures
    figure(1)
    data_hl= [Hit(1).corratehigh,Hit(1).corratelow;
        Hit(1).corrhitratehigh,Hit(1).corrhitratelow;
        Hit(1).trueposhigh,Hit(1).trueposlow;
        Hit(1).truenegtivehigh,Hit(1).truenegtivelow;
        Hit(1).missratehigh,Hit(1).missratelow;
        Hit(1).falsehithigh,Hit(1).falsehitlow];
    bar(data_hl)
    set(gca,'XTickLabel',{'CorrectRate','CorrectHitRate','TruePositive','TrueNegetive','MissRate','FalsePositive'})
    legend('High','Low')
    figure(2)
    data_lr= [Hit(1).corrateleft,Hit(1).corrateright;
        Hit(1).corrhitrateleft,Hit(1).corrhitrateright;
        Hit(1).trueposleft,Hit(1).trueposright;
        Hit(1).truenegtiveleft,Hit(1).truenegtiveright;
        Hit(1).missrateleft,Hit(1).missrateright;
        Hit(1).falsehitleft,Hit(1).falsehitright];
    bar(data_lr)
    set(gca,'XTickLabel',{'CorrectRate','CorrectHitRate','TruePositive','TrueNegetive','MissRate','FalsePositive'})
    legend('Left','Right')
    
end
