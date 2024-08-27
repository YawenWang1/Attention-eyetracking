%% settings for the staircase
MaxTrials = 120;
MaxReverals = 14;
IgnoreReversals = 4;
% Define the type of staircase
UpNum = 1;    % number of incorrect answer to go one step up
DownNum = 4;  % number of correct answers to go one step down
% Ratio of Up and Down stepsize
ratio = 0.8415;

R = [];
%% Beginning
% close all; clear all;clc;
KbName('UnifyKeyNames');
subjno=input('Please input date:', 's');
subjnum = '1';

% subjnum = input('Please input subject number:', 's');
subjname=input('Please input subject name:', 's');
a = 1; v = 0;
% settingsfortraining
Beeps
SoundEqualization_1
GenerateSoundFiles
filename = sprintf('%s_%s_%s',subjnum,subjname,subjno);
% data_directory = 'D:\Yawen\VBehavior\data\NBM_fMRI_Behavior\Aonly';
data_directory= 'C:\Users\Yawen.Wang\Documents\Matalb Script\NBM_fMRI_Behavior\data';
% data_directory=['D:\Yawen\VBehavior\data\NBM_fMRI_Behavior'];
% data_directory=['D:\Yawen\VBehavior\data\NBM_fMRI_Behavior']; % set a new directory
if ~exist(data_directory)
    mkdir(data_directory);
end

%% staircase
InitGapDur = 0.01;
StepSize   = 0.001;
MinGapDur  = 0.001;



%% Initial sound  stimuli with gap
for i = 1 : 2: length(NGap)
    gdur = round(Fs*InitGapDur);
    
    for j = idx
        %         jidx = idx(j);
        temp = EqSound(i).sound;
        temp1 = EqSound(i+1).sound;
        
        temp(1,(j-1)*BeepLength+loc:(j-1)*BeepLength+loc+gdur) = 0.3*max(temp(1,(j-1)*BeepLength+1:j*BeepLength));
        temp1(1,(j)*BeepLength+loc:(j)*BeepLength+loc+gdur) = 0.3*max(temp1(1,j*BeepLength+1:(j+1)*BeepLength));
        Hgap(3*(i-1)+j).sound =  temp;
        Hgap(3*(i-1)+j+1).sound =  temp1;
        
        
        
        temp2 = EqSound(i).sound;
        temp3 = EqSound(i+1).sound;
        temp2(2,j*BeepLength+loc:j*BeepLength+loc+gdur) = 0.3*max(temp2(2,j*BeepLength+1:(j+1)*BeepLength));
        temp3(2,(j-1)*BeepLength+loc:(j-1)*BeepLength+loc+gdur) = 0.3*max(temp3(2,(j-1)*BeepLength+1:(j)*BeepLength));
        Lgap(3*(i-1)+j).sound =  temp3;
        Lgap(3*(i-1)+j+1).sound =  temp2;

    end
end


run_experiment=1;
%% create log file 
% colheaders = {'subId','subName','blockNum','blockorder','cue','trialNum','catchTrial','soundfile','ISI','VctPostion','pressedkey','rt','rttoexpstart'}

if run_experiment
    
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
        [w,screenrect]=PsychImaging('OpenWindow',screennumber, 128,[0 0 400 600]);
        Screen('BlendFunction',w,GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
        Screen('Preference', 'DefaultFontSize', 35);
        Screen('Preference', 'DefaultFontStyle',1);
        Screen('Flip', w);
        ifi = Screen('GetFlipInterval',w);
%         Screen('BlendFunction', w, GL_ONE, GL_ONE);
        % get the center position of the screen
        [cx,cy] = RectCenter(screenrect);
%         textRect = [0; 0; 2*gaborw; 2*gaborh;];
%         dstRects(:,1) = CenterRectOnPoint(textRect',cx-34*4,cy);%1deg = 33.6 pix
%         dstRects(:,2) = CenterRectOnPoint(textRect',cx+34*4,cy);
        % for play auditory stimulus
        InitializePsychSound;
        pahandle = PsychPortAudio('Open', [], [], 0, Fs, 1);        
        keysOfInterest = zeros(1,256);
        keysOfInterest(KbName({'5%','ESCAPE',responsekeys})) = 1;
        KbQueueCreate;
        PsychHID('KbQueueCreate',-1,keysOfInterest);
        PsychHID('KbQueueStart',-1);

        
        %% fixation cross before start
        DrawFormattedText(w, 'Press Space Key To Begin', 'center', 'center', [255,255,255]);
        expstart=Screen('Flip', w);
        KbStrokeWait;
        
        
        %% experimental loops
        SDurFrames = round(SDur/ifi);CDurFrames = round(CDur/ifi); FDurFrames = round(FDur/ifi);respFrames = round(3/ifi);
        Vct = round([vchangetime(1:4) vchangetime(5:8)]./length(t)*SDurFrames) + 19;
        
        trialcounter = 0;
        for soundsti = 1 : length(HF)*length(LF)
            
        for  block  = 1 :2 %: size(Conditions,3)
            if block = 1
                tempsound = Hgap((soundsti-1)+1:(soundsti-1)+6);% initial gap with maximum gap size
            else
                tempsound = Lgap((soundsti-1)+1:(soundsti-1)+6);
            end
            
            
            revrs = []; gapdur = [];corrans = [];count=0;UpFlag=0; DownFlag=0;
            revrs_trial = [];
            trial=0;  
            blocktext = BlockCondition{block};
            DrawFormattedText(w, blocktext, 'center', 'center', 255);
            blockstartvbl = Screen('Flip',w);
            
            
            
            while length(revrs) <= MaxReverals  || trial <= MaxTrials 
                trialcounter = trialcounter + 1;
                trial = trial + 1;
                PsychHID('KbQueueStart',-1);
                
                index = randi(6,1);
                ps = sum(tempsound{index}.sound);
                tstarttime=blockstartvbl +round(ISI(trial)/ifi)*ifi
                DrawFormattedText(w, '+', 'center','center',255);
                vbl = Screen('Flip', w,tstarttime +0.5*ifi);
                
                % Now we present the isi interval with fixation point minus one frame because we presented the fixation point once already when getting a time stamp
                for frame = 1:FDurFrames - 1
                    % Draw the fixation point
                    DrawFormattedText(w, '+', 'center','center',128);
                    % Flip to the screen
                    vbl = Screen('Flip', w, vbl + 0.5 * ifi);
                end
                         
                PsychPortAudio('FillBuffer', pahandle, ps); % loads data into buffer
                PsychPortAudio('Start', pahandle,1,inf); %starts sound at infinity
                PsychPortAudio('RescheduleStart', pahandle, vbl++19*ifi+0.5*ifi, 0) %reschedules startime to
                for frame = 1 : CDurFrames
                    
                    Cue = cue{cnd(trial,2)};
                    DrawFormattedText(w, Cue, 'center','center',255);                    
                    vbl = Screen('Flip',w,vbl+0.5*ifi);
                end
                % record response
                    respStart = GetSecs;
                    
                    for frame  = 1 : respFrames
                        DrawFormattedText(w, '+', 'center','center',255);
                        vbl = Screen('Flip',w,vbl+0.5*ifi);
                    end
                    
                    %                         [secs, keyCode, deltaSecs] = KbStrokeWait(-1);
                    [KeyIsDown,firstPress] = PsychHID('KbQueueCheck',-1); %collect keyboard events in KbQueueStart was invoked
                if KeyIsDown
%                     if ~isempty(ismember(KbName(find(firstPress)),'space'))
                      if strcmp(KbName(find(firstPress)),'space')
                          exp_term=0;
                          
                          corrans(trial) = 1;
                          
                          pressedKey = 32;
                          presstime = firstPress(pressedKey) - respStart;
                          log{trialcounter}.allpressedkeyname = KbName(find(firstPress));
                          log{trialcounter}.allpressedkey = find(firstPress);
                          log{trialcounter}.allpressedkeyrt = firstPress(pressedKey) - respStart;

                       
                    end
                    if strcmp(KbName(find(firstPress)),'ESCAPE')
                        exp_term=1;
                        if exp_term
%                             
                            log{trialcounter}.subjnum = str2num(subjnum);
                            log{trialcounter.currblock = block;
                            log{trialcounter}.currtrial = trial;
                            log{trialcounter}.blockstartime = blockstartvbl-expstart;
                            log{trialcounter}.trialstartime = tstarttime - blockstartvbl;
                            log{trialcounter}.fixationstart = tstarttime - blockstartvbl + 0.5*ifi;
                            log{trialcounter}.cuestart = tstarttime - blockstartvbl +(FDurFrames-0.5)*ifi;
                            log{trialcounter}.stimulustart = tstarttime - blockstartvbl +(FDurFrames +19 -0.5)*ifi;
                            log{trialcounter}.beepos = index;
                            log{trialcounter}.correctans = corrans;
                            log{trialcounter}.key = pressedKey;
                            log{trialcounter}.rt = presstime;
                            log{trialcounter}.dur = (GetSecs-tstarttime);
                            log{trialcounter}.allpressedkeyname = KbName(find(firstPress));
                            log{trialcounter}.allpressedkey = find(firstPress);
                            log{trialcounter}.allpressedkeyrt = firstPress(pressedKey) - respStart;
                            tmp = [str2num(subjnum),block,trial,blockstartvbl-expstart,tstarttime - blockstartvbl, 0.5*ifi, (FDurFrames-0.5)*ifi, (FDurFrames+19-0.5)*ifi,index,corrans,pressedKey,presstime,(GetSecs-tstarttime)];
                            %                                     tmp = [cnd(trial,1),block,cnd(trial,2),cnd(trial,3),cnd(trial,4),cnd(trial,5),cnd(trial,6),blockstartvbl-expstart,tstarttime, 0.5*ifi, (FDurFrames-0.5)*ifi, (FDurFrames+19-0.5)*ifi,NAndG{cnd(trial,4)}.bpos,vctpos,corrans, pressedKey,presstime,(GetSecs-tstarttime)];
                            results(trialcounter) = tmp;
                            cd(data_directory)
                            save([filename,'atesting','.mat'],'results','log');
                            Priority(0);
                            ShowCursor;
                            Screen('CloseAll');
                            return;
                        end
                    end
                else
                    pressedKey = 0;
                    presstime = 0;
                    break_fixation = 1;
                    log{trialcounter}.allpressedkeyname = '0';
                    log{trialcounter}.allpressedkey = 0;
                    log{trialcounter}.allpressedkeyrt = 0;

                end
                
                gapsizes(1) = InitGapDur;
                
                if corrans(trial)
                count = count + 1; DownFlag = 1;
                gapsizes(trial+1) = gapsizes(trial);
                if count >= DownNum
                    count = 0; 
                    gapsizes(trial+1) = -StepSize;                 
                    if gapsizes(trial+1) >= InitGapDur
                        gapsizes(trial+1) = InitGapDur*0.999;
                    end
                    if UpFlag == 1
                        revrs = [revrs gapsizes(trial)];
                        revrs_trial = [revrs_trial trial];
                        
                    end
                    UpFlag = 0;
                end
                
            else
                count = 0;
                UpFlag = 1;
                gapsizes(trial+1) = gapsizes(trial) + StepSize;
                if gapsizes(trial+1) <= MinGapDur
                    gapsizes(trial+1) = MinGapDur;
                end
                if DownFlag == 1
                    revrs = [revrs gapsizes(trial)];
                    revrs_trial = [revrs_trial trial];    
                end
                DownFlag = 0;         
                end
                
                
                for a = 1 : length(tempsound)
                    newgapdur = round(Fs*gapsizes(trial+1));
                    if block == 1
                        
                        tempsound(a).sound(1,(a-1)*BeepLength+loc: (a-1)*BeepLength+loc+newgapdur) = 0.3*max(tempsound(a).sound(1,(a-1)*BeepLength+1:a*BeepLength));
                        
                    else
                        tempsound(a).sound(2,(a-1)*BeepLength+loc: (a-1)*BeepLength+loc+newgapdur) = 0.3*max(tempsound(a).sound(2,(a-1)*BeepLength+1:a*BeepLength));
                        
                    end
                end
                
                
                
                

                
                
                
                
                
                
                
                    log{trialcounter}.subjnum = cnd(trial,1);
                    log{trialcounter}.currblock = block;
                    log{trialcounter}= blockstartvbl-expstart;
                    log{trialcounter}.trialstartime = tstarttime - blockstartvbl;
                    log{trialcounter}.fixationstart = tstarttime - blockstartvbl + 0.5*ifi;
                    log{trialcounter}.cuestart = tstarttime - blockstartvbl +(FDurFrames-0.5)*ifi;
                    log{trialcounter}.stimulustart = tstarttime - blockstartvbl +(FDurFrames +19 -0.5)*ifi;
                    log{trialcounter}.beepos = index;
                    log{trialcounter}.correctans = corrans;
                    log{trialcounter}.key = pressedKey;
                    log{trialcounter}.rt = presstime;
                    log{trialcounter}.dur = (GetSecs-tstarttime);
                    
                    tmp = [str2num(subjnum),block,trial,blockstartvbl-expstart,tstarttime - blockstartvbl, 0.5*ifi, (FDurFrames-0.5)*ifi, (FDurFrames+19-0.5)*ifi,index,corrans(trial),pressedKey,presstime,(GetSecs-tstarttime)];
                    results(nTrials*(block-1)+trial,1:length(tmp)) = tmp;
                    
            end
            R(2*(soundsti-1)+block).resvr = revrs;R(2*(soundsti-1)+block).revrstrial = revrs_trial;
            T(2*(soundsti-1)+block).threshold = mean(revrs(end-10:end));
            C(2*(soundsti-1)+block).correctans  = corrans[];
           
            DrawFormattedText(w, 'Have a Rest.Press any key to start', 'center','center',255);
            vbl = Screen('Flip', w);
            blockendtime = KbStrokeWait; blockdur = blockendtime - blockstartvbl;
            bend{block} = blockendtime; bdur{block} = blockdur;
            KbStrokeWait
            WaitSecs(0.2);
            cd(data_directory)
            save(sprintf('%s%s%s',filename,'Astaircase','.mat'),'results','log','R','T','C');

     
        end
        end
        DrawFormattedText(w, 'End of the experiments,Press any key to exit', 'center','center',255);
        endexp = Screen('Flip', w);
        KbStrokeWait;
        expdur = endexp - expstart;
        save(sprintf('%s%s%s%s',filename,'Astaircase','blockinfo','.mat'),'bend','bdur','expdur','');
        Screen('CloseAll');
    catch
        Screen('CloseAll');
end
end

