%% staircase for the task of frequency change detection
% trial 1 to 10, as the filler, to aviod the pattern of key pressing.
% factors for the first 11 is given, same with all the participants.
% the 11th trial is the real start of the staircase, easiest one to detect
% Inserted the correct rejection trial
%% Beginning
close all; clear all;clc;
Sound_Freq_Change
SE_Frequency
KbName('UnifyKeyNames');
subjno=input('Please input date:', 's');
subjnum = '9';
% subjnum = input('Please input subject number:', 's');
subjname=input('Please input subject name:', 's');
filename = sprintf('%s_%s_%s',subjnum,subjname,subjno);
% data_directory = 'D:\Yawen\VBehavior\data\NBM_fMRI_Behavior\Aonly';
% data_directory= 'C:\Users\Yawen.Wang\Documents\Matalb Script\NBM_fMRI_Behavior\data';
data_directory=['D:\Yawen\VBehavior\data\NBM_fMRI_Behavior_data\Staircase']; % set a new directory
if ~exist(data_directory)
    mkdir(data_directory);
end

% load('1_Mark2_20180611Astaircase.mat');
%% settings for the staircase
MaxTrials = 120;
MaxReverals = 14;
IgnoreReversals = 4;
% Define the type of staircase
UpNum = 1;    % number of incorrect answer to go one step up
DownNum = 4;  % number of correct answers to go one step down
% Ratio of Up and Down stepsize
ratio = 0.8415;

R = [];C = [];

%% staircase

InitFactor  = .5;
InitChgHighFreq =  round(f(2)*10.^(InitFactor*0.301));
InitChgLowFreq = round(f(1)*10.^(-InitFactor*0.301));
nf(1) = InitChgLowFreq;
Chgpos     = 2;

%% Initial sound for frequency change
temp_f = [InitChgLowFreq,InitChgHighFreq]; % low frequency at 400 Hz and High frequency at 3000Hz

erb_n = 24.7*(4.37*temp_f/1000+1);           % standard equation: ERB(f) = 0.108f + 24.7 (f belongs to [100Hz to 10KHz])
b = 0.00437;
erbs_n = 21.4 * log10(1 + b*temp_f);            % ERB scale ERBS(f) = 21.4*log10(0.00437*f + 1);


%%%%%% use these values if you use the pass band filter
erbs_low_n = erbs_n - (no_erbs/2);
erbs_high_n = erbs_n + (no_erbs/2);
delta_erb_n=(erbs_high_n-erbs_low_n)./(no_freq-1);

% predefine vectors

erbs_left_n = zeros(length(temp_f),10);
erbs_right_n = zeros(length(temp_f),10);
f_left_n = zeros(length(temp_f),10);
f_right_n = zeros(length(temp_f),10);

for i = 1 : length(temp_f)
    for k = 1:(no_freq-1)/2
        
        erbs_left_n(i,k)=erbs_low_n(i)+(k-1)*delta_erb_n(i);
        erbs_right_n(i,k)=erbs_high_n(i)-(k-1)*delta_erb_n(i);
        f_left_n(i,k)=(10^(erbs_left_n(i,k)/21.4) - 1)./b;
        f_right_n(i,k)=(10^(erbs_right_n(i,k)/21.4) - 1)./b;
    end
end

f_vector_n = [f_left_n,temp_f', fliplr(f_right_n)];   %%% this is the final vector of frequency to be used
% f_vector_n1 = temp_f/1.2* ones(1,size(f_vector_n,2));
%% Generate the narrowband sound
for i = 1 : length(temp_f)
    chgfreqsound = zeros(1,nTot);
    for j = 1:no_freq
        
            amp = 1;
            freq = f_vector_n(i,j);                     % use frequency from ERB vector
%             freq = f_vector_n1(i,j);   
            phase = rand(1,1)*2*pi;                  % generate random phase for each frequency
            m(1:nRamps) = amp*RampUp .*sin(2*pi*freq * totT(1:nRamps));
            m(nRamps+1:nRamps+nSamples) = amp * sin(2*pi*freq*totT(nRamps+1:nRamps+nSamples));
            m(nRamps+nSamples+1:nTot) = amp * RampDown .* sin(2*pi*freq *totT(nRamps+nSamples+1:nTot));
        
%             m(1:nRamps) = amp*RampUp .*sin(2*pi*freq * totT(1:nRamps)+phase);
%             m(nRamps+1:nRamps+nSamples) = amp * sin(2*pi*freq*totT(nRamps+1:nRamps+nSamples)+phase);
%             m(nRamps+nSamples+1:nTot) = amp * RampDown .* sin(2*pi*freq *totT(nRamps+nSamples+1:nTot)+phase);
%         
            chgfreqsound = chgfreqsound + m;             % sum all sinusoids in signal
        
             % sum all sinusoids in signal
        
    end
    if i == 1
        chglowfreqsound =  AmpMod.*chgfreqsound;
    else
        chghighfreqsound = AmpMod.*chgfreqsound;
    end
        
    
end
chglowfreqsound = chglowfreqsound/norm(chglowfreqsound);
chghighfreqsound = chghighfreqsound/norm(chghighfreqsound);

%% Set timing for the eperiment
SDur  = 1.0; % duration of the whole stimulus
FDur  = 0.3; CDur = 1.3;
BlockCondition = {'Attend to Low Frequency',
    'Attend to High Frequency',
    'Attend to Left',
    'Attend to Right'};

cue = {'L','H','< <','> >',};
%% Initial sound  stimuli



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
        [w,screenrect]=PsychImaging('OpenWindow',screennumber, 128,[]);
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
        keysOfInterest(KbName({'5%','ESCAPE','space'})) = 1;
        KbQueueCreate(-1,keysOfInterest);
        KbQueueStart(-1);
% %         PsychHID('KbQueueCreate',-1,keysOfInterest);
%         PsychHID('KbQueueStart',-1);

        
        %% fixation cross before start
        DrawFormattedText(w, 'Press Space Key To Begin', 'center', 'center', [255,255,255]);
        expstart=Screen('Flip', w);
        KbStrokeWait;
        KbQueueFlush(-1);
%         PsychHID('KbQueueFlush',-1);

        
        %% experimental loops
        SDurFrames = round(SDur/ifi);CDurFrames = round(CDur/ifi); FDurFrames = round(FDur/ifi);respFrames = round(2/ifi);
         
        trialcounter = 0;
        for soundsti = 1 : length(EqSound)
            soundst = EqSound(soundsti).sound;
        for  block  =  1:2%size(Conditions,3)
            % Initial sound at each block
            
            
            if block == 1
               
                soundst(1,Chgpos(1)*nTot+1:(Chgpos(1)+1)*nTot) = chglowfreqsound;
            else
                soundst(2,Chgpos(1)*nTot+1:(Chgpos(1)+1)*nTot) = chghighfreqsound;
            end
            
            forder = zeros(1,11);factors= zeros(1,11);TrialIdx=ones(1,11);
            nf = zeros(1,11);
            forder(1:11) = [1 1 1 0 1 0 1 0 1 0 1];
            presetfactors = [0.5, 0.00,0.0156,0.0075,  0.125, 0.00, 0.5,0.00,0.0313,0.00,0.5];
            factors(1:11) = forder.*presetfactors;
            nf
            
%             TrialIdx(1:11) = 1;
            
            revrs = []; gapdur = [];corrans = [];count=0;UpFlag=0; DownFlag=0;
            revrs_trial = [];
            trial=0;  
            blocktext = BlockCondition{block};
            DrawFormattedText(w, blocktext, 'center', 'center', 255);
            blockstartvbl = Screen('Flip',w);
            
            
            
            while length(revrs) <= MaxReverals  && trial <= MaxTrials 
                trialcounter = trialcounter + 1;
                trial = trial + 1;
                KbQueueStart(-1);
%                 PsychHID('KbQueueStart',-1);
             
                Chgpos(trial+1) = randi([2,4],1);

                ps = mean(soundst);
                tstarttime=blockstartvbl +randi(2,1)*ifi;
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
                    
                    Cue = cue{block};
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
%                     [KeyIsDown,firstPress] = PsychHID('KbQueueCheck',-1); %collect keyboard events in KbQueueStart was invoked
                    [pressed, firstPress, firstRelease, lastPress, lastRelease] = PsychHID('KbQueueCheck',-1);
           
                    if pressed
                        %                     if ~isempty(ismember(KbName(find(firstPress)),'space'))
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
                        if strcmp(KbName(find(lastPress)),'ESCAPE')
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
                                log(trialcounter).upf = nf(trial);
                                log(trialcounter).dur = (GetSecs-tstarttime);
                                log(trialcounter).forder = forder(trial);
                                log(trialcounter).factors = factors(trial);
                                log(trialcounter).beeps = Chgpos(trial);
                                log(trialcounter).trialidx = TrialIdx(trial+1);
                                log(trialcounter).dur = (GetSecs-tstarttime);

                                log(trialcounter).allpressedkeyname = KbName(find(lastPress));
                                log(trialcounter).allpressedkey = find(lastPress);
                                log(trialcounter).allpressedkeyrt = lastPress(pressedKey) - respStart;
                                
                                tmp = [str2num(subjnum),block,trial,forder(trial),factors(trial),TrialIdx(trial+1),blockstartvbl-expstart,tstarttime - blockstartvbl, 0.5*ifi, (FDurFrames-0.5)*ifi, (FDurFrames+19-0.5)*ifi,Chgpos(trial),nf(trial),corrans(trial),pressedKey,presstime,(GetSecs-tstarttime)];
                                %                                     tmp = [cnd(trial,1),block,cnd(trial,2),cnd(trial,3),cnd(trial,4),cnd(trial,5),cnd(trial,6),blockstartvbl-expstart,tstarttime, 0.5*ifi, (FDurFrames-0.5)*ifi, (FDurFrames+19-0.5)*ifi,NAndG{cnd(trial,4)}.bpos,vctpos,corrans, pressedKey,presstime,(GetSecs-tstarttime)];
                                results(trialcounter) = tmp;
                                cd(data_directory)
                                save(sprintf('%s%s%s',filename,'Astaircase','.mat'),'results','log','R','Threshold','C','EqSound','Performance','factors','Trialidx');
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
                            % no frequency change at next trial
                            % factor will be same with the previous
                            % trial
                            nf(trial+1) = round(f(block)*10.^((-1)^(block)*0*0.301));
                            forder(trial+1) = 0;
                            factors(trial+1) = factors(trial);
%                             if corrans(trial)
%                                 count = count + 1; DownFlag = 1;
%                                 factors(trial+1) = factors(trial);
%                                 if count >= DownNum
%                                     count = 0;
%                                     factors(trial+1) = factors(trial)/2;
%                                 end
%                                 
%                             else
%                                 count = 0;UpFlag = 1;
%                                 factors(trial+1) = factors(trial)*2;
%                             end
%                             
                        else
                            
                            % frequency change
                            forder(trial+1) = 1;
                            if corrans(trial)
                                count = count + 1; DownFlag = 1;
                                factors(trial+1) = factors(trial);
                                % calculate new frequency
                                nf(trial+1) =  round(f(block)*10.^((-1)^(block)*factors(trial+1)*0.301));
                                if count >= DownNum
                                    count = 0;
%                                     factors(trial+1) = ratio*factors(trial)/2;
                                    factors(trial+1) = factors(trial)/2;
                                    nf(trial+1) =  round(f(block)*10.^((-1)^(block)*factors(trial+1)*0.301));
                                    %                     factors(trial+1) = factors(trial)/StepSize(length(revrs)+1)*ratio;
                                    if factors(trial+1) <= 0.01556
                                        factors(trial+1) =0.01556;
                                    end
                                    if UpFlag == 1
                                        revrs = [revrs factors(trial)];
                                        revrs_trial = [revrs_trial trial];
                                        
                                    end
                                    UpFlag = 0;
                                end
                                
                            else
                                count = 0;
                                UpFlag = 1;
                                factors(trial+1) = factors(trial)*2;
                                if TrialIdx(trial+1) <=0.2
                                    nf(trial+1) = nf(trial);
                                else
                                    nf(trial+1) =  round(f(block)*10.^((-1)^(block)*factors(trial+1)*0.301));
                                end
                                
                                %                 factors(trial+1) = factors(trial) + StepSize(length(revrs)+1);
                                if factors(trial+1) >= InitFactor
                                    factors(trial+1) = InitFactor;
                                end
                                if DownFlag == 1
                                    revrs = [revrs factors(trial)];
                                    revrs_trial = [revrs_trial trial];
                                end
                                DownFlag = 0;
                            end
                        end
                    else
                        nf(trial+1) = round(f(block)*10.^((-1)^(block)*factors(trial+1)*0.301));
                    end
                    
                    % generate new frequency;
                    erb_up = 24.7*(4.37*nf(trial+1)/1000+1);           % standard equation: ERB(f) = 0.108f + 24.7 (f belongs to [100Hz to 10KHz])
                    erbs_up = 21.4 * log10(1 + b*nf(trial+1));            % ERB scale ERBS(f) = 21.4*log10(0.00437*f + 1);
                    
                    
                    %%%%%% use these values if you use the pass band filter
                    erbs_low_up = erbs_up - (no_erbs/2);
                    erbs_high_up = erbs_up + (no_erbs/2);
                    delta_erb_up=(erbs_high_up-erbs_low_up)./(no_freq-1);
                    
                    % predefine vectors
                    
                    erbs_left_up = zeros(1,10);
                    erbs_right_up = zeros(1,10);
                    f_left_up = zeros(1,10);
                    f_right_up = zeros(1,10);
                    
                    
                    for k = 1:(no_freq-1)/2
                        
                        erbs_left_up(k)=erbs_low_up+(k-1)*delta_erb_up;
                        erbs_right_up(i,k)=erbs_high_up-(k-1)*delta_erb_up;
                        f_left_up(k)=(10^(erbs_left_up(k)/21.4) - 1)./b;
                        f_right_up(k)=(10^(erbs_right_up(k)/21.4) - 1)./b;
                    end
                    
                    f_vector_up = [f_left_up,nf(trial+1)', fliplr(f_right_up)];   %%% this is the final vector of frequency to be used
                    % f_vector_n1 = temp_f/1.2* ones(1,size(f_vector_n,2));
                    %% Generate the narrowband sound
                    
                    upfreqsound = zeros(1,nTot);
                    for j = 1:no_freq
                        
                        amp = 1;
                        freq = f_vector_up(j);                     % use frequency from ERB vector
%                         phase = rand(1,1)*2*pi;                  % generate random phase for each frequency
                        
                        m(1:nRamps) = amp*RampUp .*sin(2*pi*freq * totT(1:nRamps));
                        m(nRamps+1:nRamps+nSamples) = amp * sin(2*pi*freq*totT(nRamps+1:nRamps+nSamples));
                        m(nRamps+nSamples+1:nTot) = amp * RampDown .* sin(2*pi*freq *totT(nRamps+nSamples+1:nTot));
                        
                        %                     m(1:nRamps) = amp*RampUp .*sin(2*pi*freq * totT(1:nRamps)+phase);
                        %                     m(nRamps+1:nRamps+nSamples) = amp * sin(2*pi*freq*totT(nRamps+1:nRamps+nSamples)+phase);
                        %                     m(nRamps+nSamples+1:nTot) = amp * RampDown .* sin(2*pi*freq *totT(nRamps+nSamples+1:nTot)+phase);
                        
                        upfreqsound = upfreqsound + m;             % sum all sinusoids in signal
                        
                    end
%                     upfreqsound = upfreqsound;
                    
                    upfreqsound = AmpMod.*upfreqsound;
                    upfreqsound = upfreqsound/norm(upfreqsound);
                    soundst = EqSound(soundsti).sound;
                    if block == 2
                        % because of sound equalization, multiply Gain to
                        % the highfrequency. H-f is at 2nd row
                        soundst(block,(Chgpos(trial+1)-1)*nTot+1:Chgpos(trial+1)*nTot) = Gain(1).g*upfreqsound;
                    else
                        
                        soundst(block,(Chgpos(trial+1)-1)*nTot+1:Chgpos(trial+1)*nTot) = upfreqsound;
                    end

                    %record the log
                    log(trialcounter).subjnum = str2num(subjnum);
                    log(trialcounter).currblock = block;
                    log(trialcounter).blocktime= blockstartvbl-expstart;
                    log(trialcounter).trialstartime = tstarttime - blockstartvbl;
                    log(trialcounter).fixationstart = tstarttime - blockstartvbl + 0.5*ifi;
                    log(trialcounter).cuestart = tstarttime - blockstartvbl +(FDurFrames-0.5)*ifi;
                    log(trialcounter).stimulustart = tstarttime - blockstartvbl +(FDurFrames +19 -0.5)*ifi;
                    log(trialcounter).correctans = corrans(trial);
                    log(trialcounter).key = pressedKey;
                    log(trialcounter).rt = presstime;
                    log(trialcounter).upf = nf(trial);
                    log(trialcounter).forder = forder(trial);
                    log(trialcounter).factors = factors(trial);
                    log(trialcounter).beeps = Chgpos(trial);
                    log(trialcounter).trialidx = TrialIdx(trial+1);
                    log(trialcounter).dur = (GetSecs-tstarttime);
                    
                    tmp = [str2num(subjnum),block,trial,forder(trial),factors(trial),TrialIdx(trial+1),blockstartvbl-expstart,tstarttime - blockstartvbl, 0.5*ifi, (FDurFrames-0.5)*ifi, (FDurFrames+19-0.5)*ifi,Chgpos(trial),nf(trial),corrans(trial),pressedKey,presstime,(GetSecs-tstarttime)];
                    results(trialcounter,1:length(tmp)) = tmp;
                    KbQueueFlush(-1);

%                     PsychHID('KbQueueFlush',-1);
                     
            end
            Frequency(2*(soundsti-1)+block).freq = nf;
            Factors(2*(soundsti-1)+block).factors = factors;
            R(2*(soundsti-1)+block).resvr = revrs;R(2*(soundsti-1)+block).revrstrial = revrs_trial;
            Threshold(2*(soundsti-1)+block).threshold_factor = mean(revrs(end-10:end));
            Threshold(2*(soundsti-1)+block).threshold_freq = round(f(block)*10.^((-1)^(block)*mean(revrs(end-10:end))*0.301));
            C(2*(soundsti-1)+block).correctans  = corrans;
            
            DrawFormattedText(w, 'Have a Rest.Press any key to start', 'center','center',255);
            blockendtime = Screen('Flip', w);blockdur = blockendtime - blockstartvbl;
            KbStrokeWait;
            bend(block) = blockendtime; bdur(block) = blockdur;
            KbStrokeWait;
            WaitSecs(0.2);
            cd(data_directory)
            save(sprintf('%s%s%s',filename,'Astaircase','.mat'),'results','log','R','Threshold','C','EqSound','Performance','factors','TrialIdx','nf','forder');
            
            
            end
     
        
        end
        DrawFormattedText(w, 'End of the experiments,Press any key to exit', 'center','center',255);
        endexp = Screen('Flip', w);
        KbStrokeWait;
        expdur = endexp - expstart;
        cd(data_directory)
        save(sprintf('%s%s%s%s',filename,'Astaircase','blockinfo','.mat'),'bend','bdur','expdur');
        Screen('CloseAll');
    catch
        Screen('CloseAll');
end
end

plot(find(corrans==1),nf(find(corrans==1)),'g*'); hold on,plot(find(corrans==0),nf(find(corrans==0)),'r*');