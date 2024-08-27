%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sound_Freq_Change

soundequalization=0;


if soundequalization
    for i = 1 
        performance=0;correct = [];response = [];hl_order = [];
        while abs(performance)<.9
            for trial=1:10
                pause(0.1)
                order=ceil(rand*2);
                if order == 1
                    test_sound = 2*highsound;
                else
                    test_sound = 2*lowsound;
                end
                hl_order(trial) = order;
                p=audioplayer(test_sound,Fs);
                play(p);
                response(trial)=getkey;
                if order  == 1
                    if response(trial)==104; %%H
                        correct(trial)=1;
                    elseif response(trial)==108;
                        correct(trial)=-1;
                    end
                else
                    if response(trial)==104; %%H
                        correct(trial)=-1;
                    elseif response(trial)==108;
                        correct(trial)=1;
                    end
                end
            end
            performance=mean(correct)
        end
        Performance(i).p4m = performance;
        
        direction=round(performance);
        Performance(i).drct = direction;
        Performance(i).hlorder = hl_order;
        
    

        %% play the low frequency sound at 10db, and regulate the high frequency
        
        Start_level=5;High_Gain=Start_level;count=0;up_flag=0; down_flag=0; same_flag=0;
        reversal=[]; trial=0;response=[];reversal_trial=[];
        figure(i)
        numFig = sprintf('fig%s',num2str(i));
        while length(reversal)<= 14 && trial <=140%25
            trial=trial+1;
            pause(0.1)
       
%             test_sound = [10*log10(High_Gain(trial))*highsound;10*log10(10)*lowsound];
            test_sound = [High_Gain(trial)*highsound;2*lowsound];
            p=audioplayer(mean(test_sound),Fs);
            play(p);

            response(trial)=getkey;
            if direction==1  % direction will be 1 when the performance near to 1
                if response(trial)==104; %%H
                    correct(trial)=1;
                elseif response(trial)==108;
                    correct(trial)=0;
                end
            else
                if response(trial)==104; %%H
                    correct(trial)=0;
                elseif response(trial)==108;
                    correct(trial)=1;
                end
            end

            % correct(trial)=round(rand);

            step_size=1.5;
            if length(reversal)>1
                step_size=1.25;
            end
            if length(reversal)>3
                step_size=1.3;
            end
            if length(reversal)>6
                step_size=1.2;
            end
            if length(reversal)>10
                step_size=1.1;
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%             figure(i)
            High_Gain(trial+1)=High_Gain(trial);
            if correct(trial)>0
%                 plot(trial, High_Gain(trial),'g*');hold on
                count=count+1;
                if count==1
                    count=0; High_Gain(trial+1)=High_Gain(trial)/step_size;
                    up_flag=1;
                    if down_flag==1
                        reversal=[reversal High_Gain(trial)];
                        reversal_trial=[reversal_trial trial];
%                         plot(trial, High_Gain(trial),'ko');hold on
                    end
                    down_flag=0;
                    same_flag=0;
                end
            else
%                 plot(trial, High_Gain(trial),'r*');hold on
                count=0; High_Gain(trial+1)=High_Gain(trial)*step_size;
                down_flag=1;
                if up_flag==1
                    reversal=[reversal High_Gain(trial)];
                    reversal_trial=[reversal_trial trial];
%                     plot(trial,High_Gain(trial),'ko');hold on
                end
                up_flag=0;
                same_flag=0;
            end
            drawnow
        end
        Hi_gain_setting=median(reversal(4:end)); %6.9453 me
%         plot([0 trial],[Hi_gain_setting Hi_gain_setting],'r-')
        close all;
        Gain(i).g = Hi_gain_setting; %2.3203(me)
end
  %  end
else
    for i = 1
        Performance(i).p4m = 1;
        
        Performance(i).drct = 1;
        Performance(i).hlorder = [];
        Gain(i).g = 1.4815;%0.2534;%2.5278;%0.4633;%4.4071;%2.4151;%0.5551;%3.5256;%1.3058;%3.1869%0.6717%6.9453%2.3203;%1.8390;%1.2121;%170.2038;%3.620%61.4633;%1.4316;%1.0694;%1.4667;%0.4563;
    end
end

% for i = 1 : 2 : length(NGap)
%     Gain(i+1).g = Gain(i).g
% end



%%Generate the sound struct
for i  = 1
    EqSound(i).sound = [2*lowsound;highsound* Gain(i).g]; 
end

staircase = 1;
if staircase
    
    %Gaya
    temp_f = [333,4396];%[333,4075];%[287,4396];%[335,4067];%[311,4600];%[292,4136];%[311,4328];%[283,4336];%[367,4067];                         % low frequency at 400 Hz and High frequency at 3000Hz
    
    no_erbs = 4;                            % specify the bandwitdh of the noise (in ERBs units) (ues an even number)
    no_freq = 21;                           % specify the number of frequencies (use an odd number)
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
            
            chgfreqsound = chgfreqsound + m;             % sum all sinusoids in signal
            
            
            
            
            
            % sum all sinusoids in signal
            
        end
        if i == 1
            chglowfreqsound = AmpMod.*chgfreqsound;
            %                 LowSound = rawsound;
            
        else
            chghighfreqsound = AmpMod.*chgfreqsound;
            %         HighSound = rawsound;
        end
        
        
    end
    chglowfreqsound = chglowfreqsound/norm(chglowfreqsound);
    
    chghighfreqsound = chghighfreqsound/norm(chghighfreqsound);
    
    
   
    for i = 1:2
        for j = 1:3;
            if i == 1
                ttsound = EqSound(1).sound;
                ttsound(1,j*8820+1:(j+1)*8820) = 2*chglowfreqsound;
                NewSound((i-1)*3+j).sound = ttsound;
            else
                ttsound = EqSound(1).sound;
                ttsound(2,j*8820+1:(j+1)*8820) = chghighfreqsound;
                NewSound((i-1)*3+j).sound = ttsound;
                
            end
            
        end
    end
    
end

% gapduration = round(Fs*0.002)
% 
% i = 3;
% temp = EqSound(i).sound;
% 
% temp(1,(j-1)*BeepLength+loc:(j-1)*BeepLength+loc+gapduration/2) =  10*(1+cos(pi*(0:1:gapduration/2)/gapduration/2)/2).*temp(1,(j-1)*BeepLength+loc:(j-1)*BeepLength+loc+gapduration/2);
% temp(1,(j-1)*BeepLength+loc+gapduration/2+1:(j-1)*BeepLength+loc+gapduration) =  10*(1+cos(pi*(gapduration/2+1:1:gapduration)/gapduration/2)/2).*temp(1,(j-1)*BeepLength+loc+gapduration/2+1:(j-1)*BeepLength+loc+gapduration);
% p = audioplayer(temp,Fs);
% play(p)