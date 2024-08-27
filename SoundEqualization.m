% close all, clear all, close all force
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gain=[1  1  1 1];
soundcount = 0;
soundequalization=0;
for ii = 1 : length(High_Freq)
    freqmod_freq = [Slow_rate(ii),Fast_rate(ii),fliplr([Slow_rate(ii),Fast_rate(ii)])];
    ffmod = strcat(num2str(Slow_rate(ii)),'_',num2str(Fast_rate(ii)));
    for jj = 1 : length(Low_Freq)
        soundcount = soundcount + 1;
        
        Freqstr = strcat(num2str(High_Freq(ii)),'_',num2str(Low_Freq(jj)));
        freqs = repmat([High_Freq(ii),Low_Freq(jj)],1,2);
        for count=1:4
            freq=freqs(count); % FREQUENCY
            phases = zeros(1,length(t));
            W=freq.*(2*pi); %radians per sec
            freqmod=zeros(1,length(phases));
            
            freqmod_amp=freq./8;
            if freqmod_freq(1)==freqmod_freq(2)
                if count==2 | count==4
                    freqmod_orig=sin(2*pi*freqmod_freq(count)*(t+pi));
                else
                    freqmod_orig=cos(2*pi*freqmod_freq(count)*(t+pi));
                end
            else
                freqmod_orig=cos(2*pi*freqmod_freq(count)*t); %only run here.
            end
            freqmod_orig=round(freqmod_orig.*2)./3;
            freqmod=(freqmod_orig.*freqmod_amp).*(2*pi);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for time=2:T./ts %%% SIMULATION of the phase-oscillator   %%%%%%%%%%
                phases(time)= phases(time-1) + ts*(W)+freqmod(time).*ts ;
            end
            phase_trace=  (mod(phases',2*pi))';
            
            inst_freq_true=(1/ts)./ ((2*pi)./diff(phases )) ;
            powmod=ones(1,length(phase_trace)).*gain(count);
            
            center_times=find(freqmod_orig==0);
            
            beep_cnt=1;beep_idx=1;
            for j=2:length(center_times)
                if (center_times(j)-center_times(j-1))>1
                    beep_cnt=beep_cnt+1;
                end
                beep_idx(j)=beep_cnt;
            end
            all_beep_idx{count}=beep_idx;
            all_center_times{count} = center_times;
            %powmod(center_times(find(beep_idx==change_times(count))))=0;
            powmod(find(abs(freqmod_orig)>0.05))=0;
            oscillatory_sig=powmod.*cos(phase_trace);
            
            %             oscillatory_sig(center_times(find(beep_idx==3)))=0;
            oscillatory_sig=fastsmooth(oscillatory_sig,10,3,1);
            %     fn = 0.5*Fs;
            %     fd = 50;% freq-(freq/8);
            %     fu = 3000;%freq+(freq/8);
            %     [B, A] = fir1(1000,[fd fu]/fn);
            %     noise_limited = filter(B, A, randn(1,length(oscillatory_sig)));
            %     noise_limited(center_times(find(beep_idx==5)))=noise_limited(center_times(find(beep_idx==5))).*2;
            s_matrix(count,:) = oscillatory_sig;
            
            
        end
        sound_signal(soundcount).ffm = ffmod;
        sound_signal(soundcount).freq = Freqstr;
        sound_signal(soundcount).value=s_matrix;%+(noise_limited.*.2);
        sound_signal(soundcount).ct = all_center_times;
        sound_signal(soundcount).bid = all_beep_idx;
        %         close all
        %         plot(sound_signal(1:2,:)')
        %         
    end
end

% for i = 1 : size(sound_signal,1)
%     p = audioplayer(sound_signal(i,:)',Fs)
%     play(p)
%     pause(2)
% end

if soundequalization
    for i = 1 : length(sound_signal)
        performance=0;
        while abs(performance)<.9
            Hi_gain=[ones(1,5) zeros(1,5 )];
            Hi_gain=Hi_gain(randperm(length(Hi_gain)));
            for trial=1:length(Hi_gain)
                pause(0.1)
                order=ceil(rand*4 );
                if order==1
                    test_sound=([sound_signal(i).value(1,:)*(1-Hi_gain(trial)); sound_signal(i).value(2,:)*Hi_gain(trial); ]);
                elseif order==2
                    test_sound=([sound_signal(i).value(3,:)*(1-Hi_gain(trial)); sound_signal(i).value(4,:)*Hi_gain(trial)]);
                elseif order==3
                    test_sound=([sound_signal(i).value(2,:)*Hi_gain(trial);sound_signal(i).value(1,:)*(1-Hi_gain(trial))]);
                elseif order==4
                    test_sound=([sound_signal(i).value(4,:)*Hi_gain(trial); sound_signal(i).value(3,:)*(1-Hi_gain(trial))]);
                end
                %             test_sound(1,:)=fastsmooth(test_sound(1,:),10,3,1);
                %             test_sound(2,:)=fastsmooth(test_sound(2,:),10,3,1);

                p=audioplayer(mean(test_sound),Fs);
                play(p);

                response(trial)=getkey;
                if Hi_gain(trial)>0
                    if response(trial)==104; %%H
                        correct(trial)=-1;
                    elseif response(trial)==108;
                        correct(trial)=1;
                    end
                else
                    if response(trial)==104; %%H
                        correct(trial)=1;
                    elseif response(trial)==108;
                        correct(trial)=-1;
                    end
                end
            end
            performance=mean(correct)
        end
        Performance(i).p4m = performance
        direction=round(performance)
        Performance(i).drct = direction
    end




    for i = 1 : length(sound_signal)
%         Start_level=.2;
        Start_level=2;Hi_gain=Start_level;count=0;up_flag=0; down_flag=0; same_flag=0;
        reversal=[]; trial=0;response=[];reversal_trial=[];
        figure(i)
        numFig = sprintf('fig%s',num2str(i));
        while length(reversal)< 12%25
            trial=trial+1;
            pause(0.1)
            order=ceil(rand*4 );
            if order==1
                test_sound=([sound_signal(i).value(1,:); sound_signal(i).value(2,:)*Hi_gain(trial)]);
            elseif order==2
                test_sound=([sound_signal(i).value(3,:); sound_signal(i).value(4,:)*Hi_gain(trial)]);
            elseif order==3
                test_sound=([sound_signal(i).value(2,:)*Hi_gain(trial);sound_signal(i).value(1,:)]);
            elseif order==4
                test_sound=([sound_signal(i).value(4,:)*Hi_gain(trial); sound_signal(i).value(3,:)]);
            end
            %         test_sound(1,:)=fastsmooth(test_sound(1,:),10,3,1);
            %         test_sound(2,:)=fastsmooth(test_sound(2,:),10,3,1);

            p=audioplayer(mean(test_sound),Fs);
            play(p);

            response(trial)=getkey;
            if direction==1
                if response(trial)==104; %%H
                    correct(trial)=0;
                elseif response(trial)==108;
                    correct(trial)=1;
                end
            else
                if response(trial)==104; %%H
                    correct(trial)=1;
                elseif response(trial)==108;
                    correct(trial)=0;
                end
            end

            % correct(trial)=round(rand);

            step_size=1.5;
            if length(reversal)>1
                step_size=1.25;
            end
            if length(reversal)>3
                step_size=1.1;
            end
            if length(reversal)>6
                step_size=1.05;
            end
            if length(reversal)>10
                step_size=1.01;
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
            figure(i)
            Hi_gain(trial+1)=Hi_gain(trial);
            if correct(trial)>0
                plot(trial, Hi_gain(trial),'g*');hold on
                count=count+1;
                if count==1
                    count=0; Hi_gain(trial+1)=Hi_gain(trial)*1/step_size;
                    up_flag=1;
                    if down_flag==1
                        reversal=[reversal Hi_gain(trial)];
                        reversal_trial=[reversal_trial trial];
                        plot(trial, Hi_gain(trial),'ko');hold on
                    end
                    down_flag=0;
                    same_flag=0;
                end
            else
                plot(trial, Hi_gain(trial),'r*');hold on
                count=0; Hi_gain(trial+1)=Hi_gain(trial)*step_size;
                down_flag=1;
                if up_flag==1
                    reversal=[reversal Hi_gain(trial)];
                    reversal_trial=[reversal_trial trial];
                    plot(trial, Hi_gain(trial),'ko');hold on
                end
                up_flag=0;
                same_flag=0;
            end
            drawnow
        end
        Hi_gain_setting=median(reversal(10:end))
        plot([0 trial],[Hi_gain_setting Hi_gain_setting],'r-')

        Gain(i).g = Hi_gain_setting;
    end
else
    for i = 1 : length(sound_signal)
        Gain(i).g =  1.8390;%1.2121;%170.2038;%3.620%61.4633;%1.4316;%1.0694;%1.4667;%0.4563;
    end
end
    

