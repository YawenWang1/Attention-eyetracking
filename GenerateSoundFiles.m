%% Generate sounds with gap at different location
%  gain is from the sound_equalization 
%  change the number of gaps by adding and operator in for loop of k
% StimulusFolder = 'C:\Users\Yawen.Wang\Documents\Matalb Script\NBM_fMRI_Behavior\Stimulus\Sound\NGap';
% if exist(StimulusFolder,'dir')
%     cd(StimulusFolder)
% else
%     mkdir(StimulusFolder)
%     cd(StimulusFolder)
% end

 
fname = {'hslf','hfls','lfhs','lshf'};
for i = 1 : length(sound_signal)
%     tt =  split(sound_signal(i).ffm,'_');
%     tt = [{'2'}, {'4'}, {'2'}, {'4'}];
    freqmod_freq = [Slow_rate,Fast_rate,fliplr([Slow_rate,Fast_rate])];
%     tt = [{sound_signal.ffm(1)}, {sound_signal.ffm(3)}, {sound_signal.ffm(3)}, {sound_signal.ffm(1)}]; %% this matches the order we set in soundequalization     freqmod_freq = [Slow_rate(ii),Fast_rate(ii),fliplr([Slow_rate(ii),Fast_rate(ii)])];
%     freqmod_freq = [str2num(tt{1}) str2num(tt{2}) str2num(tt{2}) str2num(tt{1})];
    hslf = [sound_signal(i).value(1,:); sound_signal(i).value(2,:)*Gain(i).g];
    %     hslf = fastsmooth(hslf,10,3,1);
    filename = sprintf('%dNG_%s_%d.wav',4*(i-1)+i,'hslf',1);
    NGapSound{4*(i-1)+i}.filename = filename;
    NGapSound{4*(i-1)+i}.sound = hslf;
    NGapSound{4*(i-1)+i}.bpos = 0;
%     audiowrite(filename,hslf',Fs);
    
    
    hfls = [sound_signal(i).value(3,:); sound_signal(i).value(4,:)*Gain(i).g];
    %     hfls = fastsmooth(hfls,10,3,1);
    filename = sprintf('%dNG_%s_%d.wav',4*(i-1)+i+1,'hfls',2);
    NGapSound{4*(i-1)+i+1}.filename = filename;
    NGapSound{4*(i-1)+i+1}.sound = hfls;
    NGapSound{4*(i-1)+i+1}.bpos = 0;
    
%     audiowrite(filename,hfls',Fs);
    
    lfhs = [sound_signal(i).value(2,:)*Gain(i).g; sound_signal(i).value(1,:)];
    %     lfhs = fastsmooth(lfhs,10,3,1);
    filename = sprintf('%dNG_%s_%d.wav',4*(i-1)+i+2,'lfhs',3);
    NGapSound{4*(i-1)+i+2}.filename = filename;
    NGapSound{4*(i-1)+i+2}.sound = lfhs;
    NGapSound{4*(i-1)+i+2}.bpos = 0;
    
%     audiowrite(filename,lfhs',Fs);
    
    lshf = [sound_signal(i).value(4,:)*Gain(i).g; sound_signal(i).value(3,:)];
    %     lshf = fastsmooth(hslf,10,3,1);
    filename = sprintf('%dNG_%s_%d.wav',4*(i-1)+i+3,'lshf',4);
    NGapSound{4*(i-1)+i+3}.filename = filename;
    NGapSound{4*(i-1)+i+3}.sound = lshf;
    NGapSound{4*(i-1)+i+3}.bpos = 0;
    
%     audiowrite(filename,lshf',Fs);
end

% GStimulusFolder = 'C:\Users\Yawen.Wang\Documents\Matalb Script\NBM_fMRI_Behavior\Stimulus\Sound\Gap';
% if exist(GStimulusFolder,'dir')
%     cd(GStimulusFolder)
% else
%     mkdir(GStimulusFolder)
%     cd(GStimulusFolder)
% end
% 


for i = 1 : length(sound_signal)
%     tt = [{'2'}, {'4'}, {'2'}, {'4'}];
    freqmod_freq = [Slow_rate(ii),Fast_rate(ii),fliplr([Slow_rate(ii),Fast_rate(ii)])];
%     tt = [{sound_signal.ffm(1)}, {sound_signal.ffm(3)}, {sound_signal.ffm(3)}, {sound_signal.ffm(1)}]; %% this matches the order we set in soundequalization     freqmod_freq = [Slow_rate(ii),Fast_rate(ii),fliplr([Slow_rate(ii),Fast_rate(ii)])];
%     freqmod_freq = [str2num(tt{1}) str2num(tt{2}) str2num(tt{2}) str2num(tt{1})];

%     tt =  split(sound_signal(i).ffm,'_');
%     freqmod_freq = [str2num(tt{1}) str2num(tt{2}) str2num(tt{2}) str2num(tt{1})];
    temp_sound = sound_signal(i).value;
   
    for j = 1 : size(temp_sound,1)
        nBeeps = unique(sound_signal(i).bid{j});
        change_time = nBeeps(ceil(length(nBeeps)/2)-1:ceil(length(nBeeps)/2)+1);
    
%         if length(nBeeps) < 6
%             change_time = nBeeps(ceil(length(nBeeps)/2)-1:ceil(length(nBeeps)/2)+1);
%         else
%             change_time = nBeeps(ceil(length(nBeeps)/2)-2:ceil(length(nBeeps)/2)+1);
%         end
    
        for k = 1 : length(change_time)
            if j == 2 || j == 4 %% the high frequency sounds are placed in index 1 and 3          freqs = repmat([High_Freq(ii),Low_Freq(jj)],1,2);
                %% low frequency was at 2 and 4
                s = temp_sound(j,:)* Gain(i).g ;
                sound(1,:) = temp_sound(j-1,:);
                if j == 2
                     filename = sprintf('%dG_%s_%s_%d.wav',length(change_time)*(j-1)+k,'hslf','lf',change_time(k));
                     a = min(sound_signal(i).ct{i,j}(find(sound_signal(i).bid{i,j}==change_time(k))));
                end
                if j == 4              
                    filename = sprintf('%dG_%s_%s_%d.wav',length(change_time)*(j-1)+k,'hslf','ls',change_time(k));
                    a = min(sound_signal(i).ct{i,j}(find(sound_signal(i).bid{i,j}==change_time(k))));
                end
                    
            end

            if j == 1 || j ==3 
                s = temp_sound(j,:);
                sound(1,:) = temp_sound(j+1,:);
                if j == 1
                    filename = sprintf('%dG_%s_%s_%d.wav',length(change_time)*(j-1)+k,'lfhs','hs',change_time(k));
                    a = min(sound_signal(i).ct{i,j}(find(sound_signal(i).bid{i,j}==change_time(k))));
                end
                
                if j == 3
                    filename = sprintf('%dG_%s_%s_%d.wav',length(change_time)*(j-1)+k,'lshf','hf',change_time(k));
                    a = min(sound_signal(i).ct{i,j}(find(sound_signal(i).bid{i,j}==change_time(k))));
                end

            end
            vchangetime(4*length(change_time)*(i-1)+length(change_time)*(j-1)+k) = a;
            idx=(find(sound_signal(i).bid{j}==change_time(k)));
            idx=sound_signal(i).ct{j}(idx);
            idx=[idx(1)-50:idx(end)+50];
            s(idx)=0;
            
            %             s(sound_signal(i).ct{i,j}(find(sound_signal(i).bid{1,j}==change_time(k))))=0;
            
            %s(sound_signal(i).ct{j}([idx(1) idx(end)]))=2;
            %s(sound_signal(i).ct{j}([idx(1)-1 idx(end)+1]))=3;
            
            sound(2,:) = s; %%%%%%%% the gaps are only played on one speaker!
%             p = audioplayer(sound,Fs);
%             play(p)
%             pause
%             
            
%                          audiowrite(filename,sound',Fs);
            GapSound{4*length(change_time)*(i-1)+length(change_time)*(j-1)+k}.filename = filename;
            GapSound{4*length(change_time)*(i-1)+length(change_time)*(j-1)+k}.sound    = sound;
            GapSound{4*length(change_time)*(i-1)+length(change_time)*(j-1)+k}.bpos    = change_time(k);
            
            %             s=fastsmooth(s,10,3,1);
        end
    end
end

% % (length(NGapSound)+length(GapSound))
for i = 1 : (length(NGapSound)+length(GapSound))
    if i <= 4
        NAndG(i) = NGapSound(i);
    else
        NAndG(i) = GapSound(i-4);
    end
%     subplot(5,4,i)
%     plot(NAndG{i}.sound','.')
end

% % % % % 
% for i = 1 : length(NGapSound)
%     p = audioplayer(NGapSound{i}.sound,Fs);
%     
%     play(p);
%     pause
% 
%     NGapSound{i}.filename
% %     plot(NAndG{i}.sound')
%   
% end 
% % corr = 0;
%gaya makes all correct response at 3 and 6 

% for i =  10 : 12%length(GapSound)
%     p = audioplayer(GapSound{i}.sound,Fs);
%     play(p);
%     perform(i) = getkey;
%   
% %       pause
%     
%     if perform(i) == 104
%         if find((hgap-4) == i)
%             detectcorr(i) = 1;
%             corr = corr + 1;
%           
%             disp('Correct, you detect the gap in high frequency');
%         end
%     end
%     
%     
%     if perform(i) == 108
%         if find((lgap-4) == i)
%             detectcorr(i) = 1;
%             corr = corr + 1;
%             disp('Correct, you detect the gap in low frequency');
%         end
%     end       
%     GapSound{i}.filename
% end
% 
% 

% for i = 5 : length(NAndG)
%     p = audioplayer(NAndG{i}.sound,Fs);
%     
%     play(p);
%     NAndG{i}.filename
%     plot(NAndG{i}.sound')
%     pause
% end
% 
% for i = 1 : length(temp_sound)
%     p = audioplayer(temp_sound(i,:),Fs)
%     play(p)
%     pause
%     i
% end

% for i = 1 : length(sound_signal)
%     ttt = sound_signal(1).value;
%     for j = 1 : size(ttt,1)
%         p = audioplayer(ttt(j,:),Fs);
%         play(p)
%         pause
%     end
% end