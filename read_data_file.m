
%%tmp =  1  cnd(trial,1),
%        2  block,cnd(trial,2),
%        3  cnd(trial,3),
%        4  cnd(trial,4),
%        5  cnd(trial,5),
%        6  cnd(trial,6),
%        7  blockstartvbl,
%        8  tstarttime, 
%        9  tstarttime +0.5*ifi, 
%        10  tstarttime + (FDurFrames-0.5)*ifi, 
%        11  tstarttime + (FDurFrames+19-0.5)*ifi,
%        12  NAndG{cnd(trial,4)}.bpos,
%        13  vctpos,
%        14  corrans, 
%        15  pressedKey,
%        16  presstime,
%        17  (GetSecs-tstarttime)]
%% what are the last two numbers?

subject_number=unique(results(:,1))
blocks_numbers=unique(results(:,2));
blocks_code=unique(results(:,3));
trial_numbers=unique(results(:,4));
sound_gap=unique(results(:,5)) %% is this?
slow_sound_gap=unique(results(:,6)) %% is this?
fast_sound_gap=unique(results(:,7)) %% is this? blockstartvbl
tst_starttime=unique(results(:,8)) %% is this? tstarttime
tst_startframe=unique(results(:,9)) %% is this? tstarttime
odd_startframe=unique(results(:,10)) %% is this? tstarttime


results(find(results(:,3)==2 & results(:,16)==1),:)

results(find(results(:,3)==2 & results(:,16)==0),:)








