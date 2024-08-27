%% Settings for VisuoAudio experiment
%% set experiment paradigm
if ~exist('subjnum')
    subjnum='0';
end
NumberOfBlock = 8;
BlockCondition = {'Attend to High Frequency',
    'Attend to Low Frequency',
    'Attend to Left',
    'Attend to Right'};
% blockorder = Shuffle(repmat([1 2 3 4],1,2));
%% Response keys (optional; for no subject response use empty list)
responsekeys = {'space','ESCAPE','q'};
% responsekeys = {'UpArrow','DownArrow','ESCAPE'};
cue = {'H','L','<','>',};
%% Set timing for the eperiment
SDur  = 1.2; % duration of the whole stimulus
FDur  = 0.3; CDur = 1.5;
respwindow = 5;
% ngap = [1:1:4];hgap = [5:1:8 13:1:16];lgap = [9:1:12 17:1:20];attl = [1:1:4]; attr = [5:1:8];
ngap = [1:1:4];lgap = [8:1:10 14:1:16];hgap = [5:1:7 11:1:13];attl = [1:1:3]; attr = [4:1:6];
vchg = [1:1:6];
% a = 0;v =1;
%%catch trials & non-catch trials
BsCondition = repmat([eye(2)],18,1);
RatiOfCatch = 0.5 
nTrials = length(BsCondition); nCaTrials   = RatiOfCatch * nTrials;
ISI = round((rand(1,nTrials) + 1));
if a ==1
    blockorder = Shuffle([1,2]);
    BsCondition(find(BsCondition(:,1)),1) = repmat(hgap',nCaTrials/length(hgap),1);
    BsCondition(find(BsCondition(:,1)==0),1) = [repmat(ngap',3,1);repmat(lgap',1,1)];
    BsCondition(find(BsCondition(:,2)),2) = repmat(lgap',nCaTrials/length(lgap),1);
    BsCondition(find(BsCondition(:,2)==0),2) = [repmat(ngap',3,1);repmat(lgap',1,1)];
end
if v == 1
    blockorder = Shuffle([3,4]);
    BsCondition(find(BsCondition(:,1)),1) = repmat(attl',nCaTrials/length(attl),1);
    BsCondition(find(BsCondition(:,2)),2) = repmat(attr',nCaTrials/length(attr),1);

end
%% settings for visual stimulus
ngabors = 2;
% Phase of underlying sine grating in degrees:
phase = 0;
% Spatial constant of the exponential "hull"
sc = 10.0;
% Frequency of sine grating:
sf = .05;
% Contrast of grating:
contrast = 10;
% Aspect ratio width vs. height:
aspectratio = 1.0;
si = 32; gaborw = 4*si+1; gaborh = 4*si + 1;

%% Set parameters for the auditory stimulus
T=SDur;Fs=44100;ts=1/Fs; t=0:ts:T;
change_times=[500 500];
% High_Freq = [2500 2000 1000];Low_Freq  = [200 400 500 600];
High_Freq = [1500];Low_Freq  = [400];
% Slow_rate = [2 3 4]; Fast_rate = Slow_rate * 2;
Slow_rate = [3]; Fast_rate = Slow_rate * 2;

for i = 1 : length(blockorder)
        bs = BsCondition;
        Conditions(:,:,i) = [str2num(subjnum)*ones(size(bs,1),1) blockorder(i)*ones(size(bs,1),1) [1:nTrials]' bs(randperm(size(bs,1)),:) ];
end
for i = 1 : length(blockorder)
        bs = BsCondition;
        Conditions01((i-1)*nTrials+1:(i-1)*nTrials+nTrials,:) = [str2num(subjnum)*ones(size(bs,1),1) blockorder(i)*ones(size(bs,1),1) [1:nTrials]' bs(randperm(size(bs,1)),:) ];
end




