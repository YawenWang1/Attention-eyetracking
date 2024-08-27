%% Settings for VisuoAudio experiment
%% set experiment paradigmConditions
if ~exist('subjnum')
    subjnum='0';
end
NumberOfBlock = 8;
BlockCondition = {'Attend to Low Frequency',
    'Attend to High Frequency',
    'Attend to Left',
    'Attend to Right'};
% blockorder = Shuffle(repmat([1 2],1,2));
blockorder = Shuffle(repmat([1 2 3 4],1,2));
%% Response keys (optional; for no subject response use empty list)
responsekeys = {'space','ESCAPE','q'};
% responsekeys = {'UpArrow','DownArrow','ESCAPE'};
cue = {'L','H','<+<','>+>'};
%% Set timing for the eperiment
SDur  = 1.0; % duration of the whole stimulus
CDur = 2.6;
% CDur = 1.3;
% FDur  = 0.3;
% respwindow = 2; %gaya 3
% ngap = [1:1:4];lgap = [5:8];hgap = [9:12];attl = [5:8]; attr = [9:12];

%catch trials & non-catch trials
RatiOfCatch = 0.2;
BsCondition = repmat([eye(4);zeros(1,4)],20,1);
nTrials = size(BsCondition,1);
nCaTrials   = RatiOfCatch * nTrials;
phases = repmat(floor(linspace(0,180,4)),1,25)';
Orientations = repmat([floor(linspace(45,135,4));floor(linspace(135,45,4))]',25,1);
BsCondition = [BsCondition phases Orientations];

lsoundid = [repmat([1;2;3],6,1);[2;3]];
hsoundid = [repmat([4;5;6],6,1);[5;6]];
locofchange = [repmat([2;3;4],6,1);[3;4]];
iti = [0.5*ones(50,1);ones(50,1)];
iti = iti(randperm(nTrials));
for i = 1: length(blockorder)
    tempcond = [str2num(subjnum)*ones(size(BsCondition,1),1) blockorder(i)*ones(size(BsCondition,1),1) [1:nTrials]' BsCondition(randperm(size(BsCondition,1)),:)];
    tempcond(find(tempcond(:,5)),5) = lsoundid(randperm(nCaTrials));
    tempcond(find(tempcond(:,4)),4) = hsoundid(randperm(nCaTrials));
    tempcond(find(tempcond(:,6)),6) = locofchange(randperm(nCaTrials));
    tempcond(find(tempcond(:,7)),7) = locofchange(randperm(nCaTrials));

    Conditions(:,:,i)     = tempcond;
end

for i = 1 : length(blockorder)
    AllConditions((i-1)*nTrials+1:i*nTrials,:) = Conditions(:,:,i);
    
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
contrast = 20;
% Aspect ratio width vs. height:
aspectratio = 1.0;


bend = 0;
bdur = 0;
expdur = 0;
bstart = 0;

