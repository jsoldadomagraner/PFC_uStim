% Fit Naive Bayes decoder (Gaussian or Poisson) to spike count data

%% Paths and specs

run ../addpaths

[datapath, statspath, ~] = addpaths;

% analyse correct trials
datafolder = 'correct';

monkeys    = {'Wa','Sa'};

% naive bayes model (Gaussian or Poisson)
NBmethod  = 'Poisson';

% cross-validation method (LOOCV or Kfold)
CVmethod  = 'LOOCV';

% CV folds (for CVmethod = Kfold)
folds    = 10;

% trial period to analyse (all trial)
trialperiod = [-800,800]; %ms

% bin spike counts
binWidth    = 50; %ms

%% Fit decoders to spiking data from each session

for m=1:length(monkeys)
    
monkey = monkeys{m}; 

[filenames,~,~,stimorder] = datafiles(monkey);


for f=1:length(filenames)
    
filename = filenames{f};

fprintf('loading file %s \n',filename)

statsfile = [filename '_' NBmethod(1) 'NBdecoder.mat'];


if ~exist([statspath 'decoder/' statsfile],'file')

load([datapath datafolder '/' filename],'spiketrain')
    
switch monkey
    case 'Sa'
        load([datapath datafolder '/' filename],'behavior')
        spiketrain = paddtrialswithnans(spiketrain,1,behavior);
end

% change order of uStim conditions to have control condition first
spiketrain = spiketrain(stimorder,:,:,:);

% remove shorted channels
chidx  = selectchannels(monkey,datapath,filename);

% Format data (pre-uStim and post-uStim binned spikecounts)
[spikeTrainsBin01, spikeTrainsBin02, ~] = ...
    binspiketrains(spiketrain,chidx,binWidth,trialperiod);

binspikecounts = cat(4,spikeTrainsBin01,spikeTrainsBin02);
[N,S,D,T,~]    = size(binspikecounts);

% Fit multi-class NB classifiers
% fit a classifier independently at each time point and for each uStim
% condition
Y     = cell(S,T);
Yhat  = cell(S,T);
Yhat0 = cell(S,T);

prob  = cell(T,1);
prob0 = cell(T,1);

accuracy  = cell(T,S);
accuracy0 = cell(T,S);

for t=1:T

% format data and create labels
% binspikecounts (chan,uStim,angle,time,trial)    
for i=1:S
    X{i} = []; Y{i,t} = [];
    for j=1:D
        angledat = reshape(binspikecounts(:,i,j,t,:),N,[]);
        nanidx   = ~isnan(sum(angledat,1));
        angledat = angledat(:,nanidx);
        X{i}   = [X{i} ; angledat'];
        Y{i,t} = [Y{i,t} ; ones(size(angledat,2),1)*j];
    end
end

for i=1:S
    disp('fitting Naive Bayes classifier')
    [Yhat{i,t},Error,prob{t,i}] = NaiveBayesClassifier(X{i}',Y{i,t}, ...
        NBmethod,CVmethod,folds,binWidth);
    
    disp('fitting Naive Bayes classifier control')
    idx = randperm(size(X{i},1)); % shuffled data
    [Yhat0{i,t},genError0,prob0{t,i}] = NaiveBayesClassifier(X{i}(idx,:)',Y{i,t}, ...
        NBmethod,CVmethod,folds,binWidth);
        
    Error = Error*100;
    accuracy{t,i} = 100-Error;
    
    genError0 = genError0*100;
    accuracy0{t,i} = 100-genError0;
end

end

save([statspath 'decoder/' statsfile], 'accuracy','accuracy0','folds',...
    'Y','Yhat','Yhat0','prob','prob0')

else
    disp('decoder already fit')
end

end
end