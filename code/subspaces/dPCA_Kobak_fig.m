% Plot dPCA model results for an example session 
% following Kobak et al. 2016 style

%% Paths and specs

monkey   = 'Wa';
session  = 8;

savefig  = false;

[datapath, statspath, figurepath] = addpaths;

% analyse correct trials
datafolder = 'correct';

% use mean optimal latent dimensionality of the FA model
numdims    = 4;

%% Load data

[filenames,~,~,stimorder] = datafiles(monkey);

filename = filenames{session};

fprintf('file %s \n',filename)

load([datapath datafolder '/' filename],'spiketrain','params')

% trial period to analyse (pre and post uStim)
trialperiod = [-300,800]; %ms

switch monkey
    case 'Sa'
        % in monkey Sa, padd short trials with nans
        load([datapath datafolder '/' filename],'behavior')
        spiketrain = paddtrialswithnans(spiketrain,1,behavior);
end

% change order of uStim conditions to have control condition first
spiketrain = spiketrain(stimorder,:,:,:);
params.cond_uStimChan = params.cond_uStimChan(stimorder);

% remove shorted channels
chidx    = selectchannels(monkey,datapath,filename);

binWidth = params.binSize_spikerate*1000;

% compute binned spike counts (pre-uStim and post-uStim)
[spikeTrainsBin01, spikeTrainsBin02, ~] = ...
    binspiketrains(spiketrain,chidx,binWidth,trialperiod);

[N,S,D,T01,E] = size(spikeTrainsBin01);
taupre        = 1:T01;

firingRates01 = spikeTrainsBin01/(binWidth/1000);
firingRatesAverage01 = nanmean(firingRates01,5);

firingRates02 = spikeTrainsBin02/(binWidth/1000);
firingRatesAverage02 = nanmean(firingRates02,5);

firingRatesAverage = cat(4,firingRatesAverage01,firingRatesAverage02);
T = size(firingRatesAverage,4);
time = 1:T+4; % add uStim period

%% Fit dPCA model to mean firing rates

% dPCA marginalizations
margNames = {'Stimulus', 'Angle', 'Condition-independent(Time)', 'S/A Interaction'};
combinedParams = {{1, [1 3]}, {2, [2 3]}, {3}, {[1 2], [1 2 3]}};
margColours    = [[0.9290 0.6940 0.1250]*256; 187 20 25; 100 100 100; 46 139 87]/256;
numComp   = [numdims numdims numdims numdims];

[W,V,whichMarg] = dpca(firingRatesAverage, numComp, ...
    'combinedParams', combinedParams);

explVar = dpca_explainedVariance(firingRatesAverage, W, V, ...
    'combinedParams', combinedParams);


%% Plot Kobak et al. dPCA figures

[h1, h2]=dpca_plot_uStim(firingRatesAverage, W, V, taupre, ... 
    @dpca_plot_default_uStim, ...
    'explainedVar', explVar, ...
    'marginalizationNames', margNames, ...
    'marginalizationColours', margColours, ...
    'whichMarg', whichMarg, ...
    'time', time);

if savefig
    saveas(h1,[figurepath 'subspaces/' 'dPCA_Kobak_fig1'],'svg')
    saveas(h2,[figurepath 'subspaces/' 'dPCA_Kobak_fig2'],'svg')
end
