% Fit FA model to each session and identify the optimal latent 
% dimensionality using cross-validation

%% Paths and specs

run ../addpaths

[datapath, statspath, figurepath] = addpaths;

% analyse correct trials
datafolder = 'correct';

monkeys    = {'Wa','Sa'};

% fit FA only to control (no-uStim) data
onlycontrol = true;

savefig     = false;

% pick the bin size carefully! too small and the covariance won't be
% properly stimated, so the latents will be dominated by single units.
% <50ms might be problematic, depending on the firing rate.
binWidth = 50; %ms

% fit FA models with different latent dimensionalities
nLatents = 1:15;

% trial period to analyse (pre and post uStim)
trialperiod = [-300,800]; %ms

if onlycontrol; condlab = 'nouStim'; else condlab = 'all'; end

optdim = [];

for m = 1:length(monkeys)

monkey = monkeys{m};

[filenames,~,~,stimorder] = datafiles(monkey);


for f=1:length(filenames)
    
filename = filenames{f};
fprintf('loading file %s \n',filename)

statsfile = [filename '_FAmodel_' condlab '.mat'];

if ~exist([statspath 'subspaces/FA_model_fits/' statsfile],'file')
    

load([datapath datafolder '/' filename],'spiketrain','params')

% change order of uStim conditions to have control condition first
spiketrain = spiketrain(stimorder,:,:,:);
params.cond_uStimChan = params.cond_uStimChan(stimorder);

if onlycontrol; spiketrain = spiketrain(1,:,:,:);end

switch monkey
    case 'Sa'
        % in monkey Sa, padd short trials with nans
        load([datapath datafolder '/' filename],'behavior')
        behavior = behavior(stimorder,:,:,:);
        if onlycontrol; behavior = behavior(1,:,:,:);end
        spiketrain = paddtrialswithnans(spiketrain,1,behavior);
end

% remove shorted channels
chidx  = selectchannels(monkey,datapath,filename);

% compute binned spike counts (pre-uStim, post-uStim and all times)
[spikeTrainsBin01, spikeTrainsBin02, spikeTrainsBin0] = ...
    binspiketrains(spiketrain,chidx,binWidth,trialperiod);

y  = cat(4,spikeTrainsBin01, spikeTrainsBin02); % uStim period excluded
y0 = spikeTrainsBin0;                           % uStim period included

% re-format data for plotting functions
y  = permute(y,[1 4 5 2 3]); % (numch,time,trials,uS,ang)
y0 = permute(y0,[1 4 5 2 3]);

dim     = size(y); dim0 = size(y0);
spikes  = nan([dim(1:3) prod(dim(4:5))]); %(numch,time,trials,prod(cond))
spikes0 = nan([dim0(1:3) prod(dim0(4:5))]);

c = 1;
for c1 = 1:dim(4)
    for c2 = 1:dim(5)
        spikes(:,:,:,c)  = y(:,:,:,c1,c2);
        spikes0(:,:,:,c) = y0(:,:,:,c1,c2);
        c = c+1;
    end
end
y  = spikes;
y0 = spikes0;

% save data and specs
FAdata          = y;
FAdata0         = y0;
FAinfo.params   = params;
FAinfo.datafile = filename;
FAinfo.binWidth = binWidth;

% fit FA model to data with uStimp excluded due to uStim artifacts
dims = size(y);
y    = reshape(y,dims(1),[]);

idx = isnan(sum(y,1));
y   = y(:,~idx);


%% Fit FA model

disp('fitting FA model...')

diags = crossvalidate_fa(y,'zDimList',nLatents);

% Identify optimal latent dimensionality
LL = [diags.sumLL]; % cross-validated log likelihood
PE = [diags.sumPE]; % cross-validated prediction error
[~,maxLL] = max(LL);
[~,minPE] = min(PE);
CVstats.nLatents  = nLatents;
CVstats.LL        = LL;
CVstats.PE        = PE;
CVstats.bestDimLL = nLatents(maxLL);
CVstats.bestDimPE = nLatents(minPE);
FAparams          = [diags.estParams];

save([statspath 'subspaces/FA_model_fits/' statsfile], ...
    'FAdata','FAdata0','FAinfo','FAparams','CVstats')

else
    disp('stats already computed')
    load([statspath 'subspaces/FA_model_fits/' statsfile], ...
    'CVstats')
    optdim  = [optdim CVstats.bestDimLL];
end

end
end

disp('done')


%%  Plot optimal dimensionality

h1 = figure; hold on
xx=1; yy=1;
pos = get(h1,'position');
set(h1,'position',[pos(1:2) pos(3)*xx pos(4)*yy])

numdims = mean(optdim);

histogram(optdim)
yLim = ylim;
plot([numdims numdims],[yLim(1) yLim(2)],'k','LineWidth',3)
xlabel('dimensionality')
ylabel('sessions')
title('FA optimal dimensionality')

fprintf('Mean optimal dimensionality = %.0f \n', numdims)

if savefig
    saveas(h1,[figurepath 'subspaces/' 'FAoptdim_' condlab],'svg')
end

