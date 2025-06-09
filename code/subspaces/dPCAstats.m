% Compute memory subspace stats (dPCA model)

%% Paths and specs

run ../addpaths

[datapath, statspath, ~] = addpaths;

% analyse correct trials
datafolder = 'correct';

monkeys    = {'Wa','Sa'};

% use mean optimal latent dimensionality of the FA model
numdims    = 4;

% compute stats for the first 2 latent dims for example 2D trajectories
stats2D    = false;

dimlab  = [num2str(numdims) 'Dmodel'];
if stats2D; dimlab = [dimlab '_2Dstats'];end

% fit dPCA model to each session using 5-fold cross-validation
folds = 5;

% trial period to analyse (pre and post uStim)
trialperiod = [-300,800]; %ms


disp('computing dPCA stats')

for m = 1:length(monkeys)

monkey   = monkeys{m};

statsfile = [monkey '_dPCAstats_' dimlab '.mat'];    
  
if ~exist([statspath 'subspaces/' statsfile],'file')

[filenames,~,~,stimorder] = datafiles(monkey);

for f=1:length(filenames)

filename = filenames{f};
fprintf('file %s \n',filename)

load([datapath datafolder '/' filename],'spiketrain','params')

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

firingRates01 = spikeTrainsBin01/(binWidth/1000);
firingRatesAverage01 = nanmean(firingRates01,5);

firingRates02 = spikeTrainsBin02/(binWidth/1000);
firingRatesAverage02 = nanmean(firingRates02,5);

firingRatesAverage = cat(4,firingRatesAverage01,firingRatesAverage02);
T = size(firingRatesAverage,4);

%% load FA model to compute variance and marginalized variance using dPCA functions

FApath   = [statspath 'subspaces/FA_model_fits/'];
FAmodel0 = [filename '_FAmodel_nouStim'];

% load FA model fitted to no uStim data
load([FApath FAmodel0],'FAparams')

L = FAparams(numdims).L; % non-orthonormalized loadings

% orthogonalize L
nLatents = size(L,2); 
nVars = size(L,1); 
minVl = min(nLatents, nVars); 
[U, Sigma, ~] = svd(L); 
% The SVD almost gives a unique decomposition of L, but we still need to
% account for the the sign on the columns of U.  
for lI = 1:minVl
    uColSign = sign(sum(U(:,lI))); 
   if uColSign < 0
        U(:,lI) = -1*U(:,lI); 
        Sigma(lI,lI) = -1*Sigma(lI,lI);
    end
end
L = U(:,1:minVl);

%Use uStim trials only (which were held-out in the estimated FA model) 
X     = firingRatesAverage(:,2:S,:,:); 
Xdims = size(X);
idx   = isnan(squeeze(sum(X,1:length(Xdims)-1)));
X(:,:,:,idx) = [];
XuStim = X;

% Compute variance in each dimension and marginalize var
% Define parameter grouping for dPCA function
% {'Stimulus', 'Angle', 'Condition-independent(Time)', 'S/A Interaction'}
combinedParams = {{1, [1 3]}, {2, [2 3]}, {3}, {[1 2], [1 2 3]}};

explVar_FA = dpca_explainedVariance(XuStim, L, L, ...
    'combinedParams', combinedParams);

FAvar{f} = explVar_FA;

%% Cross-validation for dPCA model estimation

cvpar = cvpartition(E,'KFold',folds);

trialreconstruct = nan(N,S-1,D,T,E,numdims);

for cv=1:folds
    
idx1 = cvpar.training(cv); % training indexes
idx2 = ~idx1;              % test indexes

% training data
firingRates01 = spikeTrainsBin01(:,:,:,:,idx1)/(binWidth/1000); 
firingRatesAverage01 = nanmean(firingRates01,5);

firingRates02 = spikeTrainsBin02(:,:,:,:,idx1)/(binWidth/1000);
firingRatesAverage02 = nanmean(firingRates02,5);

firingRatesAverage = cat(4,firingRatesAverage01,firingRatesAverage02);
T = size(firingRatesAverage,4);

firingRatesAverage_train = firingRatesAverage;
firingRates01_train      = firingRates01;
firingRates02_train      = firingRates02;

% test data
firingRates01_test = spikeTrainsBin01(:,:,:,:,idx2)/(binWidth/1000);
firingRates02_test = spikeTrainsBin02(:,:,:,:,idx2)/(binWidth/1000); 
firingRates_test   = cat(4,firingRates01_test,firingRates02_test);

%% Fit dPCA to mean firing rates from training data

% W is the decoder, V is the encoder (ordered by explained variance),
% whichMarg is an array that tells you which component comes from which
% marginalization

numComp = [numdims numdims numdims numdims];
 
[W,V,whichMarg] = dpca(firingRatesAverage_train, numComp, ...
    'combinedParams', combinedParams);


%% Reconstruct test trials

% Find memory dimensions (marginalization 2)
marg = 2;
mg   = find(whichMarg==marg);
mg   = mg(1:numdims);

Wmem = W(:,mg);
Vmem = V(:,mg);

% Reconstructed test trials based on dPCA model fitted to training data
testsz = size(firingRates_test);
for si=2:S %uStim
    for di=1:D % angle
        for ti=1:testsz(4)
            for i=1:size(Vmem,2)
                trialreconstruct(:,si-1,di,ti,idx2,i) = ...
                    Vmem(:,i)*Wmem(:,i)'*squeeze(firingRates_test(:,si,di,ti,:));
            end
        end
    end
end

%% Project training and test data onto the memory subspace

A = find(whichMarg==marg);

A = A(1:numdims);

la = length(A);
componentsToPlot = 1:length(whichMarg);

% training data
Yplotstd{1} = firingRates01_train;
Yplotstd{2} = firingRates02_train;

for d=1:2
    dataDim = size(Yplotstd{d});
    X = reshape(Yplotstd{d},dataDim(1),[])';
    Xcen = bsxfun(@minus, X, nanmean(X));
    Z = Xcen * W;
    Zfull = reshape(Z(:,componentsToPlot)', [length(componentsToPlot) dataDim(2:end)]);
    
    for i=1:la
        compstdadist_train{d,i} = Zfull(A(i),:,:,:,:); 
    end
end

% test data
Yplotstd{1} = firingRates01_test;
Yplotstd{2} = firingRates02_test;

for d=1:2
    dataDim = size(Yplotstd{d});
    X = reshape(Yplotstd{d},dataDim(1),[])';
    Xcen = bsxfun(@minus, X, nanmean(X));
    Z = Xcen * W;
    Zfull = reshape(Z(:,componentsToPlot)', [length(componentsToPlot) dataDim(2:end)]);
    
    for i=1:la
        compstdadist_test{d,i} = Zfull(A(i),:,:,:,:); 
        % gather all test trial projections for all cv folds here
        compstdadist{d,i}(:,:,:,:,idx2) = Zfull(A(i),:,:,:,:);
    end
end

% reformat data
% xmang {angle}(dims,time,trials,u))
xmang = cell(1,4);
for a=1:D
    i=1;
    temp = squeeze(cat(4,compstdadist_train{1,i}(1,:,a,:,:),compstdadist_train{2,i}(1,:,a,:,:)));
    xmang{a} = zeros([3,size(temp)]);
    xmang{a}(i,:,:,:) = temp;
    for i=2:la
        temp = squeeze(cat(4,compstdadist_train{1,i}(1,:,a,:,:),compstdadist_train{2,i}(1,:,a,:,:)));
        xmang{a}(i,:,:,:) = temp; %(dims,uStim,time,trials)
    end
xmang{a} = permute(xmang{a},[1 3 4 2]); %(dims,time,trials,uStim)
end
xmang_train = xmang;

xmang = cell(1,4);
for a=1:D
    i=1;
    temp = squeeze(cat(4,compstdadist_test{1,i}(1,:,a,:,:),compstdadist_test{2,i}(1,:,a,:,:)));
    xmang{a} = zeros([3,size(temp)]);
    xmang{a}(i,:,:,:) = temp;
    for i=2:la
        temp = squeeze(cat(4,compstdadist_test{1,i}(1,:,a,:,:),compstdadist_test{2,i}(1,:,a,:,:)));
        xmang{a}(i,:,:,:) = temp; %(dims,uStim,time,trials)
    end
xmang{a} = permute(xmang{a},[1 3 4 2]); %(dims,time,trials,uStim)
end
xmang_test = xmang;


%% Decoding in the memory subspace

% Based on distance of held-out uStim trials to the centroid of the 
% no uStim distributions for each target angle condition

% compute stats in all dimensions
dimstats = la;
%focus on the first dimensions only
if stats2D; dimstats = 2; end

[~,uStimp,~,tpre,~] = analysisperiod(trialperiod);
taupre = 1:tpre;
tpost  = tpre+1;
uStimp = nan(1,length(uStimp));

%xmang{ang}(ln,ti,tr,length(uStime);
tr  = size(xmang_test{1},3);
Tau = {1:tpre,tpost:T};
uStimptr = repmat(uStimp,tr,1);
for a=1:D
    for co=2:S
        for d=1:dimstats
            for p=1:2
                tau = Tau{p};
                controlpos = squeeze(nanmean(xmang_train{a}(d,tau,:,1),3));
                uStimposdist    = squeeze(xmang_test{a}(d,tau,:,co));
                mcontprepost{p} = controlpos;
                if tr>1; uStimposdist = uStimposdist';end
                uStrprepost{p}  = uStimposdist;
            end
        meancont(d,a,:) = cat(2,mcontprepost{1},uStimp,mcontprepost{2});
        uStrials(d,a,co,:,:) = cat(2,uStrprepost{1},uStimptr,uStrprepost{2});
        end
    end
end

angmeancont = nanmean(meancont,3); %(dims,ang,t)
for a=1:D
    for co=2:S
        dist = zeros(D,tr,size(uStrials,5));
        for tri=1:tr
        for t=1:size(uStrials,5)
            trial = squeeze(uStrials(:,a,co,tri,t));
            %compute distance wrt control vectors
            for ang=1:D
                dist(ang,tri,t) = norm(angmeancont(:,ang)-trial);
            end
        end
        end
        % assign class based on minimum distance to centroids
        [~,distmin] = min(dist,[],1);
        distmin = squeeze(distmin);
        idx = distmin == a;
        dPCAdec{f}(a,co-1,idx2,:) = idx;
    end
end

W_cv(:,:,cv) = W;
V_cv(:,:,cv) = V;
whichMarg_cv(:,cv) = whichMarg;

clear meancont uStrials 

end

VdPCA{f} = V_cv;
WdPCA{f} = W_cv;
whichMargdPCA{f} = whichMarg_cv;
CVpartition{f}   = cvpar;

% Statistics such as variance of the held-out data explained on each
% CV-fold are too noisy to compute independently on each fold, since only a
% few trials are held out per fold.
% Instead, compute statistics on the projections of all held-out data
% points, noting that this will also be a noisy estimate, given that the
% projections are made on slightly different subspaces (found based on the
% training data for each CV-fold)

%% Percent variance explained, cross-validated (based on held out data reconstructions for all trials)

XuStim_hat = squeeze(nanmean(trialreconstruct,5));

% Compute variance in each dimension and marginalize var
explVar_dPCA = dpca_explainedVariance_CVtest(XuStim, XuStim_hat, ...
    'combinedParams', combinedParams);

dPCAvar{f}   = explVar_dPCA;


%% Gather projections of all test trials onto the memory subspace

xmang = cell(1,4);
for a=1:D
    i=1;
    temp = squeeze(cat(4,compstdadist{1,i}(1,:,a,:,:),compstdadist{2,i}(1,:,a,:,:)));
    xmang{a} = zeros([3,size(temp)]);
    xmang{a}(i,:,:,:) = temp;
    for i=2:la
        temp = squeeze(cat(4,compstdadist{1,i}(1,:,a,:,:),compstdadist{2,i}(1,:,a,:,:)));
        xmang{a}(i,:,:,:) = temp; %(dims,uStim,time,trials)
    end
xmang{a} = permute(xmang{a},[1 3 4 2]); %(dims,time,trials,uStim)
end

dims  = 1:dimstats;

%% Are uStim and no uStim distributions statistically different?

% run test and compute p-values
pv    = multivartest(xmang,dims);

sigthreshold = 0.05;

%pv(Ang,uStim,ldims,time)
hy    = pv<sigthreshold;
hpre  = hy(:,:,:,1:tpre);
hpost = hy(:,:,:,tpost:end);
dimsh = size(hy);
dimsh(end) = length(uStimp);
uStimph = nan(dimsh);
hy = cat(4,hpre,uStimph,hpost);

dPCAst(f,:,:,:,:) = hy;

%% Distance between uStim and no uStim trajectories in no uStim CI units

z = 1.96; % confidence interval

for a=1:D
    for co=2:S
        for d=1:dimstats
            for p=1:2
                controlpos = squeeze(nanmean(compstdadist{p,d}(1,1,a,:,:),5));
                uStimpos   = squeeze(nanmean(compstdadist{p,d}(1,co,a,:,:),5));
                stdcontrol = squeeze(nanstd(compstdadist{p,d}(1,1,a,:,:),0,5));
                trials     = sum(squeeze(~isnan(compstdadist{p,d}(1,1,a,:,:))),2);
                CIcontrol  = z*stdcontrol./sqrt(trials);
                mdistprepost{p} = (uStimpos - controlpos)./CIcontrol;
            end
        normdimstime(d,:)  = cat(2,mdistprepost{1}',uStimp,mdistprepost{2}');
        end
        for t=1:size(normdimstime,2)
            normtime(t) = norm(normdimstime(:,t));
        end
        dPCAcist(f,a,co-1,:)  = normtime;
    end
end

clear compstdadist W_cv V_cv whichMarg_cv hy_cv

end
    save([statspath 'subspaces/' statsfile],'dPCAst','dPCAcist','dPCAdec', ...
        'VdPCA','WdPCA','whichMargdPCA','CVpartition','dPCAvar', 'FAvar')
else
    disp('stats already computed')
end

clear dPCAst dPCAcist dPCAdec VdPCA WdPCA whichMargdPCA CVpartition dPCAvar FAvar

end
