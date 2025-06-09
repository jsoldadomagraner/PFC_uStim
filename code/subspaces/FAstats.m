% Compute dominant subspace stats (FA model)

%% Paths and specs

run ../addpaths

[datapath, statspath, ~] = addpaths;

% analyse correct trials
datafolder = 'correct';

monkeys    = {'Wa','Sa'};

numdims    = 4; % mean optimal latent dimensionality

% compute stats for the first 2 latent dims for example 2D trajectories
stats2D    = false;

dimlab  = [num2str(numdims) 'Dmodel'];
if stats2D; dimlab = [dimlab '_2Dstats'];end

smoothlndata = false; % smooth latents

%% Compute stats

disp('computing FA stats')

for m = 1:length(monkeys)

monkey    = monkeys{m};
statsfile = [monkey '_FAstats_' dimlab '.mat'];    
    
if ~exist([statspath 'subspaces/' statsfile],'file')
    
[filenames,~,~,~] = datafiles(monkey);

for i=1:length(filenames)
    
filename = filenames{i};

fprintf('file %s \n',filename)

FAmodel  = [filename '_FAmodel_all'];
FAmodel0 = [filename '_FAmodel_nouStim'];

FApath   = [statspath 'subspaces/FA_model_fits/'];

% load FA data from both uStim and no uStim conditions. 
% FAdata contains the data pre and post uStim
load([FApath FAmodel],'FAdata')
% load FA model fitted to no uStim data
load([FApath FAmodel0],'FAinfo','FAparams')

uStimA = FAinfo.params.cond_uStimAmp;
uStime = FAinfo.params.cond_uStimChan;
angle  = FAinfo.params.cond_targetAngle;

[ch,ti,tr,co] = size(FAdata);
y  = reshape(FAdata,ch,[]);

params.d   = FAparams(numdims).d;
params.c   = FAparams(numdims).L; % non-orthonormalized
params.psi = FAparams(numdims).Ph;

[postMeans, ~]  = fastfa_estep(y,FAparams(numdims));
[~,postMeans]   = orthonormalizeFAMdl(params,postMeans.mean');
postMeans       = postMeans'; %sorted by decreasing shared var

ln  = size(postMeans,1);
xm  = reshape(postMeans,ln,ti,tr,co);

trialperiod = [-300,800]; %ms
[~,uStimp,~,tpre,~] = analysisperiod(trialperiod);
taupre = 1:tpre;
tpost  = tpre+1;
uStimp = nan(1,length(uStimp));

if smoothlndata
    xmpre  = smoothdata(xm(:,taupre,:,:),2,'gaussian',6);
    xmpost = smoothdata(xm(:,taupre(end)+1:end,:,:),2,'gaussian',6);
else
    xmpre  = xm(:,taupre,:,:);
    xmpost = xm(:,taupre(end)+1:end,:,:);
end
xm  = cat(2,xmpre,xmpost);

tau     = 1:size(xm,2);
ti      = length(tau);
taupost = tau(taupre(end)+1:end);

c = 1;
xmang  = cell(1,4); xmang(:) = {zeros(ln,ti,tr,length(uStime))};
for c1 = 1:length(uStime)
    for c2 = 1:length(uStimA)
        for c3 = 1:length(angle)
            xmang{c3}(:,:,:,c1)  = xm(:,:,:,c);
            c = c+1;
        end
    end
end

% compute stats in all dimensions
dimstats = ln;

%focus on the first dimensions only
if stats2D; dimstats = 2; end

dims  = 1:dimstats;

%% Are uStim and no uStim distributions statistically different?

% Since the FA model is fit to no uStim data, the latent trajectories
% inferred using this model for uStim conditions are computed using 
% held-out data.

% run test and compute p-values
pv    = multivartest(xmang,dims);

sigthreshold = 0.05;

%pv(ang,uStim-1,ldims,time), ldims=1 for multivariate tests
hy    = pv<sigthreshold;
hpre  = hy(:,:,:,1:tpre);
hpost = hy(:,:,:,tpost:end);
dimsh = size(hy);
dimsh(end) = length(uStimp);
uStimph    = nan(dimsh);
hy = cat(4,hpre,uStimph,hpost);

FAst(i,:,:,:,:) = hy;

%% Distance between uStim and no uStim trajectories in no uStim CI units

z = 1.96; % confidence interval

%xmang{ang}(ln,ti,tr,length(uStime);
Tau = {1:tpre,tpost:ti};
for a=1:length(angle)
    for co=2:length(uStime)
        for d=1:dimstats
            for p=1:2
                tau = Tau{p};
                controlpos = squeeze(nanmean(xmang{a}(d,tau,:,1),3));
                uStimpos   = squeeze(nanmean(xmang{a}(d,tau,:,co),3));
                stdcontrol = squeeze(nanstd(xmang{a}(d,tau,:,1),0,3));
                trials     = sum(squeeze(~isnan(xmang{a}(d,tau,:,1))),2)';
                CIcontrol  = z*stdcontrol./sqrt(trials);
                mdistprepost{p} = (uStimpos - controlpos)./CIcontrol;
            end
        normdimstime(d,:) = cat(2,mdistprepost{1},uStimp,mdistprepost{2});
        end
        for t=1:size(normdimstime,2)
            normtime(t) = norm(normdimstime(:,t));
        end
        FAcist(i,a,co-1,:)  = normtime;
    end
end


%% Decoding in the dominant subspace

% Based on distance of held-out uStim trials to the centroid of the 
% no uStim distributions for each target angle condition

%xmang{ang}(ln,ti,tr,length(uStime);
Tau = {1:tpre,tpost:ti};
uStimptr = repmat(uStimp,tr,1);
for a=1:length(angle)
    for co=2:length(uStime)
        for d=1:dimstats
            for p=1:2
                tau = Tau{p};
                controlpos = squeeze(nanmean(xmang{a}(d,tau,:,1),3));
                uStimposdist   = squeeze(xmang{a}(d,tau,:,co));
                mcontprepost{p} = controlpos;
                uStrprepost{p}  = uStimposdist';
            end
        meancont(d,a,:)      = cat(2,mcontprepost{1},uStimp,mcontprepost{2});
        uStrials(d,a,co,:,:) = cat(2,uStrprepost{1},uStimptr,uStrprepost{2});
        end
    end
end

Tau = {1:tpre,tpost:ti};
angmeancont = nanmean(meancont,3); %(dims,ang,t)
for a=1:length(angle)
    for co=2:length(uStime)
        dist = zeros(length(angle),tr,size(normdimstime,2));
        for tri=1:tr
        for t=1:size(normdimstime,2)
            trial = squeeze(uStrials(:,a,co,tri,t));
            %compute distance wrt control vectors
            for ang=1:length(angle)
                dist(ang,tri,t) = norm(angmeancont(:,ang)-trial);
            end
        end
        end
        % assign class based on minimum distance to centroids
        [~,distmin] = min(dist,[],1);
        distmin = squeeze(distmin);
        idx = distmin == a;
        FAdec{i}(a,co-1,:,:) = idx;
    end
end
clear uStrials

end
    save([statspath 'subspaces/' statsfile],'FAst','FAcist','FAdec')

else
    disp('stats already computed')
end

clear FAst FAcist FAdec

end