% Compute alignment between dominant (FA) and memory (dPCA) subspaces

%% Paths and specs

run ../addpaths

[~, statspath, ~] = addpaths;

monkeys    = {'Wa','Sa'};

% use mean optimal latent dimensionality of the FA model
numdims = 4;

dimlab  = [num2str(numdims) 'Dmodel'];

FApath  = [statspath 'subspaces/FA_model_fits/'];

%% Compute subspace angle

for m = 1:length(monkeys)

monkey    = monkeys{m};

statsfile = [monkey '_dPCAstats_' dimlab '.mat'];    
load([statspath 'subspaces/' statsfile],'WdPCA','whichMargdPCA')

[filenames,~,~,~] = datafiles(monkey);

ln     = length(filenames);
subang = zeros(ln,1);

for f=1:ln

fprintf('%s \n', filenames{f})

FAmodel0 = [filenames{f} '_FAmodel_nouStim'];

load([FApath FAmodel0],'FAparams')

Marg     = whichMargdPCA{f};
W        = WdPCA{f};
FAsub    = FAparams(numdims).L;

% compute angle estimates for each dPCA cv-fold
for cv =1:size(W,3)
    
    mg = find(Marg(:,cv)==2);
    mg = mg(1:numdims);
    dPCAsub = squeeze(W(:,mg,cv));

    % returns the smallest angle across al subspace dimensions
    subang(f,cv) = subspace_min(dPCAsub,FAsub);

end
end

subang = subang*(180/pi);

fprintf('\n monkey %s \n', monkey)

% mean angle across sessions and cv-folds
subangmean = mean(subang,'all');
subangstd  = std(subang,0,'all');
fprintf('FA-dPCA min subspace angle = %.0f±%.0f \n \n', ...
  subangmean,subangstd)

end

