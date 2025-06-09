% Dominant (FA) and memory (dPCA) subspaces explained variance

%% Paths and specs
run ../addpaths

[~, statspath, figurepath] = addpaths;

monkeys = {'Wa','Sa'};

numdims = 4; % mean optimal latent dimensionality

dimlab  = [num2str(numdims) 'Dmodel'];


%% Explained marginalized and total variance of held-out data

%  * explVar.totalVar             - total variance
%  * explVar.totalMarginalizedVar - total variance in each marginalization
%  * explVar.componentVar         - variance of each component (%)
%  * explVar.margVar              - variance of each component in each
%                                   marginalization (%), (marg x components)
%  * explVar.margVarsplit         - marginalized variance of each component in each
%                                   marginalization (%), (marg x components)

for m=1:length(monkeys)
    
    monkey = monkeys{m};
    
    disp(monkey)

    statsfile = [monkey '_dPCAstats_' dimlab '.mat'];
    load([statspath 'subspaces/' statsfile],'dPCAvar', 'FAvar','whichMargdPCA')

    F  = length(whichMargdPCA);
    dPCAv = []; dPCAmv = [];
    FAv   = []; FAmv = [];
    for f = 1:F
        dPCAv(:,:,f)  = dPCAvar{f}.margVar;
        dPCAmv(:,:,f) = dPCAvar{f}.margVarsplit;
        FAv(:,:,f)    = FAvar{f}.margVar;
        FAmv(:,:,f)   = FAvar{f}.margVarsplit;
    end

    % marginalized variance
    dPCAm  = squeeze(mean(dPCAmv,3));
    FAm    = squeeze(mean(FAmv,3));

    n = prod(size(dPCAmv,3));
    dPCAstd  = 1.96*squeeze(std(dPCAmv,0,3))/sqrt(n);
    n = prod(size(FAmv,3));
    FAstd    = 1.96*squeeze(std(FAmv,0,3))/sqrt(n);
    
    disp('% of marginalized variance captured by each subspace on each marginalization')
    
    disp('dominant subspace')
    fprintf('uStim var = %.0f ± %.0f \n',sum(FAm(1,:)),sum(FAstd(1,:)))
    fprintf('mem var = %.0f ± %.0f \n',sum(FAm(2,:)),sum(FAstd(2,:)))
    fprintf('time var = %.0f ± %.0f \n',sum(FAm(3,:)),sum(FAstd(3,:)))
    
    disp('memory subspace')
    fprintf('uStim var = %.0f ± %.0f \n',sum(dPCAm(1,:)),sum(dPCAstd(1,:)))
    fprintf('mem var = %.0f ± %.0f \n',sum(dPCAm(2,:)),sum(dPCAstd(2,:)))
    fprintf('time var = %.0f ± %.0f \n \n',sum(dPCAm(3,:)),sum(dPCAstd(3,:)))

    
    % total variance        
    dPCAm  = squeeze(mean(dPCAv,3));
    FAm    = squeeze(mean(FAv,3));

    n = prod(size(dPCAv,3));
    dPCAstd  = 1.96*squeeze(std(dPCAv,0,3))/sqrt(n);
    n = prod(size(FAv,3));
    FAstd    = 1.96*squeeze(std(FAv,0,3))/sqrt(n);
    
    disp('% of total variance captured by each subspace on each marginalization')

    disp('dominant subspace')
    fprintf('uStim var = %.0f ± %.0f \n',sum(FAm(1,:)),sum(FAstd(1,:)))
    fprintf('mem var = %.0f ± %.0f \n',sum(FAm(2,:)),sum(FAstd(2,:)))
    fprintf('time var = %.0f ± %.0f \n',sum(FAm(3,:)),sum(FAstd(3,:)))
    
    disp('memory subspace')
    fprintf('uStim var = %.0f ± %.0f \n',sum(dPCAm(1,:)),sum(dPCAstd(1,:)))
    fprintf('mem var = %.0f ± %.0f \n',sum(dPCAm(2,:)),sum(dPCAstd(2,:)))
    fprintf('time var = %.0f ± %.0f \n \n',sum(dPCAm(3,:)),sum(dPCAstd(3,:)))
end
