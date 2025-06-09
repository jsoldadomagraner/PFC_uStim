function explVar = dpca_explainedVariance_CVtest(Xfull, Xfull_hat, varargin)

% explVar = dpca_explainedVariance(X, W, V) computes various measures of
% explained variance and returns them in a structure explVar. X is the data
% matrix, W is the decoder matrix, V is the encoder matrix. Returned values:
%
%  * explVar.totalVar             - total variance
%  * explVar.totalMarginalizedVar - total variance in each marginalization
%  * explVar.componentVar         - variance of each component (%)
%  * explVar.margVar              - variance of each component in each marginalization (%)
%  * explVar.margVarsplit         - marginalized variance of each component in each marginalization (%)
%  * explVar.cumulativePCA        - cumulative variance of the PCA components (%)
%  * explVar.cumulativeDPCA       - cumulative variance of the dPCA components (%)
%
% [...] = dpca(..., 'PARAM1',val1, 'PARAM2',val2, ...) 
% specifies optional parameter name/value pairs:
%
%  'combinedParams' - cell array of cell arrays specifying 
%                     which marginalizations should be added up together,
%                     e.g. for the three-parameter case with parameters
%                           1: stimulus
%                           2: decision
%                           3: time
%                     one could use the following value:
%                     {{1, [1 3]}, {2, [2 3]}, {3}, {[1 2], [1 2 3]}}.


% default input parameters
options = struct('combinedParams', [], ...   
                 'X_trial',        [], ...
                 'numOfTrials',    [], ...
                 'Cnoise',         []);

% read input parameters
optionNames = fieldnames(options);
if mod(length(varargin),2) == 1
	error('Please provide propertyName/propertyValue pairs')
end
for pair = reshape(varargin,2,[])    % pair is {propName; propValue}
	if any(strcmp(pair{1}, optionNames))
        options.(pair{1}) = pair{2};
    else
        error('%s is not a recognized parameter name', pair{1})
	end
end

% centering
X = Xfull(:,:);
Xfull = bsxfun(@minus, Xfull, mean(X,2));
X = bsxfun(@minus, X, mean(X,2));

% marginalizing
Xmargs = dpca_marginalize(Xfull, 'combinedParams', options.combinedParams, 'ifFlat', 'yes');

% total variance
explVar.totalVar = sum(sum(X.^2));

% total marginalized variance
for i=1:length(Xmargs)
    explVar.totalMarginalizedVar(i) = sum(Xmargs{i}(:).^2);
end

% dims
dims = size(Xfull_hat,5);

% PCA explained variance
[~,S,~] = svd(X', 'econ');
S = diag(S);
S = S(1:dims);
explVar.cumulativePCA = cumsum(S.^2'/ explVar.totalVar * 100);

% dPCA explained variance
for i=1:dims
    
    Xhati   = squeeze(Xfull_hat(:,:,:,:,i));
    Xhatimx = Xhati(:,:);
    Xhati   = bsxfun(@minus, Xhati, mean(Xhatimx,2));
    Xhatimx = bsxfun(@minus, Xhatimx, mean(Xhatimx,2));
    
    Xhat(:,:,i) = Xhatimx;

    explVar.cumulativeDPCA(i) = 100 - sum(sum((X - squeeze(sum(Xhat(:,:,1:i),3))).^2)) / explVar.totalVar * 100;    
    explVar.componentVar(i) = 100 - sum(sum((X - Xhatimx).^2)) / explVar.totalVar * 100;    
   
    % marginalizing
    Xmargshati = dpca_marginalize(Xhati, 'combinedParams', options.combinedParams, 'ifFlat', 'yes');
    for j=1:length(Xmargs)
        ZZ = Xmargs{j} - Xmargshati{j};
        explVar.margVar(j,i) = (explVar.totalMarginalizedVar(j) - sum(ZZ(:).^2)) / explVar.totalVar * 100;    
        explVar.margVarsplit(j,i) = (explVar.totalMarginalizedVar(j) - sum(ZZ(:).^2)) / explVar.totalMarginalizedVar(j) * 100;
    end
end

% if dPCA expl var = PCA expl var, then probably the function was supplied
% with PCA axes instead of dPCA encoder/decoder. In this case it will not
% return the dPCA explained variance.
if max(abs(explVar.cumulativePCA-explVar.cumulativeDPCA)) < 1e-10
    explVar.cumulativeDPCA = [];
end


end
