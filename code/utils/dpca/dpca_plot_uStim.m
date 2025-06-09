function [h1,h2]=dpca_plot_uStim(Xfull, W, V, taupre, plotFunction, varargin)

% dpca_plot(X, W, V, plotFunction, ...) 
% produces a plot of the dPCA results. X is the data matrix, W and V
% are decoder and encoder matrices, plotFunction is a
% pointer to to the function that plots one component (see dpca_plot_default()
% for the template)

% dpca_plot(..., 'PARAM1',val1, 'PARAM2',val2, ...) 
% specifies optional parameter name/value pairs:
%
%  'whichMarg'              - which marginalization each component comes
%                             from. Is provided as an output of the dpca()
%                             function.
%
%  'time'                   - time axis
%
%  'timeEvents'             - time-points that should be marked on each subplot
%
%  'ylims'                  - array of y-axis spans for each
%                             marginalization or a single value to be used
%                             for each marginalization
%
%  'componentsSignif'       - time-periods of significant classification for each
%                             component. See dpca_signifComponents()
%
%  'timeMarginalization'    - if provided, it will be shown on top, and 
%                             irrespective of significance (because
%                             significant classification is not assessed for 
%                             time components)
%
%  'legendSubplot'          - number of the legend subplot
%
%  'marginalizationNames'   - names of each marginalization
%
%  'marginalizationColours' - colours for each marginalization
%
%  'explainedVar'           - structure returned by the dpca_explainedVariance
%
%  'numCompToShow'          - number of components to show on the explained
%                             variance plots (default = 15)
%
%  'X_extra'                - data array used for plotting that can be larger
%                             (i.e. have more conditions) than the one used
%                             for dpca computations
%  'showNonsignificantComponents'
%                           - display non-signficant components when there
%                             are fewer significant components than
%                             subplots

% Soldado-Magraner J. made modifications for custom plotting, 2025.

% default input parameters
options = struct('time',           [], ...   
                 'whichMarg',      [], ...
                 'timeEvents',     [], ...
                 'ylims',          [], ...
                 'componentsSignif', [], ...
                 'timeMarginalization', [], ...
                 'legendSubplot',  [], ...
                 'marginalizationNames', [], ...
                 'marginalizationColours', [], ...
                 'explainedVar',   [], ...
                 'numCompToShow',  20, ...
                 'X_extra',        [], ...
                 'showNonsignificantComponents', false);

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

% can't show more than there is
numCompToShow = min(options.numCompToShow, size(W,2));

X = Xfull(:,:)';
Xcen = bsxfun(@minus, X, mean(X));
XfullCen = bsxfun(@minus, Xfull, mean(X)');
N = size(X, 1);
dataDim = size(Xfull);
Z = Xcen * W;
%!!
%Z = bsxfun(@times, Z, 1./std(Z, [], 1));
%!!

toDisplayMargNames = 0;
    
componentsToPlot = [1 2 3];
Zfull = reshape(Z(:,componentsToPlot)', [length(componentsToPlot) dataDim(2:end)]);

if ~isempty(options.X_extra)
    XF = options.X_extra(:,:)';
    XFcen = bsxfun(@minus, XF, mean(X));
    ZF = XFcen * W;
    %!!
    %ZF = bsxfun(@times, ZF, 1./std(ZF, [], 1));
    %!!
    dataDimFull = size(options.X_extra);
    Zfull = reshape(ZF(:,componentsToPlot)', [length(componentsToPlot) dataDimFull(2:end)]);
end

tpre  = taupre(end);
tpost = tpre+1;

uStimp = [nan nan nan nan];
hpre  = Zfull(:,:,:,1:tpre);
hpost = Zfull(:,:,:,tpost:end);
dimsh = size(Zfull);
dimsh(end) = length(uStimp);
uStimph = nan(dimsh);
Zfull   = cat(4,hpre,uStimph,hpost);

set(0,'DefaultAxesFontSize',18)

h1 = figure; hold on
xx=2; yy=1;
pos = get(h1,'position');
set(h1,'position',[pos(1:2) pos(3)*xx pos(4)*yy])

a = 1; b = 3;
tl = gobjects(a,b);
p  = tiledlayout(a,b);
p.TileSpacing = 'compact';
p.Padding     = 'compact';

% y-axis spans
if isempty(options.ylims)
    options.ylims = max(abs(Zfull(:))) * 1.1;
end
if length(options.ylims) == 1
    if ~isempty(options.whichMarg)
        options.ylims = repmat(options.ylims, [1 max(options.whichMarg)]);
    end
end

% plotting all components as subplots
for c = 1:length(componentsToPlot)
    cc = componentsToPlot(c);
    nexttile
    
    if ~isempty(options.componentsSignif)
        signifTrace = options.componentsSignif(cc,:);
    else
        signifTrace = [];
    end
    
    if ~isempty(options.explainedVar)
        thisVar = options.explainedVar.componentVar(cc);
    else
        thisVar = [];
    end
    
    if ~isempty(options.whichMarg)
        thisYlim = options.ylims(options.whichMarg(cc));
        thisMarg = options.whichMarg(cc);
    else
        thisYlim = options.ylims;
        thisMarg = [];
    end
        
    dim = size(Xfull);
    cln = {c};
    for i=2:length(dim)
        cln{i} = ':';
    end

    % plot individual components using provided function
    plotFunction(Zfull(cln{:}), options.time, [-thisYlim thisYlim], ...
        thisVar, cc, options.timeEvents, ...
        signifTrace, thisMarg,c)
   
    % time period to plot
    trialperiod = [-300,750]; %ms
    [analyp,~,trialtimes,tpre,~] = analysisperiod(trialperiod);
    trialtimes = trialtimes(analyp);
    taupre     = 1:tpre;
    T          = length(trialtimes);
    
    ax = gca;
    ax.XTick = 3:4:T;
    ax.XTickLabel = string(trialtimes(ax.XTick));
    
    if c>1; ax.YTick = [];end
    
    yLim = ylim; xLim = xlim;
    xlim([xLim(1)-1 xLim(2)+1])
    plot([taupre(end)+1 taupre(end)+1],[yLim(1) yLim(2)], ...
        'Color','k','LineWidth',2)
    plot([taupre(end)+4 taupre(end)+4],[yLim(1) yLim(2)],...
        'Color','k','LineWidth',2)
    
    xlabel('time (ms)')

    if c==1; ylabel('Normalized firing rate (sp/s)');end

    if toDisplayMargNames && ~isempty(options.marginalizationNames)
        xx = xlim;
        yy = ylim;
        text(xx(1)+(xx(2)-xx(1))*0.1, yy(2)-(yy(2)-yy(1))*0.1, options.marginalizationNames(thisMarg))
    end
end 

% colours for marginalizations
if isempty(options.marginalizationColours)
    if ~isempty(options.explainedVar)
        L = length(options.explainedVar.totalMarginalizedVar);
        options.marginalizationColours = lines(L);
    elseif ~isempty(options.whichMarg)
        L = length(unique(options.whichMarg));
        options.marginalizationColours = lines(L);
    else
        options.marginalizationColours = [];
    end
end

h2 = figure; hold on
xx=1.2; yy=1;
pos = get(h2,'position');
set(h2,'position',[pos(1:2) pos(3)*xx pos(4)*yy])

% bar plot with projected variances
if ~isempty(options.explainedVar)
    hold on
    axis([0 numCompToShow+1 0 15])
    ylabel('Component variance (%)')
    xlabel('Component')
    b = bar(options.explainedVar.margVar(:,1:numCompToShow)' , 'stacked', 'BarWidth', 0.75);
    
    % fix for newer Matlab versions, May 2018, thanks to Tucker Fisher
    for idx = 1:numel(b)
        b(idx).FaceColor = options.marginalizationColours(idx,:);
    end    
    % end fix
    
    caxis([1 length(options.marginalizationColours)+256])
end

% pie chart
if ~isempty(options.explainedVar)
%     [left bottom width height]
    axes('position', [0.34 0.6 0.3 0.3])
    
    if isfield(options.explainedVar, 'totalMarginalizedVar_signal')
        d = options.explainedVar.totalMarginalizedVar_signal / options.explainedVar.totalVar_signal * 100;
       
        % In some rare cases the *signal* explained variances can be
        % negative (usually around 0 though); this means that the
        % corresponding marginalization does not carry [almost] any signal.
        % In order to avoid confusing pie charts, we set those to zero and
        % rescale the others to sum to 100%.
        if ~isempty(find(d<0, 1))
            d(d<0) = 0;
            d = d/sum(d)*100;
        end
    else
        d = options.explainedVar.totalMarginalizedVar / options.explainedVar.totalVar * 100;
    end
    
    % Rounding such that the rounded values still sum to 100%. Using
    % "largest remainder method" of allocation
    roundedD = floor(d);
    while sum(roundedD) < 100
        [~, ind] = max(d-roundedD);
        roundedD(ind) = roundedD(ind) + 1;
    end
    
    if ~isempty(options.marginalizationNames)
        for i=1:length(d)
            margNamesPerc{i} = [' '];
            margNamesPercd = [options.marginalizationNames{i} ' ' num2str(roundedD(i)) '%'];
            disp(margNamesPercd)
        end
    else
        for i=1:length(d)
            margNamesPerc{i} = [num2str(roundedD(i)) '%'];
        end
    end
    pieg = pie(d, ones(size(d)), margNamesPerc);
    idx = 0;
    for i = 1:length(pieg) 
        if mod(i,2)~=0
            idx = idx+1;
            pieg(i).FaceColor = options.marginalizationColours(idx,:);
        else
            pieg(i).FontSize = 12;
        end
    end
end
