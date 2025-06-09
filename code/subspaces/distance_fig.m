% Plot distance between uStim and no uStim activity in the
% dominant (FA) and memory (dPCA) subspaces

%% Paths and specs
run ../addpaths

[~, statspath, figurepath] = addpaths;

monkeys = {'Wa','Sa'};

numdims = 4; % mean optimal latent dimensionality

dimlab  = [num2str(numdims) 'Dmodel'];

savefig = false;

% time period to plot
trialperiod = [-300,750]; %ms
[analyp,uStimp,trialtimes,tpre,~] = analysisperiod(trialperiod);
trialtimes = trialtimes(analyp);
taupre     = 1:tpre;
tpost      = tpre+1;
T          = length(trialtimes);

% plotting
set(0,'DefaultAxesFontSize',22)
colours   = lines(5);
marker    = {'sq','d'};
linestyle = {'-','-'};
legendlab = {'dominant','memory'};
col = colours([5 4],:);
ms  = 8;


for m = 1:length(monkeys)

monkey   = monkeys{m};

statsfile = [monkey '_FAstats_' dimlab '.mat'];    
load([statspath 'subspaces/' statsfile],'FAst','FAcist')

statsfile = [monkey '_dPCAstats_' dimlab '.mat'];    
load([statspath 'subspaces/' statsfile],'dPCAst','dPCAcist')

switch m
    case 1
        shift1 = 0.1;
        shift2 = 0.05;
    case 2
        shift1 = 0.11;
        shift2 = 0.05;
end

avgdims = [1,2,3];
yLabel = {'normalized distance';'between uStim and no uStim activity'};
if m==1; yLim   = [0.5 5.5]; yTicks = [0 1 2 3 4 5];end
if m==2; yLim   = [0.5 5.5]; yTicks = [0 1 2 3 4 5];end


%%  Distance between uStim and no uStim trajectories

h1 = figure; hold on
xx=1.6; yy=1.2;
pos = get(h1,'position');
set(h1,'position',[pos(1:2) pos(3)*xx pos(4)*yy])

% dominant subspace
dims = size(FAcist);
FAp  = squeeze(nanmean(FAcist,avgdims)); % variability across F and Ang
u    = 1.96*squeeze(nanstd(FAcist,0,avgdims))./sqrt(prod(dims(avgdims)));
u    = u';
l    = u;

shadedErrorBar(1:tpre,FAp(1:tpre),[u(1:tpre); l(1:tpre)],...
    'lineProps',{'LineStyle','-','Color',col(1,:),'linewidth',1});
shadedErrorBar(tpost+4:T,FAp(tpost+4:end),[u(tpost+4:end);u(tpost+4:end)], ...
    'lineProps',{'LineStyle','-','Color',col(1,:),'linewidth',1});
a=plot(FAp,'Linewidth',3,'MarkerSize',ms, 'Marker',marker{m},...
         'MarkerFaceColor',col(1,:),'Color',col(1,:),'LineStyle',linestyle{m});

% memory subspace
dims   = size(dPCAcist);
dPCAp  = squeeze(nanmean(dPCAcist,avgdims));
u      = 1.96*squeeze(nanstd(dPCAcist,0,avgdims))./sqrt(prod(dims(avgdims)));
u      = u';
l      = u;
   
shadedErrorBar(1:tpre,dPCAp(1:tpre),[u(1:tpre);l(1:tpre)],...
    'lineProps',{'LineStyle','-','Color',col(2,:),'linewidth',1});
shadedErrorBar(tpost+4:T,dPCAp(tpost+4:end),[u(tpost+4:end);l(tpost+4:end)], ...
    'lineProps',{'LineStyle','-','Color',col(2,:),'linewidth',1});
b=plot(dPCAp,'Linewidth',3,'MarkerSize',ms,'Marker',marker{m},...
    'MarkerFaceColor',col(2,:),'Color',col(2,:),'LineStyle',linestyle{m});
    
axis tight
ax = gca;

ylim(yLim)
ax.YTick = yTicks;
ax.XTick = 3:4:T;
ax.XTickLabel = string(trialtimes(ax.XTick));
xlabel('time (ms)')
ylabel({'normalized distance'; 'between uStim and'; 'no uStim activity'})

yLim = ylim; xLim = xlim;
xlim([xLim(1)-1 xLim(2)+1])
plot([taupre(end)+1 taupre(end)+1],[yLim(1) yLim(2)], ...
    'Color','k','LineWidth',2)
plot([taupre(end)+4 taupre(end)+4],[yLim(1) yLim(2)],...
    'Color','k','LineWidth',2)
sig05 = ones(1,T)*1;
sig05(tpre+1:tpre+4) = nan;
plot(sig05,'--','Color','k','LineWidth',2)

if savefig
    saveas(h1,[figurepath 'subspaces/' monkey '_FA_dPCA_distance'],'svg')
end


%% Percent of experiments with significant uStim modulation

h2 = figure; hold on
xx=1.6; yy=1.2;
pos = get(h2,'position');
set(h2,'position',[pos(1:2) pos(3)*xx pos(4)*yy])

% in the dominant and memory subspaces
FAp   = squeeze(nanmean(FAst,[1,2,3]))*100;
dPCAp = squeeze(nanmean(dPCAst,[1,2,3]))*100;

plot(1:T,FAp,'Linewidth',3,'MarkerSize',ms,'Marker',marker{m},...
    'MarkerFaceColor',col(1,:),'Color',col(1,:),'LineStyle',linestyle{m})
plot(1:T,dPCAp,'Linewidth',3,'MarkerSize',ms,'Marker',marker{m},...
    'MarkerFaceColor',col(2,:),'Color',col(2,:),'LineStyle',linestyle{m})

axis tight
ylim([0 100])

yLim = ylim; xLim = xlim;
xlim([xLim(1)-1 xLim(2)+1])
plot([taupre(end)+1 taupre(end)+1],[yLim(1) yLim(2)], ...
    'Color','k','LineWidth',2)
plot([taupre(end)+4 taupre(end)+4],[yLim(1) yLim(2)],...
    'Color','k','LineWidth',2)

ax = gca;
ax.XTick = 3:4:T;
ax.XTickLabel = string(trialtimes(ax.XTick));

xlabel('time (ms)')
ylabel({'% of uStim experiments'; 'with significant'; 'activity modulation'})

disp(monkey)
disp('% of experiments with significant uStim modulation')
fprintf('in the dominant subspace: tpost=%.0f, tend=%.0f \n', ...
    FAp(tpost+4), FAp(end))
fprintf('in the memory subspace: tpost=%.0f, tend=%.0f \n', ...
    dPCAp(tpost+4), dPCAp(end))

if savefig
    saveas(h2,[figurepath 'subspaces/' monkey '_FA_dPCA_percentsig'],'svg')
end

end
