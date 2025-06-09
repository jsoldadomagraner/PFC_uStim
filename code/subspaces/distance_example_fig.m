% Plot distance between uStim and no uStim activity in the
% dominant (FA) and memory (dPCA) subspaces
% Example session
% Distance computed considering the first two dimensions only

%% Paths and specs
run ../addpaths

[~, statspath, figurepath] = addpaths;

monkeys = {'Wa'};

session = 8;
angles  = [1,4]; % 45 and 315
uScond  = 3;

numdims = 4; % mean optimal latent dimensionality

dimlab  = [num2str(numdims) 'Dmodel_2Dstats'];

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
col  = colours([5 4],:);
ms   = 8;
lw   = 3;
barw = 4;

for m = 1:length(monkeys)

monkey   = monkeys{m};

statsfile = [monkey '_FAstats_' dimlab '.mat'];    
load([statspath 'subspaces/' statsfile],'FAst','FAcist')

statsfile = [monkey '_dPCAstats_' dimlab '.mat'];    
load([statspath 'subspaces/' statsfile],'dPCAst','dPCAcist')

YLim   = [0 6.1];
yTicks = [0 1 2 3 4 5 6];
s      = uScond;

for a = angles

    h1 = figure; hold on
    xx=1.6; yy=1.2;
    pos = get(h1,'position');
    set(h1,'position',[pos(1:2) pos(3)*xx pos(4)*yy])

for i=1:2    
    
    if i==1
        % Distance (F,ang,uStim,1,time)
        dist = squeeze(FAcist(session,a,s,:));   % 95% CI
        % Significant modulation
        sigm = squeeze(FAst(session,:,:,:,:));
    else
        % Distance (F,ang,uStim,1,time)
        dist = squeeze(dPCAcist(session,a,s,:)); % 95% CI
        sigm = squeeze(dPCAst(session,:,:,:,:));
    end
    
    plot(dist,'Linewidth',lw,'MarkerSize',ms, 'Marker',marker{m},...
      'MarkerFaceColor',col(i,:),'Color',col(i,:),...
      'LineStyle',linestyle{m});

    ylim(YLim)
    yLim = YLim;

    % significance bar
    idx1    = 1:size(dist,1);
    idx2    = 1:0.5:size(dist,1);
    bar_uS  = ones(1,length(idx2))*(yLim(2)+i*0.15);
    idx2nan = find(squeeze(sigm(a,s,:))==0);
    [~,idx1nan,~] = intersect(idx2,idx2nan);
    bar_uS(idx1nan)=nan;
    bar_uS(tpre*2:tpost*2+3*2)=nan;
    bar_uS(tpre*2)=nan;
    bar_uS(tpost*2+3*2)=nan;
    plot(idx2,bar_uS,'Color',col(i,:),'Linewidth',barw)

end    

axis tight
ax = gca;
yLim = ylim;
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
    saveas(h1,[figurepath 'subspaces/' monkey ...
        '_FA_dPCA_distance_example_ang' num2str(a)],'svg')
end

end

end
