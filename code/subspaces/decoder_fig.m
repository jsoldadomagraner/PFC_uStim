% Plot decoding accuracy in the dominant (FA) and memory (dPCA) subspaces

%% Paths and specs
run ../addpaths

[~, statspath, figurepath] = addpaths;

monkeys = {'Wa','Sa'};

numdims = 4; % mean optimal latent dimensionalit

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
load([statspath 'subspaces/' statsfile],'FAdec')

statsfile = [monkey '_dPCAstats_' dimlab '.mat'];    
load([statspath 'subspaces/' statsfile],'dPCAdec')


% FAdec{F}(ang,uStim,trials,time)
F    = length(FAdec);
dims = size(FAdec{1});
T    = dims(4);
Acc  = zeros([F 1 dims([2,4])]);
Samples = zeros(F,T);
Counts  = zeros(F,T);
for f=1:F
    Acc(f,:,:,:) = nanmean(FAdec{f},[1,3])*100;
    Samples(f,:) = sum(~isnan(FAdec{f}),[1,2,3]);
    Counts(f,:)  = nansum(FAdec{f},[1,2,3]);
end

avgdims = [1,2,3];
yLim    = [15 63];

h1 = figure; hold on
xx=1.4; yy=1.2;
pos = get(h1,'position');
set(h1,'position',[pos(1:2) pos(3)*xx pos(4)*yy])

% dominant subspace
dims    = size(Acc);
meanAcc = squeeze(nanmean(Acc,avgdims)); 
u       = 1.96*squeeze(nanstd(Acc,0,avgdims))./sqrt(prod(dims(avgdims)));
u       = u';
l       = u;
meanAcc(tpre+1:tpre+4) = nan;

disp(monkey)
fprintf('FA decoding accuracy: tpre=%.0f±%.0f tpost=%.0f±%.0f, tend=%.0f±%.0f \n', ...
    meanAcc(tpre), u(tpre), meanAcc(tpost+4), u(tpost+4), ... 
    meanAcc(end), u(end))

% significance level FA
n_FA  = sum(Samples,1); % samples
k_FA  = sum(Counts,1);  % observed counts
p0_FA = ones(1,T)*0.25;

for t=1:T
    pval_FA(t)=0;
    for i=k_FA(t):n_FA(t)
        pval_FA(t) = pval_FA(t) + binopdf(i,n_FA(t),p0_FA(t));
    end
end
pval_FA(tpre+1:tpre+4) = nan;


% Plotting
shadedErrorBar(1:tpre,meanAcc(1:tpre),[u(1:tpre); l(1:tpre)],...
    'lineProps',{'LineStyle','-','Color',col(1,:),'linewidth',1});
shadedErrorBar(tpost+3:T,meanAcc(tpost+3:end),[u(tpost+3:end);u(tpost+3:end)], ...
    'lineProps',{'LineStyle','-','Color',col(1,:),'linewidth',1});
a=plot(meanAcc,'Linewidth',3,'MarkerSize',ms, 'Marker',marker{m},...
         'MarkerFaceColor',col(1,:),'Color',col(1,:),'LineStyle',linestyle{m});


% memory subspace     
dims = size(dPCAdec{1});
T    = dims(4);
Acc  = zeros([F 1 dims([2,4])]);
Samples = zeros(F,T);
Counts  = zeros(F,T);
for f=1:F
    Acc(f,:,:,:) = nanmean(dPCAdec{f},[1,3])*100;
    Samples(f,:) = sum(~isnan(dPCAdec{f}),[1,2,3]);
    Counts(f,:)  = nansum(dPCAdec{f},[1,2,3]);
end

dims     = size(Acc);
meanAcc = squeeze(nanmean(Acc,avgdims));
u        = 1.96*squeeze(nanstd(Acc,0,avgdims))./sqrt(prod(dims(avgdims)));
u        = u';
l        = u;
meanAcc(tpre+1:tpre+4) = nan;

fprintf('dPCA decoding accuracy: tpre=%.0f±%.0f tpost=%.0f±%.0f, tend=%.0f±%.0f \n', ...
    meanAcc(tpre), u(tpre), meanAcc(tpost+4), u(tpost+4), ... 
    meanAcc(end), u(end))

% significance level dPCA
n_dPCA  = sum(Samples,1); % samples
k_dPCA  = sum(Counts,1);  % observed counts
p0_dPCA = ones(1,T)*0.25;

for t=1:T
    pval_dPCA(t)=0;
    for i=k_dPCA(t):n_dPCA(t)
        pval_dPCA(t) = pval_dPCA(t) + binopdf(i,n_dPCA(t),p0_dPCA(t));
    end
end
pval_dPCA(tpre+1:tpre+4) = nan;

   
shadedErrorBar(1:tpre,meanAcc(1:tpre),[u(1:tpre);l(1:tpre)],...
    'lineProps',{'LineStyle','-','Color',col(2,:),'linewidth',1});
shadedErrorBar(tpost+3:T,meanAcc(tpost+3:end),[u(tpost+3:end);l(tpost+3:end)], ...
    'lineProps',{'LineStyle','-','Color',col(2,:),'linewidth',1});
b=plot(meanAcc,'Linewidth',3,'MarkerSize',ms,'Marker',marker{m},...
    'MarkerFaceColor',col(2,:),'Color',col(2,:),'LineStyle',linestyle{m});

chance = ones(1,T)*25;
chance(tpre+1:tpre+4) = nan;
plot(chance,'--k','Linewidth',3)

axis tight

ylim(yLim)

% significance bar
siglevel = 0.05;
barw = 4;
for i=1:2
    if i==1
        sigidx  = pval_FA <= siglevel;
    else
        sigidx  = pval_dPCA <= siglevel;
    end
    idx1    = 1:length(meanAcc);
    idx2    = 1:0.5:length(meanAcc);
    bar_uS  = ones(1,length(idx2))*(yLim(2)+i*1);
    idx2nan = find(squeeze(sigidx)==0);
    [~,idx1nan,~] = intersect(idx2,idx2nan);
    bar_uS(idx1nan)=nan;
    bar_uS(tpre*2:tpost*2+3*2)=nan;
    bar_uS(tpre*2)=nan;
    bar_uS(tpost*2+3*2)=nan;
    plot(idx2,bar_uS,'Color',col(i,:),'Linewidth',barw)
end

ylim([yLim(1) yLim(2)+i*1])
ax = gca;
ax.XTick = 3:4:T;
ax.XTickLabel = trialtimes(ax.XTick);
xlabel('time (ms)')
ylabel({'target angle'; 'classification accuracy (%)'})

yLim = ylim; xLim = xlim;
xlim([xLim(1)-1 xLim(2)+1])
plot([taupre(end)+1 taupre(end)+1],[yLim(1) yLim(2)], ...
    'Color','k','LineWidth',2)
plot([taupre(end)+4 taupre(end)+4],[yLim(1) yLim(2)],...
    'Color','k','LineWidth',2)

if savefig
    saveas(h1,[figurepath 'subspaces/' monkey '_FA_dPCA_decoding'],'svg')
end

end
