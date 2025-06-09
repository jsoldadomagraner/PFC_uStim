%% FRs analysis: FRs stats

%% Paths and specs

run ../addpaths

[datapath, statspath, figurepath] = addpaths;

monkeys   = {'Wa','Sa'};

% analyse correct trials
datafolder = 'correct';

savefig    = false;


%% Plotting

h1=figure; hold on
xx=1.5; yy=1;
pos = get(h1,'position');
set(h1,'position',[pos(1:2) pos(3)*xx pos(4)*yy])

h2=figure; hold on
xx=1.5; yy=1;
pos = get(h2,'position');
set(h2,'position',[pos(1:2) pos(3)*xx pos(4)*yy])

a = 1; b = 2;
p = tiledlayout(a,b);
p.Padding = 'compact';
p.TileSpacing = 'compact';

h3=figure; hold on
xx=1; yy=1;
pos = get(h3,'position');
set(h3,'position',[pos(1:2) pos(3)*xx pos(4)*yy])

set(0,'DefaultAxesTitleFontWeight','normal');
set(0,'DefaultAxesFontSize',20)
fs  = 22; % font size for multipanel plots

marker    = {'sq','d'};
linestyle = {'-','-'};
legendlab = {'monkey W','monkey S'};
ms    = 8;
alpha = 0.2;
lw    = 3;
gry   = 0.6;
colm  = [[0 0 0];[gry gry gry]];

%% Compute FR stats

rebound_units   = cell(1,length(monkeys));
excited_units   = cell(1,length(monkeys));
inhibited_units = cell(1,length(monkeys));

for m = 1:length(monkeys)

monkey   = monkeys{m};

[filenames,numuStim,numang,stimorder] = datafiles(monkey);

% trial period to analyse
trialperiod = [-300,750]; %ms
[timep,uStimp,trialtimes,tpre,tpost] = analysisperiod(trialperiod);
timesms = trialtimes(timep);
T       = length(timep);

% uStim-induced modulation, Tuning strength (TS) and 
% uStim-induced changes in TS.
uStimmodtot         = nan(length(filenames),numuStim-1,length(timep));
uStimmod            = cell(1,length(filenames));
deltatuningstrength = cell(1,length(filenames));
tuningstrength0     = cell(1,length(filenames));

uStimmodang     = cell(1,length(filenames));
uStimmodangsign = cell(1,length(filenames));

% Assess prevalence of rebound effects after uStim, as well as
% uStim-induced excitation and inhibition
rebound_units{m}   = nan(length(filenames),numuStim-1);
excited_units{m}   = nan(length(filenames),numuStim-1);
inhibited_units{m} = nan(length(filenames),numuStim-1);

for f=1:length(filenames)
    
filename = filenames{f};

fprintf('loading file %s \n',filename)
load([datapath datafolder '/' filename],'spikerate', 'params')

if strcmp(monkey,'Sa')
    % in monkey Sa, padd short trials with nans
    load([datapath datafolder '/' filename],'behavior')
    binWidth  = params.binSize_spikerate*1000;
    spikerate = paddtrialswithnans(spikerate,2,behavior,binWidth);
end

% change order of uStim conditions to have control condition first
spikerate = spikerate(stimorder,:,:,:);

% remove shorted channels
chidx = selectchannels(monkey,datapath,filename);
numch = length(chidx);

trials = size(spikerate,4);
angle  = 1:numang; %target angles

% uStim modulation
% p-value < 0.05  (*). 
% p-value < 0.01  (**) 
% p-value < 0.001 (***)
pthresh = 0.05; % significance threshold

FRs  = cell(numch,numang);
FRs0 = cell(numch,numang);

uStimmod{f}            = nan(numch,numuStim-1,length(timep));
deltatuningstrength{f} = nan(numch,numuStim-1,length(timep));
tuningstrength0{f}     = nan(numch,1,length(timep));

uStimmodang{f}     = nan(numch,numuStim-1,numang,length(timep));
uStimmodangsign{f} = nan(numch,numuStim-1,numang,length(timep));

% Compute tuning strength
for ch = 1:numch
    
   ci  = chidx(ch);
    
   aidx = 1; %uStim amplitude 
   % Compute channel FRs for earch target angle
   for uidx = 1:numuStim
       FRch = nan(numang,trials,length(timep));
        for gidx = angle
            FRs{ch,gidx}  = nan(trials,length(timep));
            for r = 1:trials
                if ~isempty(spikerate{uidx,aidx,gidx,r})
                    FRs{ch,gidx}(r,:) = spikerate{uidx,aidx,gidx,r}(ci,timep);
                end
            end
            FRch(gidx,:,:)= FRs{ch,gidx};
        end     
        
        % Compute uStim effect across all angles
            if uidx==1
                uS0 = reshape(FRch,numang*trials,[]);
            else
                for t=1:T
                    uStim = reshape(FRch(:,:,t),numang*trials,[]);
                    cont  = uS0(:,t);
                    [pval,~] = ranksum(uStim,cont);
                    uStimmod{f}(ch,uidx-1,t) = pval;
                end
            end
        
            % Compute uStim effect per angle and direction of modulation
            % (excitation vs. inhibition) to estimate rebound effects
            if uidx==1
                uS0a = FRch;
            else
                for gidx = angle
                for t=1:T
                    uStim = squeeze(FRch(gidx,:,t));
                    cont  = squeeze(uS0a(gidx,:,t));
                    [pval,~] = ranksum(uStim,cont);
                    uStimmodang{f}(ch,uidx-1,gidx,t) = pval;
                    
                    uStimmean = nanmean(uStim);
                    contmean  = nanmean(cont);
                    uStimmodangsign{f}(ch,uidx-1,gidx,t) = ...
                        sign(uStimmean - contmean);
                end
                end
            end
        FRmean     = nanmean(FRch,2);
        
        % tuning strength
        TSalltr = squeeze(max(FRmean,[],1)-min(FRmean,[],1));

        if uidx==1
            TScont = TSalltr;
            tuningstrength0{f}(ch,1,:) = TScont;
        else            
            TSuStim = TSalltr;
            for t=1:T
                uSmod    = uStimmod{f}(ch,uidx-1,t)<pthresh;
                meandiff = TSuStim(t)-TScont(t);
                if uSmod==1
                    deltatuningstrength{f}(ch,uidx-1,t) = meandiff;
                end
            end
        end
   end
  
end
uStimmodtot(f,:,:) = mean(uStimmod{f} < pthresh,1);

% Estimate % of units with rebound effects 

% Here we use non-smoothed data. In Figure 2A, for visualization purposes,
% we have smoothed the firing rates and computed significance stats (bars)
% on the smoothed traces. This explains discrepancies with stats computed
% here (e.g., when considering the non-smoothed FRs, Unit 1 in Fig. 2A does
% not show significant rebound effects for target angle condition 45, 
% although it does for target angle 315).

% significant early excitation
idxe1 = uStimmodang{f}(:,:,:,tpost:tpost+1) < pthresh;
idxe2 = uStimmodangsign{f}(:,:,:,tpost:tpost+1) == +1;
idxee = idxe1 & idxe2;
win   = 2;                     % in at least two consecutive time bins
idxee = nansum(movsum(idxee,win,4)==win,4) >= 1; 

% significant early inhibition
idxi1 = uStimmodang{f}(:,:,:,tpost:tpost+1) < pthresh;
idxi2 = uStimmodangsign{f}(:,:,:,tpost:tpost+1) == -1;
idxi  = idxi1 & idxi2;
win   = 2;                     % in at least two consecutive time bins
idxi  = nansum(movsum(idxi,win,4)==win,4) >= 1; 

% significant later excitation
idxe1 = uStimmodang{f}(:,:,:,tpost+2:end) < pthresh;
idxe2 = uStimmodangsign{f}(:,:,:,tpost+2:end) == +1;
idxe  = idxe1 & idxe2;
win   = 2;                     % in at least two consecutive time bins
idxe  = nansum(movsum(idxe,win,4)==win,4) >= 1; 

% Rebound = union of early inhibition and later excitation 
idx   = idxi & idxe;

% consider presence of rebound effects if this occur in at least one target
% angle condition
idx = sum(idx,3) >= 1;
rebound_units{m}(f,:) = squeeze(mean(idx,1))*100;

% same for early inhibited units
idx = sum(idxi,3) >= 1;
inhibited_units{m}(f,:) = squeeze(mean(idx,1))*100;

% same for early excited units
idx = sum(idxee,3) >= 1;
excited_units{m}(f,:) = squeeze(mean(idx,1))*100;

end

%% Fig 1: % of units modulated by uStim over time

uStimmodmean = squeeze(nanmean(100*uStimmodtot,[1,2]));
uStimmodstd  = squeeze(nanstd(100*uStimmodtot,0,[1,2]));

uStimmodmean(uStimp) = nan;
uStimmodstd(uStimp)  = nan;

stats{1} = uStimmodmean;
stats{2} = uStimmodstd;

figure(h1)
shadedErrorBar(1:tpre,stats{1}(1:tpre),stats{2}(1:tpre), ...
    'lineProps',{'LineStyle','-','Color',colm(m,:),'linewidth',1}, ...
    'patchSaturation',alpha);
shadedErrorBar(tpost:T,stats{1}(tpost:end),stats{2}(tpost:end), ...
    'lineProps',{'LineStyle','-','Color',colm(m,:),'linewidth',1}, ...
    'patchSaturation',alpha);

l1(m)=plot(stats{1},'k','Linewidth',3,'Color',colm(m,:), 'MarkerSize',ms, ...
    'Marker',marker{m},'MarkerFaceColor',colm(m,:),'LineStyle',linestyle{m});

ax = gca;
ax.XTick = 3:4:T;
ax.XTickLabel = string(timesms(ax.XTick));
xlim([0,T+1])
xlabel('time (ms)')
ylabel('% of modulated units')

%% Fig 2: delta tuning strenght (uStim-no uStim) at tpost and tend

figure(h2)
tau = [tpost,T];
for i=1:2
    if m==1
        tl(i)=nexttile; hold on
    else
        nexttile(i)
    end
    dTS = [];
    for  f=1:length(filenames)
        dTSt = deltatuningstrength{f}(:,:,tau(i));
        dTS = [dTS;dTSt(:)]; %(ch,uidx-1,t)
    end
    if i==1
        l2(m)=histogram(dTS,'Normalization','probability', ...
            'DisplayStyle','stairs');
    else
        l2(m)=histogram(dTS,10,'Normalization','probability', ...
            'DisplayStyle','stairs');
    end
    l2(m).EdgeColor = colm(m,:);
    l2(m).LineWidth = lw;
    histmean(m,i) = nanmean(dTS);
    [~,pv(m,i)]   = ttest(dTS,0,'tail','left');
    ax = gca;
    ax.YTick = [0 0.05 0.1 0.15];
    if i==2; ax.YTick = [];end
end

%% Fig 3: tuning strenght in control conditions at tpost

figure(h3)
tau = tpost;
for i=1
    dTS = [];
    for  f=1:length(filenames)
        dTSt = tuningstrength0{f}(:,:,tau(i));
        dTS = [dTS;dTSt(:)]; %(ch,uidx-1,t)
    end
    l3(m)=histogram(dTS,7,'BinWidth',5,'Normalization','probability', ...
        'DisplayStyle','stairs');
    l3(m).EdgeColor = colm(m,:);
    l3(m).LineWidth = lw;
    histmean0(m,i)  = nanmean(dTS);
    histstd0(m,i)   = nanstd(dTS);
end


end
linkaxes(tl,'xy')

figure(h1)
yLim = ylim;
plot([tpre+1 tpre+1],[yLim(1) yLim(2)], ...
    'Color','k', 'Linewidth',1)
plot([tpost-1 tpost-1],[yLim(1) yLim(2)],...
    'Color','k', 'Linewidth',1)
legend(l1,legendlab,'box','off', ...
    'FontSize',14,'Location','NorthWest')


figure(h2)
for i=1:2
    nexttile(i);hold on   
    if i==1
    yLim = ylim;
    ylim([yLim(1) yLim(2)+0.01])
    end
    plot([0 0],[yLim(1) yLim(2)],'--','Color','k', 'Linewidth',1)
    for m=1:length(monkeys)
    plot(histmean(m,i),yLim(2)+0.01,'v','Markersize',14,'Color',colm(m,:), ...
        'MarkerFaceColor',colm(m,:))   
    end
end
xlabel(p,'\Delta tuning strength (uStim-no uStim) (sp/s)','FontSize',fs)
ylabel(p,'fraction of units','FontSize',fs)

disp('Is there a significant decrease in tuning strength on average across experiments for uStim-modulated units?')
fprintf('monkey W, tpost p = %.3f, tend p = %.3f \n',pv(1,1),pv(1,2))
fprintf('monkey S, tpost p = %.3f, tend p = %.3f \n',pv(2,1),pv(2,2))


figure(h3)
yLim = ylim;
ylim([yLim(1) yLim(2)+0.01])
for m=1:length(monkeys)
    plot(histmean0(m,1),yLim(2)+0.01,'v','Markersize',14,'Color',colm(m,:), ...
        'MarkerFaceColor',colm(m,:))
end
xlabel('tuning strength (sp/s)','FontSize',fs)
ylabel('fraction of units','FontSize',fs)

disp('mean tuning strength across units and sessions')
fprintf('monkey W = %.2f ± %.2f \n',histmean0(1,1),histstd0(1,1))
fprintf('monkey S = %.2f ± %.2f \n',histmean0(2,1),histstd0(2,1))


%% Rebound units stats

rbu     = [rebound_units{1}(:);rebound_units{2}(:)];
rebmean = mean(rbu);
rebstd  = std(rbu);

excu    = [excited_units{1}(:);excited_units{2}(:)];
excmean = mean(excu);
excstd  = std(excu);

inu    = [inhibited_units{1}(:);inhibited_units{2}(:)];
inmean = mean(inu);
instd  = std(inu);

fprintf('Percentage of excited units: %.0f ± %.0f (mean ± std) \n',  ...
    excmean,excstd)
fprintf('Percentage of inhibited units: %.0f ± %.0f (mean ± std) \n',  ...
    inmean,instd)
fprintf('Percentage of units with rebound effects: %.0f ± %.0f (mean ± std) \n',  ...
    rebmean,rebstd)

% Plotting
gap = 0.01;

% plot early excited units histogram
h4=figure; hold on
xx=1.2; yy=1.2;
pos = get(h4,'position');
set(h4,'position',[pos(1:2) pos(3)*xx pos(4)*yy])

l=histogram(excu,'Normalization','probability','DisplayStyle','stairs');
l.EdgeColor = 'k';
l.LineWidth = lw;

yLim = ylim;
ylim([yLim(1) yLim(2)+gap])
plot(excmean,yLim(2)+gap/2,'v','Markersize',14,'Color','k', ...
    'MarkerFaceColor','k')

title('uStim-induced early excitation')
xlabel('% of units with early excitatory effects','FontSize',fs)
ylabel('fraction of uStim experiments','FontSize',fs)


% plot early inhibited units histogram
h5=figure; hold on
xx=1.2; yy=1.2;
pos = get(h5,'position');
set(h5,'position',[pos(1:2) pos(3)*xx pos(4)*yy])

l=histogram(inu,'Normalization','probability','DisplayStyle','stairs');
l.EdgeColor = 'k';
l.LineWidth = lw;

yLim = ylim;
ylim([yLim(1) yLim(2)+gap])
plot(inmean,yLim(2)+gap/2,'v','Markersize',14,'Color','k', ...
    'MarkerFaceColor','k')

title('uStim-induced early inhibition')
xlabel('% of units with early inhibitory effects','FontSize',fs)
ylabel('fraction of uStim experiments','FontSize',fs)


% plot rebound units histogram
h6=figure; hold on
xx=1.2; yy=1.2;
pos = get(h6,'position');
set(h6,'position',[pos(1:2) pos(3)*xx pos(4)*yy])
    
l=histogram(rbu,'Normalization','probability','DisplayStyle','stairs');
l.EdgeColor = 'k';
l.LineWidth = lw;

yLim = ylim; xLim = xlim;
ylim([yLim(1) yLim(2)+gap])
plot(rebmean,yLim(2)+gap/2,'v','Markersize',14,'Color','k', ...
    'MarkerFaceColor','k')

title('uStim-induced rebound excitation')
xlabel('% of units with rebound effects','FontSize',fs)
ylabel('fraction of uStim experiments','FontSize',fs)

if savefig
    figurepath = [figurepath 'FRs/'];
    saveas(h1,[figurepath 'FRstats_uStimmod'],'svg')
    saveas(h2,[figurepath 'FRstats_tuningchanges'],'svg')
    saveas(h3,[figurepath 'FRstats_tuning'],'svg')
    saveas(h4,[figurepath 'FRstats_excitatory'],'svg')
    saveas(h5,[figurepath 'FRstats_inhibitory'],'svg')
    saveas(h6,[figurepath 'FRstats_rebound'],'svg')
end
