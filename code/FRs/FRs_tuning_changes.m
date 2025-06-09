%% FRs analysis: uStim-induced tuning changes, example units

%% Paths and specs

run ../addpaths

[datapath, statspath, figurepath] = addpaths;

% analyse correct trials
datafolder = 'correct';

monkey   = 'Wa';
session  = 8;

savefig  = false;
savename = ['FRs_tuning_changes_' monkey];

%% Load file 

[filenames,numuStim,numang,stimorder] = datafiles(monkey);

filename = filenames{session};

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
chidx     = selectchannels(monkey,datapath,filename);

% trial period to analyse
trialperiod = [-300,750]; %ms
[timep,uStimp,trialtimes,tpre,tpost] = analysisperiod(trialperiod);
timesms = trialtimes(timep);
T       = length(timep);

trials  = size(spikerate,4);
angle   = 1:numang; %target angles
    

%%  Plotting

% select channels to plot
chidxn = [30 56];

% select uStim conditions to plot for each channel.
Uidx   = [4 2];

h = figure;
xx=1.3; yy=1.5;
pos = get(h,'position');
set(h,'position',[pos(1:2) pos(3)*xx pos(4)*yy])

a = 2; b = 2;
p = tiledlayout(a,b);
p.Padding = 'compact';

[colang,~,~] = plottingspecs;

set(0,'DefaultAxesFontSize',14)
fs  = 15;

win = 4; %gaussian smoothing window
lw  = 2;
alpha = 0.2;

numch = length(chidxn);

% FRs in no uStim conditions
FRs0 = cell(numch,numang);

for ch = 1:numch
    
    ci  = chidxn(ch);
    
    tl(ch)=nexttile;
    hold on
    
    FR0 = [];
    aidx  = 1; %uStim amplitude
    for gidx = angle
        FRs0{ch,gidx} = nan(trials,length(timep));
        for r = 1:trials
            % control trials FRs
            if ~isempty(spikerate{1,aidx,gidx,r})
                FRs0{ch,gidx}(r,:) = spikerate{1,aidx,gidx,r}(ci,timep);
            end
        end
    meanFRs0 = nanmean(FRs0{ch,gidx},1);
    stdFRs0  = nanstd(FRs0{ch,gidx},0,1)/sqrt(trials);
    
    meanFRs0 = smoothdata(meanFRs0,'gaussian',win);
    stdFRs0  = smoothdata(stdFRs0,'gaussian',win);
    
    % compute tuning streght post smoothing only for visualization example
    FR0(gidx) = meanFRs0(tpost);

    shadedErrorBar(1:T,meanFRs0,stdFRs0, ...
    'lineProps',{'LineStyle','-','Color',colang(gidx,:),'linewidth',lw}, ...
    'patchSaturation',alpha);
    end
    ax = gca;
    ax.XTick = 3:4:T;
    ax.XTickLabel = [];
    axis tight

    % tuning strength
    TS = squeeze((max(FR0)-min(FR0)));

end

% uStim modulation
FRs  = cell(numch,numang);
FRs0 = cell(numch,numang);

for ch = 1:numch
    
    ci  = chidxn(ch);
    
    tl(ch+numch)=nexttile;
    hold on
    
    Uidxch = Uidx(ch);
    
    FR = [];
    for gidx = angle
        for uidx = Uidxch
            FRs{ch,gidx}  = nan(trials,length(timep));
            FRs0{ch,gidx} = nan(trials,length(timep));
            for r = 1:trials
                % uStim trials FRs
                if ~isempty(spikerate{uidx,aidx,gidx,r})
                    FRs{ch,gidx}(r,:) = spikerate{uidx,aidx,gidx,r}(ci,timep);
                    FRssmooth(r,:) = FRs{ch,gidx}(r,:);
                    FRssmooth(r,1:tpre)    = smoothdata(FRs{ch,gidx}(r,1:tpre),'gaussian',win);
                    FRssmooth(r,tpost:end) = smoothdata(FRs{ch,gidx}(r,tpost:end),'gaussian',win);
                end
                % no uStim trials FRs
                if ~isempty(spikerate{1,aidx,gidx,r})
                    FRs0{ch,gidx}(r,:) = spikerate{1,aidx,gidx,r}(ci,timep);
                    FRs0smooth(r,:) = FRs0{ch,gidx}(r,:);
                    FRs0smooth(r,1:tpre)    = smoothdata(FRs0{ch,gidx}(r,1:tpre),'gaussian',win);
                    FRs0smooth(r,tpost:end) = smoothdata(FRs0{ch,gidx}(r,tpost:end),'gaussian',win);
                end
            end
            meanFRs = nanmean(FRs{ch,gidx},1);
            stdFRs  = nanstd(FRs{ch,gidx},0,1)/sqrt(trials);
            
            meanFRs(uStimp) = nan;
            stdFRs(uStimp)  = nan;
            
            % uStim modulation, statistical significance
            % compute significance after smoothing for visualization
            for t=1:length(timep)
                [pval(ch,gidx,t),~] = ranksum(FRs0smooth(:,t),FRssmooth(:,t));
            end
            pval(ch,gidx,uStimp) = nan;
                                    
            meanFRs(1:tpre) = smoothdata(meanFRs(1:tpre),'gaussian',win);
            stdFRs(1:tpre)  = smoothdata(stdFRs(1:tpre),'gaussian',win);
            
            meanFRs(tpost:end) = smoothdata(meanFRs(tpost:end),'gaussian',win);
            stdFRs(tpost:end)  = smoothdata(stdFRs(tpost:end),'gaussian',win);
            
            % compute tuning streght post smoothing only for visualization example
            FR(gidx) = meanFRs(tpost);
            
            shadedErrorBar(1:tpre,meanFRs(1:tpre),stdFRs(1:tpre), ...
                'lineProps',{'LineStyle','-','Color',colang(gidx,:),'linewidth',lw}, ...
                'patchSaturation',alpha);
            shadedErrorBar(tpost:T,meanFRs(tpost:end),stdFRs(tpost:end), ...
                'lineProps',{'LineStyle','-','Color',colang(gidx,:),'linewidth',lw}, ...
                'patchSaturation',alpha);
            
        end
    end
    ax = gca;
    if ch==1 || ch==2
        ax.XTick = 3:4:T;
        ax.XTickLabel = string(timesms(ax.XTick));
    else
        ax.XTickLabel = [];
    end
    axis tight
    
    % tuning strength
    TS = squeeze((max(FR)-min(FR)));
end

for i=1:numch
    linkaxes(tl([i,numch+i]),'xy')
end

% plot significance bar
siglevel = 0.05;
barw     = 3;
for i=3:length(tl)
    nexttile(i);
    yLim = ylim;
    if i==3;offset=1.5;else;offset=1;end
    for gidx = angle
        idx1    = 1:size(meanFRs,2);
        idx2    = 1:0.5:size(meanFRs,2);
        bar     = ones(1,length(idx2))*(yLim(2)+gidx*offset);
        sigidx  = squeeze(pval(i-2,gidx,:)) <= siglevel;
        idx2nan = find(sigidx==0);
        [~,idx1nan,~] = intersect(idx2,idx2nan);
        bar(idx1nan)=nan;
        plot(idx2,bar,'Color',colang(gidx,:),'Linewidth',barw)
        ylim([yLim(1) yLim(2)+gidx*offset])
    end
end

for i=1:length(tl)
    nexttile(i);
    ax = gca;
    yLim = ylim;
    len = length(ax.YTick);
    ax.YTick = ax.YTick(1:ceil(len/4):len);
    plot([tpre+1 tpre+1],[yLim(1) yLim(2)], ...
        'Color','k', 'Linewidth',1)
    plot([tpost-1 tpost-1],[yLim(1) yLim(2)],...
        'Color','k', 'Linewidth',1)
end

xlabel(p,'time (ms)','FontSize',fs)
ylabel(p,'firing rates (sp/s)','FontSize',fs)

if savefig
    saveas(h,[figurepath 'FRs/' savename],'svg')
end

