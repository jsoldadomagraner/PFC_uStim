%% FRs analysis: FR tuning, example unit

%% Paths and specs

run ../addpaths

[datapath, statspath, figurepath] = addpaths;

% analyse correct trials
datafolder = 'correct';

monkey   = 'Wa';
session  = 8;

savefig  = false;
savename = ['FRs_tuning_' monkey];


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

% trial period to plot (use all trial in this case)
trialtimes = -800:50:800;
timep      = find(trialtimes>=-800 & trialtimes<=800);
T          = length(timep);
uStimp     = find(trialtimes>=0 & trialtimes<200);
[~,uStimp,~] = intersect(timep,uStimp);
% time steps before and after uStim
tpre  = uStimp(1)-1;
tpost = uStimp(end)+1;

trials  = size(spikerate,4);
angle   = 1:numang; %target angles


%% Plotting

% select channel to plot
chidxn = 30;
numch  = length(chidxn);

[colang,~,~] = plottingspecs;

set(0,'DefaultAxesFontSize',14)
fs  = 18;
win = 4; %gaussian smoothing window
lw  = 3;
alpha = 0.2;
yLim  = [0 40];

ch = 1;
ci = chidxn;

% figure size
xx=1; yy=1;

% PSTHs
FRs0  = cell(numch,numang);
FR    = [];
aidx  = 1;
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
FR(gidx) = meanFRs0(tpost);


h(gidx)=figure;
pos = get(h(gidx),'position');
set(h(gidx),'position',[pos(1:2) pos(3)*xx pos(4)*yy])
hold on

shadedErrorBar(1:T,meanFRs0,stdFRs0, ...
'lineProps',{'LineStyle','-','Color',colang(gidx,:),'linewidth',lw}, ...
'patchSaturation',alpha);

ax = gca;
axis tight
ax.XTick = 3:4:T;
ax.XTickLabel = trialtimes(ax.XTick);
ylim(yLim)

plot([tpre+1 tpre+1],[yLim(1) yLim(2)], ...
    'Color','k', 'Linewidth',1)
plot([tpost-1 tpost-1],[yLim(1) yLim(2)],...
    'Color','k', 'Linewidth',1)

plot([3 3],[yLim(1) yLim(2)], '--',...
    'Color','k', 'Linewidth',1)
plot([5 5],[yLim(1) yLim(2)],'--',...
    'Color','k', 'Linewidth',1)
plot([7 7],[yLim(1) yLim(2)],'--',...
    'Color','k', 'Linewidth',1)
plot([T-1 T-1],[yLim(1) yLim(2)],'--',...
    'Color','k', 'Linewidth',1)

xlabel('time (ms)','FontSize',fs)
ylabel('FRs (sp/s)','FontSize',fs)

set(gca,'Visible','off')

if savefig
    saveas(h(gidx),[figurepath 'FRs/' savename num2str(gidx)],'svg')
end

end

TS = squeeze((max(FR)-min(FR)));

fprintf('tuning strength at tpost = %.0f \n',TS)

