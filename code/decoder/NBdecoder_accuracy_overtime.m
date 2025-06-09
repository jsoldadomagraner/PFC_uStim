%% Decoder: plot classification accuracy over time for one example session

%% Paths and specs

run ../addpaths

[~, statspath, figurepath] = addpaths;
 
monkey   = 'Wa';
session  = 2;

savefig  = false;
savename = ['accuracy_overtime_' monkey];

% naive bayes model
NBmethod = 'Poisson';

% load decoder stats
[filenames,~,~,~] = datafiles(monkey);

filename = filenames{session};

fprintf('loading decoder stats for %s \n',filename)

statsfile = [filename '_' NBmethod(1) 'NBdecoder.mat'];

load([statspath 'decoder/' statsfile],'accuracy')


%% Compute stats

[T,S] = size(accuracy);

maxfolds   = zeros(T,S);
accuracyf  = cell(T,S);
accuracyf0 = cell(T,S);

for t=1:T
    for i=1:S
        maxfolds(t,i) = length(accuracy{t,i});
    end
end

maxfolds     = max(maxfolds,[],'all');
accsessions  = nan(T,S,maxfolds);

for t=1:T
    for i=1:S
        numfolds = length(accuracy{t,i});
        accsessions(t,i,1:numfolds) = accuracy{t,i};
    end
end

% Mean accuracy across folds for the control condition
accmean(:,1)  = nanmean(accsessions(:,1,:),[2,3]);
accmean0(:,1) = ones(T,1)*25;

% compute posterior estimates with confidence intervals
for t=1:T
    dat = accsessions(t,1,:)/100;
    dat = dat(~isnan(dat));
    [accpost(t,1), l(t,1), u(t,1), ~, ~] = bernoulliPost(dat(:));
end

% Mean accuracy across stimulations and folds
for s=2:S
accmean(:,s)  = nanmean(accsessions(:,s,:),[2,3]);
accmean0(:,s) = ones(T,1)*25;

% compute posterior estimates with confidence intervals
    for t=1:T
        dat = accsessions(t,s,:)/100;
        dat = dat(~isnan(dat));
        [accpost(t,s), l(t,s), u(t,s), ~, ~] = bernoulliPost(dat(:));
    end
end
accpost = accpost*100; l=(accmean-l*100); u=(u*100-accmean);

% Significance testing, binomial test
siglevel = 0.05;

% control
n_cont  = sum(~isnan(accsessions(:,1,:)),[2,3]);  % samples
k_cont  = nansum(accsessions(:,:,1,:)/100,[2,3]); % observed counts
p0_cont = ones(1,T)*0.25;                         % chance accuracy

for t=1:T
    pval_cont(t)=0;
    % test prob that p_cont > p0 (accuracy is bigger than chance)
    for i=k_cont(t):n_cont(t)
        pval_cont(t) = pval_cont(t) + binopdf(i,n_cont(t),p0_cont(t));
    end
end

% uStim
for s=2:S
    n_uS  = sum(~isnan(accsessions(:,s,:)),[2,3]); % samples
    k_uS  = nansum(accsessions(:,s,:)/100,[2,3]);  % observed counts
    p0_uS = ones(1,T)*0.25;                        % chance accuracy
    
    for t=1:T
        pval_uS(t,s-1)=0;
        % test prob that p_uS > p0 (accuracy is bigger than chance)
        for i=k_uS(t):n_uS(t)
            pval_uS(t,s-1) = pval_uS(t,s-1) + binopdf(i,n_uS(t),p0_uS(t));
        end
    end
end

%% Plots

smoothaccuracies = true;

[angcol,~,~,stimcol] = plottingspecs;
set(0,'DefaultAxesFontSize',22)

h=figure; hold on
xx=1.6; yy=1.2;
pos = get(h,'position');
set(h,'position',[pos(1:2) pos(3)*xx pos(4)*yy])


% trial period to analyse (all trial)
trialperiod = [-800,750]; %ms
[~,uStimp,trialtimes,tpre,~] = analysisperiod(trialperiod);
uStimp      = nan(1,length(uStimp));
tpost       = tpre+1;
targeton    = 5; % -600ms
targetoff   = 7; % -700ms

dimsh    = size(accmean);
dimsh(1) = length(uStimp);
uStimp1  = nan(dimsh);

dimsh    = size(pval_uS);
dimsh(1) = length(uStimp);
uStimp2  = nan(dimsh);

barw = 4;

% smoothing window size. 
% The standard deviation of the Gaussian for smoothdata is fixed to be 
% 1/5th of the total window width.
win  = 4; % 4 time points * 50ms each, 200 ms total, std window = 40 ms. 

% plot inferred posterior mean instead of empirical estimate of accuracy
accmeanuS  = cat(1,accpost(1:tpre,:),uStimp1,accpost(tpost:end,:));
accmeanuS0 = cat(1,accmean0(1:tpre,:),uStimp1,accmean0(tpost:end,:));
accsteuSu  = cat(1,u(1:tpre,:),uStimp1,u(tpost:end,:));
accsteuSl  = cat(1,l(1:tpre,:),uStimp1,l(tpost:end,:));

if smoothaccuracies
accmeanuS(1:tpre,:)  = smoothdata(accmeanuS(1:tpre,:),1,'gaussian',win);
accmeanuS0(1:tpre,:) = smoothdata(accmeanuS0(1:tpre,:),1,'gaussian',win);
accsteuSu(1:tpre,:)  = smoothdata(accsteuSu(1:tpre,:),1,'gaussian',win);
accsteuSl(1:tpre,:)  = smoothdata(accsteuSl(1:tpre,:),1,'gaussian',win);

accmeanuS(tpost+4:end,:)  = smoothdata(accmeanuS(tpost+4:end,:),1,'gaussian',win);
accmeanuS0(tpost+4:end,:) = smoothdata(accmeanuS0(tpost+4:end,:),1,'gaussian',win);
accsteuSu(tpost+4:end,:)  = smoothdata(accsteuSu(tpost+4:end,:),1,'gaussian',win);
accsteuSl(tpost+4:end,:)  = smoothdata(accsteuSl(tpost+4:end,:),1,'gaussian',win);
end

pval_cont = cat(2,pval_cont(1:tpre),uStimp,pval_cont(tpost:end));
pval_uS   = cat(1,pval_uS(1:tpre,:),uStimp2,pval_uS(tpost:end,:));

for s =1:S
shadedErrorBar(1:tpre,accmeanuS(1:tpre,s),[accsteuSu(1:tpre,s)'; accsteuSl(1:tpre,s)'],...
    'lineProps',{'LineStyle','-','Color',stimcol(s,:),'linewidth',1});
shadedErrorBar(tpost+3:T+4,accmeanuS(tpost+3:end,s),[accsteuSu(tpost+3:end,s)';accsteuSl(tpost+3:end,s)'], ...
    'lineProps',{'LineStyle','-','Color',stimcol(s,:),'linewidth',1});
p1(s) = plot(accmeanuS(:,s),'.-k','Linewidth',3,'MarkerSize',16);
p1(s).Color = stimcol(s,:);
if s==1;p2 = plot(accmeanuS0(:,s),'--k','Linewidth',3,'MarkerSize',16);end
end
ax = gca;
axis tight
yLim = ylim; xLim = xlim;
xlim([xLim(1)-0.5 xLim(2)+0.5])
ax.XTick = 5:4:T+4;
ax.XTickLabel = string(trialtimes(ax.XTick));
xlabel('time (ms)');
ylabel({'target angle';'classification accuracy (%)'})

yLim = ylim;

% plot significance bars
idx1     = 1:size(accmeanuS,1);
idx2     = 1:0.5:size(accmeanuS,1);
bar_cont = ones(1,length(idx2))*(yLim(2)+1);
sigidx   = pval_cont <= siglevel;
idx2nan  = find(sigidx==0);
[~,idx1nan,~] = intersect(idx2,idx2nan);
bar_cont(idx1nan)=nan;
bar_cont(tpre*2)=nan;
bar_cont(tpost*2+3*2)=nan;
plot(idx2,bar_cont,'Color',stimcol(1,:),'Linewidth',barw)

for s=1:S-1
idx1    = 1:size(accmeanuS,1);
idx2    = 1:0.5:size(accmeanuS,1);
bar_uS  = ones(1,length(idx2))*(yLim(2)+1+s);
sigidx  = pval_uS(:,s) <= siglevel;
idx2nan = find(sigidx==0);
[~,idx1nan,~] = intersect(idx2,idx2nan);
bar_uS(idx1nan)=nan;
bar_uS(tpre*2)=nan;
bar_uS(tpost*2+3*2)=nan;
plot(idx2,bar_uS,'Color',stimcol(s+1,:),'Linewidth',barw)
end

% plot trial events
ylim([10,yLim(2)+1+s])
yLim = ylim; xLim = xlim;
xlim([xLim(1)-0.5 xLim(2)+0.5])
plot([tpre(end)+1 tpre(end)+1],[yLim(1) yLim(2)], ...
    'Color','k', 'Linewidth',3)
plot([tpre(end)+4 tpre(end)+4],[yLim(1) yLim(2)],...
    'Color','k', 'Linewidth',3)
plot([targeton targeton],[yLim(1) yLim(2)], '--',...
    'Color','k', 'Linewidth',2)
plot([targetoff targetoff],[yLim(1) yLim(2)], '--',...
    'Color','k', 'Linewidth',2)

if savefig
    saveas(h,[figurepath 'decoder/' savename],'svg')
end
