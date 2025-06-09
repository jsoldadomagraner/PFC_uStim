%% Decoder: plot classification accuracy across sessions, scatterplots

%% Paths and specs

run ../addpaths

[~, statspath, figurepath] = addpaths;

savefig   = false;
savename  = 'accuracy_scatterplots';

monkeys   = {'Wa','Sa'};

% naive bayes model
NBmethod  = 'Poisson';


for m=1:length(monkeys)
    
monkey   = monkeys{m};

[filenames,numuStim,numang,~] = datafiles(monkey);

% trial period to analyse (all trial)
trialperiod = [-800,750]; %ms
[timep,uStimp,~,tpre,~] = analysisperiod(trialperiod);
tpost = tpre+1;
T     = length(timep)-length(uStimp);

S = numuStim;
A = numang;
F = length(filenames);

maxfolds    = zeros(F,T,S);
accuracyf   = cell(F,T,S);
labels      = cell(F,T,S);
maxfoldsang = zeros(F,T,S,A);

% load decoder stats
for f=1:F
    
filename = filenames{f};

fprintf('loading decoder stats for %s \n',filename)

statsfile = [filename '_' NBmethod(1) 'NBdecoder.mat'];

load([statspath 'decoder/' statsfile],'accuracy','Y')

for i=1:S
    for t=1:T
        maxfolds(f,t,i)   = length(accuracy{t,i});
        accuracyf{f,t,i}  = accuracy{t,i};
        
        labels{f,t,i} = Y{i,t};
        for a=1:A
            idx = labels{f,t,i}==a;
            maxfoldsang(f,t,i,a) = sum(idx);
        end
        
    end
end

end

%% Compute decoding accuracy across sessions

maxfolds       = max(maxfolds,[],'all');
accsessions    = nan(T,F,S,maxfolds);

maxfoldsang    = max(maxfoldsang,[],'all');
accsessionsang = nan(T,F,S,A,maxfoldsang);

for f=1:F
    for t=1:T
        for i=1:S
            numfolds = length(accuracyf{f,t,i});
            accsessions(t,f,i,1:numfolds)  = accuracyf{f,t,i};
            for a=1:A
                idx = labels{f,t,i}==a;
                numfoldsa = sum(idx);
                accsessionsang(t,f,i,a,1:numfoldsa) = accuracyf{f,t,i}(idx);
            end
        end
    end
end

accmean    = zeros(F,S,T);
accmeanang = zeros(F,S,A,T);

for f=1:F
    for i=1:S
        % Mean accuracy across sessions and folds
        accmean(f,i,:)  = nanmean(squeeze(accsessions(:,f,i,:)),2);
    for a=1:A
        accmeanang(f,i,a,:)  = nanmean(squeeze(accsessionsang(:,f,i,a,:)),2);
    end
    end
end

Acc{m}    = accmean;
Accang{m} = accmeanang;
Exp{m}    = [F,S];

Accsess{m}    = accsessions;
Accsessang{m} = accsessionsang;

end


%% Plot

set(0,'DefaultAxesTitleFontWeight','normal');
set(0,'DefaultAxesFontSize',20)

% example session points to highlight
exsess = 2;

periods   = [tpre, tpost, T];
periodlab = {'t_{pre}','t_{post}','t_{end}'};

[angcol,markers,~,stimcol] = plottingspecs;

stimcol = stimcol(2:end,:);
gry2    = [0.7 0.7 0.7];

monkeysleg = {'monkey W','monkey S'};

fs = 22; % font size
sp = 2;  % make data breathe

h0=figure; hold on
xx = 2; yy=1;
pos = get(h0,'position');
set(h0,'position',[pos(1:2) pos(3)*xx pos(4)*yy])

pl = tiledlayout(1,3);
pl.Padding = 'compact';

for t = 1:length(periods)
    
tau = periods(t);    

fprintf('\n Period %i \n',t)

acccont{1} = repmat(squeeze(Acc{1}(:,1,tau)),[1 Exp{1}(2)-1]); %(F,S,T);
acccont{2} = repmat(squeeze(Acc{2}(:,1,tau)),[1 Exp{2}(2)-1]);

accuS{1} = squeeze(Acc{1}(:,2:end,tau));
accuS{2} = squeeze(Acc{2}(:,2:end,tau));

X = [acccont{1}(:) ; acccont{2}(:)];
Y = [accuS{1}(:) ; accuS{2}(:)];

% Two-tailed paired Wilcoxon rank sum test
% Test of the null hypothesis that the two samples have the same 
% distribution, against the alternative that they do not.
[p,hy]   = ranksum(X,Y);
fprintf('rank sum test p=%.2f \n',p)

% Linear regression
disp('linear regression')
R = [X, ones(size(X))];
[b,~,~,~,stats] = regress(Y,R);
fprintf('slope=%.2f, offset=%.2f, p=%.2e \n', b(1), b(2), stats(3))
Yhat = R*b;

diffm = mean(Y-X);
fprintf('accuracy difference = %.2f  (uStim - no uStim) \n',diffm)

% Binomial test for each uStim pattern
siglev = 0.05;

Pvall = [];
Pvalh = [];
for m=1:length(monkeys)
    
    F = size(Accsess{m},2);
    S = size(Accsess{m},3);
    
    pvall = zeros(F,S-1);
    pvalh = zeros(F,S-1);
    
    for f=1:F
        
    n = sum(~isnan(Accsess{m}(tau,f,1,:)));
    k = nansum(Accsess{m}(tau,f,1,:)/100);

    n_cont  = n; % samples
    k_cont  = k; % observed counts
    p0      = k_cont/n_cont; % no uStim accuracy
        
    for s=1:S-1
        n = sum(~isnan(Accsess{m}(tau,f,s+1,:)));
        k = nansum(Accsess{m}(tau,f,s+1,:)/100);

        n_uS  = n; % samples
        k_uS  = k; % observed counts

        pvall(f,s)=0;
        % test prob that p_uS < p0 (accuracy is lower in uStim)
        for i=0:k_uS
            pvall(f,s) = pvall(f,s) + binopdf(i,n_uS,p0);
        end
        % test prob that p_uS > p0 (accuracy is bigger in uStim)
        pvalh(f,s)=0;
        for i=k_uS:n_uS
            pvalh(f,s) = pvalh(f,s) + binopdf(i,n_uS,p0);
        end
    end
    
    end
    Pvall = [Pvall pvall(:)'];
    Pvalh = [Pvalh pvalh(:)'];
    Pvallm{m} = pvall<siglev;
    Pvalhm{m} = pvalh<siglev;
end

sigpvall = 100*sum(Pvall<siglev)/length(Pvall);
sigpvalh = 100*sum(Pvalh<siglev)/length(Pvalh);

disp('% of uStim patterns with accuracy in uStim significantly lower than in no-uStim')
fprintf('Percent of patterns = %.2f \n',sigpvall)

disp('% of uStim patterns with accuracy in uStim significantly higher than in no-uStim')
fprintf('Percent of patterns = %.2f \n',sigpvalh)


% Scatter plots
% figure(h0)
tl(t) = nexttile; hold on

for m=1:length(monkeys)
        
acccontrol = acccont{m};
accuStim   = accuS{m};

plot(acccontrol(:),accuStim(:),'o','MarkerFaceColor',gry2, ...
    'MarkerEdgeColor',gry2,'Marker',markers{m})

% color cases with significant increases/decreases in accuracy wrt control
plot(acccontrol(Pvallm{m}),accuStim(Pvallm{m}),'o','MarkerFaceColor','k', ...
    'MarkerEdgeColor','k','Marker',markers{m})
plot(acccontrol(Pvalhm{m}),accuStim(Pvalhm{m}),'o','MarkerFaceColor','k', ...
    'MarkerEdgeColor','k','Marker',markers{m})

if m==1
    plot(acccontrol(exsess,1:3),accuStim(exsess,1:3),'o',...
        'MarkerEdgeColor',gry2,'MarkerSize',16)
    
    idx = Pvallm{m}(exsess,1:3);
    a   = acccontrol(exsess,1:3);
    b   = accuStim(exsess,1:3);
    plot(a(idx),b(idx),'o','MarkerEdgeColor','k','MarkerSize',16)
    idx = Pvalhm{m}(exsess,1:3);
    plot(a(idx),b(idx),'o','MarkerEdgeColor','k','MarkerSize',16)
    
    for s=1:3
    plot(acccontrol(exsess,s),accuStim(exsess,s),'o', ...
        'MarkerFaceColor',stimcol(s,:), ...
        'MarkerEdgeColor',stimcol(s,:),'Marker',markers{m})
    end
end

end
axis tight

title(periodlab{t});

% make axis equal
xLim = xlim;
yLim = ylim;
minLim = min([xLim yLim]);
maxLim = max([xLim yLim]);
xlim([minLim maxLim])
ylim([minLim maxLim])

end

linkaxes(tl,'xy')
xLim = xlim;
yLim = ylim;
minLim = min([xLim yLim]);
maxLim = max([xLim yLim]);

for i=1:length(tl)
    nexttile(i)
    % plot unity line
    plot([minLim maxLim],[minLim maxLim],'k-','LineWidth',2)
    % plot chance levels
    plot([xLim(1) xLim(2)],[25 25],'--k','LineWidth',1)
    plot([25 25],[yLim(1) yLim(2)],'--k','LineWidth',1)
    % make data breathe
    xLim = xlim;
    yLim = ylim;
    xlim([xLim(1)-sp, xLim(2)+sp])
    ylim([yLim(1)-sp, yLim(2)+sp])

if i==1
    % dummy plot for legend
    l1 = plot(nan,nan,'o','MarkerFaceColor','k', ...
        'MarkerEdgeColor','k','Marker',markers{1});
    l2 = plot(nan,nan,'o','MarkerFaceColor','k', ...
        'MarkerEdgeColor','k','Marker',markers{2});
    legend([l1,l2],monkeysleg,'box','off','Location','NorthWest', ...
        'FontSize',14,'Color','k')
end

end

for i=1:length(tl)
    nexttile(i)
    ax = gca;
    ax.XTick = ax.XTick([1,3,5]);
    ax.YTick = ax.YTick([1,3,5]);
    if i~=1
        ax.YTick = [];
    end
end

xlabel(pl,'classification accuracy, no uStim (%)','FontSize',fs)
ylabel(pl,{'classification accuracy,'; 'uStim (%)'},'FontSize',fs)

if savefig
    saveas(h0,[figurepath 'decoder/' 'scatter'],'svg')
end

%% Effect of uStim electrode tuning on decoding accuracy

% In particular, we plot the change in accuracy (uStim - no uStim) as a 
% function of the angular distance between the vector representing the 
% spatial tuning of the stimulated site and the vector representing the 
% target angle presented to the monkey. 
% Thus, we ask whether the difference between the target angle and the 
% tuning present at the stimulated location was impactful for decoding.
load([statspath 'FRs/' 'uStime_tuning.mat'],'targetdist')

h1=figure; hold on
xx = 2; yy=1.1;
pos = get(h1,'position');
set(h1,'position',[pos(1:2) pos(3)*xx pos(4)*yy])

clear tl pl

pl = tiledlayout(1,3);
pl.Padding = 'compact';

periods   = [tpre, tpost, T];
periodlab = {'t_{pre}','t_{post}','t_{end}'};

disp(' ')
disp('uStime tuning effect on accuracy')

for t = 1:length(periods)
    
tau = periods(t);    

fprintf('\n Period %i \n',t)

acccont{1} = repmat(Accang{1}(:,1,:,tau),[1 Exp{1}(2)-1 1 1]); %(F,S,A,T);
acccont{2} = repmat(Accang{2}(:,1,:,tau),[1 Exp{2}(2)-1 1 1]);

accuS{1} = Accang{1}(:,2:end,:,tau);
accuS{2} = Accang{2}(:,2:end,:,tau);

X = [targetdist{1}(:) ; targetdist{2}(:)];
Y = [accuS{1}(:) - acccont{1}(:) ; accuS{2}(:) - acccont{2}(:)];

idx = isnan(X);
X = X(~idx);
Y = Y(~idx);

[rho,pval]   = corr([X Y]);

fprintf('Pearson correlation r=%.2f \n',rho(1,2))
fprintf(['p = ' num2str(pval(1,2),'%.2f') ' \n'])

% line fit
R = [X, ones(size(X))];
B = R \ Y;
Yhat = R * B;

% Scatter plots
tl(t) = nexttile; hold on

for m=1:length(monkeys)
    for a = 1:A        
    acccontrol = acccont{m}(:,:,a);
    accuStim   = accuS{m}(:,:,a);

    deltaacc   = accuStim - acccontrol;
    targdis    = targetdist{m}(:,:,a);

    plot(targdis(:),deltaacc(:),'o','MarkerFaceColor',angcol(a,:), ...
        'MarkerEdgeColor',angcol(a,:),'Marker',markers{m})
    end
end
axis tight

% line fit
plot(X,Yhat,'k-','Linewidth',2)

title(periodlab{t});

end

linkaxes(tl,'xy')
xLim = xlim;
yLim = ylim;
minLim = min([xLim yLim]);
maxLim = max([xLim yLim]);

for i=1:length(tl)
    nexttile(i)
    % make data breathe
    xLim = xlim;
    yLim = ylim;
    gapx = 5; gapy = 1; 
    xlim([xLim(1)-gapx, xLim(2)+gapx])
    ylim([yLim(1)-gapy, yLim(2)+gapy])
end

for i=1:length(tl)
    nexttile(i)
    ax = gca;
    if i~=1
        ax.YTick = [];
    end
end

title(pl,'uStim electrode tuning direction & classification accuracy', 'FontSize', fs);
xlabel(pl,'Angle to uStim electrode tuning vector (\circ)','FontSize',fs)
ylabel(pl,{'\Delta classification accuracy (%)'; '(uStim - no uStim)'},'FontSize',fs)

if savefig
    saveas(h1,[figurepath 'decoder/' 'scatter_tuning'],'svg')
end


%% Histograms with (no uStim - uStim) distribution difference

for t = 1:length(periods)
    
tau = periods(t);    

acccont{1} = repmat(squeeze(Acc{1}(:,1,tau)),[1 Exp{1}(2)-1]); %(F,S,T);
acccont{2} = repmat(squeeze(Acc{2}(:,1,tau)),[1 Exp{2}(2)-1]);

accuS{1} = squeeze(Acc{1}(:,2:end,tau));
accuS{2} = squeeze(Acc{2}(:,2:end,tau));

X = [acccont{1}(:) ; acccont{2}(:)];
Y = [accuS{1}(:) ; accuS{2}(:)];

h=figure; hold on
xx = 1; yy=1;
pos = get(h,'position');
set(h,'position',[pos(1:2) pos(3)*xx pos(4)*yy])

hi = histogram(X-Y,'DisplayStyle','stairs');
hi.EdgeColor = 'k';
hi.LineWidth = 2;
yLim = ylim;
ylim([yLim(1) yLim(2)+4])
yLim = ylim;
diffm = mean(X-Y);
plot([0 0],[yLim(1) yLim(2)-2],'k','Linewidth',5)
plot(diffm,yLim(2),'v','Markersize',24,'Color','k', ...
    'MarkerFaceColor','k')
axis off

if savefig
    saveas(h,[figurepath 'decoder/' ...
        'scatter_hist_t' num2str(t)],'svg')
end
end