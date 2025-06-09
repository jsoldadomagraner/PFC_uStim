%% Behavior analysis: errors

%% Paths and specs

run ../addpaths

[datapath, statspath, figurepath] = addpaths;

monkeys    = {'Wa','Sa'};

savefig    = false;
savename   = 'errors';

% Take errors from complete trials only (exclude fix breaks and false
% starts, which are essencially fix breaks after go cue).
% Target breaks count as completed trials, but we won't consider them as 
% errors, given that monkeys correctly remembered the target location.
% Thus as errors, we will consider only "target misses". 
% 154, trialtype = 'error trials (broke target)'
% 157, trialtype = 'error trials (no choice/miss target/wrong target)'
% 160, trialtype = 'error trials (false start)'
errcode1 = 154; % Include target breaks as correctly completed trials
errcode2 = 157; % Target misses due to 1. not performing a choice after go 
                % cue, 2. saccading but missing the target or 3. selecting
                % the wrong target

% Plotting
[col_ang, markers, dispName] = plottingspecs;


%% Compute stats

statsfile = [statspath 'behavior/' 'errstats.mat'];
if ~exist(statsfile,'file')

disp('Computing error stats...')

Errorcount     = cell(1,length(monkeys));
Correctcount   = cell(1,length(monkeys));
Breaktargcount = cell(1,length(monkeys));

for m=1:length(monkeys)

monkey = monkeys{m};

[filenames,numuStim,numang,stimorder] = datafiles(monkey);

fn = length(filenames);

errorcount     = zeros(fn,numuStim,numang);
breaktargcount = zeros(fn,numuStim,numang);
correctcount   = zeros(fn,numuStim,numang);

for i=1:fn
    
filename   = [filenames{i} '_err'];

fprintf('loading file %s \n',filename)
datafolder = 'error';
load([datapath datafolder '/' filename],'behavior')

% change order of uStim conditions to have control condition first
behavior(:,:,:,:)  = behavior(stimorder,:,:,:);

trials = size(behavior,4);

%% Count error trials and target breaks

am = 1; % uStim amplitude is fixed
for u = 1:numuStim
    for an = 1:numang
        for tr = 1:trials
            if ~isempty(behavior{u,am,an,tr})
                codes = behavior{u,am,an,tr}.codes;
                if sum(codes(:,1)==errcode1,'all')~=0
                    breaktargcount(i,u,an) = breaktargcount(i,u,an) + 1;
                end
                if sum(codes(:,1)==errcode2,'all')~=0
                    errorcount(i,u,an) = errorcount(i,u,an) + 1;
                end
            end
        end
    end
end

%% Count correct trials

filename   = filenames{i};
datafolder = 'correct';
load([datapath datafolder '/' filename],'spikerate')

% change order of uStim conditions to have control condition first
spikerate(:,:,:,:)  = spikerate(stimorder,:,:,:);

trials   = size(spikerate,4);

am = 1; 
for u = 1:numuStim
    for an = 1:numang
        for tr = 1:trials
            if ~isempty(spikerate{u,am,an,tr})
                correctcount(i,u,an) = correctcount(i,u,an) + 1;
            end
        end
    end
end

end

Errorcount{m}     = errorcount;
Correctcount{m}   = correctcount;
Breaktargcount{m} = breaktargcount;

end
    save(statsfile,'Errorcount','Correctcount','Breaktargcount')
else
    load(statsfile)
end

%% Compute error fraction

for m=1:length(monkeys)
    
[~,numuStim, numang] = size(Errorcount{m});

dims = 3;

% Compute errors pulling across all target angle conditions
% in uStim and no-uStim conditions
Err{m,2} = sum(Errorcount{m}(:,2:end,:),dims)./(sum(Errorcount{m}(:,2:end,:),dims) + ...
    sum(Breaktargcount{m}(:,2:end,:),dims) + sum(Correctcount{m}(:,2:end,:),dims));
Err{m,1} = sum(Errorcount{m}(:,1,:),dims)./(sum(Errorcount{m}(:,1,:),dims) + ...
    sum(Breaktargcount{m}(:,1,:),dims) + sum(Correctcount{m}(:,1,:),dims));
Err{m,1} = repmat(Err{m,1},[1 numuStim-1 1]);

% Compute errors per target angle conditions
% in uStim and no-uStim conditions
Err_ang{m,2} = Errorcount{m}(:,2:end,:)./(Errorcount{m}(:,2:end,:) + ...
    Breaktargcount{m}(:,2:end,:) + Correctcount{m}(:,2:end,:));
Err_ang{m,1} = Errorcount{m}(:,1,:)./(Errorcount{m}(:,1,:) + ...
    Breaktargcount{m}(:,1,:) + Correctcount{m}(:,1,:));
Err_ang{m,1} = repmat(Err_ang{m,1},[1 numuStim-1 1]);

end

%% Test whether uStim and no uStim error distributions are statistically different

X = [Err{1,1}(:) ; Err{2,1}(:)]; % no uStim
Y = [Err{1,2}(:) ; Err{2,2}(:)]; % uStim

% Two-tailed paired Wilcoxon rank sum test
% Test of the null hypothesis that the two samples have the same 
% distribution, against the alternative that they do not.
[p,hy]   = ranksum(X,Y);
fprintf('rank sum test p=%.2f \n \n',p)

% Linear regression
disp('linear regression')
% Compute p-values of the full model based on the F-statistic (H0: reduced 
% model, just intercept; HA: full model, intercept + slope).
R = [X, ones(size(X))];
[b,~,~,~,stats] = regress(Y,R);
fprintf('slope=%.2f, offset=%.2f, p=%.2e \n', b(1), b(2), stats(3))


%% Plot errors (target misses) in uStim vs. no uStim

h=figure; hold on
pos = get(h,'position');
set(h,'position',[pos(1:2) pos(3)*1 pos(4)*1])

% font size
fs = 20;

for m=1:length(monkeys)
    
fracerrcontrol = Err{m,1};
fracerruStim   = Err{m,2};

plot(fracerrcontrol(:),fracerruStim(:),'o','MarkerFaceColor','k', ...
    'MarkerEdgeColor','k','Marker',markers{m})

end

fprintf('\n')

title('Fraction of target misses', 'FontSize', fs);
ylabel({'fraction of';'target misses, uStim'}, 'FontSize',fs);
xlabel('fraction of target misses, no uStim', 'FontSize', fs);

axis tight

% plot unity line
xLim = xlim;
yLim = ylim;
minLim = min([xLim yLim]);
maxLim = max([xLim yLim]);
plot([minLim maxLim],[minLim maxLim],'k-','LineWidth',2)
% make data breathe
gap = 0.04;
xlim([xLim(1)-gap, xLim(2)+gap])
ylim([yLim(1)-gap, yLim(2)+gap])
% make axis equal
xLim = xlim;
yLim = ylim;
minLim = min([xLim yLim]);
maxLim = max([xLim yLim]);
xlim([minLim maxLim])
ylim([minLim maxLim])

ax = gca;
ax.XTick = [0 0.2 0.4];
ax.YTick = [0 0.2 0.4];

if savefig
    saveas(h,[figurepath 'behavior/' savename],'svg')
end


%% Plot delta error fraction as a function of the tuning of the uStime

% In particular, we plot the change in errors (uStim - no uStim) as a 
% function of the angular distance between the vector representing the 
% spatial tuning of the stimulated site and the vector representing the 
% target angle presented to the monkey. 
% Thus, we ask whether the difference between the target angle and the 
% tuning present at the stimulated location was impactful for behavior.

load([statspath 'FRs/' 'uStime_tuning.mat'],'targetdist')

X = [targetdist{1}(:) ; targetdist{2}(:)];
Y = [Err_ang{1,2}(:) - Err_ang{1,1}(:) ; Err_ang{2,2}(:) - Err_ang{2,1}(:)];

idx = isnan(X);
X = X(~idx);
Y = Y(~idx);

[rho,pval]   = corr([X Y]);
strP = {['\rho = ' num2str(rho(1,2),'%.2f')], ['p = ' num2str(pval(1,2),'%.2f')]};

% line fit
R = [X, ones(size(X))];
B = R \ Y;
Yhat = R * B;

% plotting
h=figure; hold on
pos = get(h,'position');
set(h,'position',[pos(1:2) pos(3)*1.05 pos(4)*1])

for m=1:length(monkeys)
    for an = 1:numang
        Errmeancontrol = Err_ang{m,1}(:,:,an);
        ErrmeanuStim   = Err_ang{m,2}(:,:,an);
        
        deltaerr = (ErrmeanuStim - Errmeancontrol);
        targdis  = targetdist{m}(:,:,an);
        
        plot(targdis(:),deltaerr(:),'o','MarkerFaceColor',col_ang(an,:), ...
            'MarkerEdgeColor',col_ang(an,:),'Marker',markers{m})
    end
end

plot(X,Yhat,'k-','Linewidth',2)

title({'uStim electrode tuning direction', '& target misses'}, 'FontSize', fs);
ylabel({'\Delta fraction of target misses', '(uStim - no uStim) '}, 'FontSize',fs);
xlabel('Angle to uStim elec. tuning vector (\circ)', 'FontSize', fs);

axis tight

% make data breathe
xLim = xlim;
yLim = ylim;
gapx = 10; gapy = 0.04;
xlim([xLim(1)-gapx, xLim(2)+gapx])
ylim([yLim(1)-gapy, yLim(2)+gapy])

dim = [.7 .72 .1 .1]; %[x y w h]
annotation('textbox',dim,'String',strP,'Linestyle','none','FontSize',12);

if savefig
    saveas(h,[figurepath 'behavior/' savename '_tuning'],'svg')
end