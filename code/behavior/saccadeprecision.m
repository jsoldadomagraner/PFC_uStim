%% Behavior analysis: saccade precision

%% Paths and specs

run ../addpaths

[datapath, statspath, figurepath] = addpaths;

monkeys    = {'Wa','Sa'};

savefig    = false;
savename   = 'saccade';

% analyse correct trials
code       = 150;
datafolder = 'correct';

% Plotting
[col_ang, markers, dispName] = plottingspecs;

%% Compute stats

statsfile = [statspath 'behavior/' 'saccadestats.mat'];

if ~exist(statsfile,'file')

disp('Computing saccade stats...')

Sacc     = cell(length(monkeys),2,2);
Sacc_ang = cell(length(monkeys),2,2);

for m =1:length(monkeys)

monkey = monkeys{m};

[filenames,numuStim,numang,stimorder] = datafiles(monkey);

fn = length(filenames);

saccangmean = nan(fn,numuStim,numang);
saccangsem  = nan(fn,numuStim,numang);
saccmean    = nan(fn,numuStim,1);
saccsem     = nan(fn,numuStim,1);

for i=1:fn
    
filename = filenames{i};

fprintf('loading file %s \n',filename)
load([datapath datafolder '/' filename],'behavior')

% change order of uStim conditions to have control condition first
behavior(:,:,:,:)  = behavior(stimorder,:,:,:);

trials   = size(behavior,4);

SaccCell = nan(numuStim,numang,trials);

%% Compute saccade precision (saccade end distance to target)

am = 1; % uStim amplitude is fixed
for u = 1:numuStim
    for an = 1:numang
        for tr = 1:trials
            if ~isempty(behavior{u,am,an,tr})
                codes     = behavior{u,am,an,tr}.codes;
                trialinfo = behavior{u,am,an,tr}.trialinfo;
                eyes      = behavior{u,am,an,tr}.eyedat;
                nstime    = behavior{u,am,an,tr}.times;
                              
                eyeXYZ  = eye2deg(eyes,trialinfo);
                eyeXTmp = eyeXYZ(1,:);
                eyeYTmp = eyeXYZ(2,:);

                switch monkey
                    case 'Wa'
                        [targetX, targetY] = pol2cart(deg2rad(trialinfo.saccadeDir), ...
                            pix2deg(trialinfo.saccadeLength, trialinfo.screenDistance, trialinfo.pixPerCM));
                    case 'Sa'
                        [targetX, targetY] = pol2cart(deg2rad(trialinfo.angle), ...
                            pix2deg(trialinfo.distance, trialinfo.screenDistance, trialinfo.pixPerCM));
                end

                sacEndTime   = codes(codes(:,1)==code,2);
                sacEndNSTIME = min(nstime(nstime>sacEndTime));
                sacEndX = eyeXTmp(nstime==sacEndNSTIME);
                sacEndY = eyeYTmp(nstime==sacEndNSTIME);

                precision = sqrt((sacEndX-targetX)^2+(sacEndY-targetY)^2);
                SaccCell(u,an,tr) = precision;
            end
        end
    end
end

% per angle mean saccade distance
saccangmean(i,:,:) = nanmean(SaccCell,3);
saccangsem(i,:,:)  = nanstd(SaccCell,0,3)./sqrt(sum(~isnan(SaccCell),3));

% mean saccade distance across target angle conditions
saccmean(i,:,:) = nanmean(SaccCell,[2,3]);
saccsem(i,:,:)  = nanstd(SaccCell,0,[2,3])./sqrt(sum(~isnan(SaccCell),[2,3]));

end

% saccade distance averaged across target angles
saccmeanuStim   = saccmean(:,2:end,:);
saccmeancontrol = saccmean(:,1,:);
saccmeancontrol = repmat(saccmeancontrol,[1 numuStim-1 1]);

saccsemuStim   = saccsem(:,2:end,:);
saccsemcontrol = saccsem(:,1,:);
saccsemcontrol = repmat(saccsemcontrol,[1 numuStim-1 1]);

Sacc{m,1,1} = saccmeancontrol; %degrees
Sacc{m,1,2} = saccmeanuStim;
Sacc{m,2,1} = saccsemcontrol;
Sacc{m,2,2} = saccsemuStim;


% saccade distance per angle
saccmeanuStim   = saccangmean(:,2:end,:);
saccmeancontrol = saccangmean(:,1,:);
saccmeancontrol = repmat(saccmeancontrol,[1 numuStim-1 1]);

saccsemuStim   = saccangsem(:,2:end,:);
saccsemcontrol = saccangsem(:,1,:);
saccsemcontrol = repmat(saccsemcontrol,[1 numuStim-1 1]);

Sacc_ang{m,1,1} = saccmeancontrol; %degrees
Sacc_ang{m,1,2} = saccmeanuStim;
Sacc_ang{m,2,1} = saccsemcontrol;
Sacc_ang{m,2,2} = saccsemuStim;

end
    save(statsfile,'Sacc','Sacc_ang')
else
    load(statsfile)
end

%% Test whether uStim and no uStim saccades are statistically different

X = [Sacc{1,1,1}(:) ; Sacc{2,1,1}(:)]; % no uStim
Y = [Sacc{1,1,2}(:) ; Sacc{2,1,2}(:)]; % uStim

% Two-tailed paired t-test
% Test of the null hypothesis that the two samples come from distributions 
% with equal means, against the alternative that they do not.
[hy,p] = ttest(X,Y);
fprintf('t-test p=%.2f \n \n',p)

% Linear regression
disp('linear regression')
% Compute p-values of the full model based on the F-statistic (H0: reduced 
% model, just intercept; HA: full model, intercept + slope).
R = [X, ones(size(X))];
[b,~,~,~,stats] = regress(Y,R);
fprintf('slope=%.2f, offset=%.2f, p=%.2e \n', b(1), b(2), stats(3))


%% Plot mean saccade precision in uStim vs. no uStim

h=figure; hold on
pos = get(h,'position');
set(h,'position',[pos(1:2) pos(3)*1 pos(4)*1])

% font size
fs = 20;

for m=1:length(monkeys)
    
saccmeancontrol = Sacc{m,1,1};
saccmeanuStim   = Sacc{m,1,2};

saccsemcontrol = Sacc{m,2,1};
saccsemuStim   = Sacc{m,2,2};

    plot(saccmeancontrol(:),saccmeanuStim(:),'o','MarkerFaceColor','k', ...
        'MarkerEdgeColor','k','Marker',markers{m})
    
end

title("Mean saccade precision", 'FontSize', fs);
ylabel({'distance from target,';'uStim (\circ)'}, 'FontSize',fs);
xlabel('distance from target, no uStim (\circ)', 'FontSize', fs);

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
ax.XTick = [1 1.2 1.4];
ax.YTick = [1 1.2 1.4];

% dummy plot for legend
l1 = plot(nan,nan,'o','MarkerFaceColor','k', ...
    'MarkerEdgeColor','k','Marker',markers{1});
l2 = plot(nan,nan,'o','MarkerFaceColor','k', ...
    'MarkerEdgeColor','k','Marker',markers{2});
legend([l1,l2],{'monkey W','monkey S'},'box','off','Location','NorthWest')


if savefig
    saveas(h,[figurepath 'behavior/' savename],'svg')
end


%% Plot delta saccade accuracy as a function of the tuning of the uStime

% In particular, we plot the change in accuracy (uStim - no uStim) as a 
% function of the angular distance between the vector representing the 
% spatial tuning of the stimulated site and the vector representing the 
% target angle presented to the monkey. 
% Thus, we ask whether the difference between the target angle and the 
% tuning present at the stimulated location was impactful for behavior.

load([statspath 'FRs/' 'uStime_tuning.mat'],'targetdist')

X = [targetdist{1}(:) ; targetdist{2}(:)];
Y = [Sacc_ang{1,1,2}(:) - Sacc_ang{1,1,1}(:) ; Sacc_ang{2,1,2}(:) - Sacc_ang{2,1,1}(:)];

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
set(h,'position',[pos(1:2) pos(3)*1 pos(4)*1])

[~,numuStim, numang] = size(Sacc_ang{m});

for m=1:length(monkeys)
    for an = 1:numang
        Saccmeancontrol = Sacc_ang{m,1,1}(:,:,an);
        SaccmeanuStim   = Sacc_ang{m,1,2}(:,:,an);
        
        deltaprec = (SaccmeanuStim - Saccmeancontrol);
        targdis   = targetdist{m}(:,:,an);
        
        plot(targdis(:),deltaprec(:),'o','MarkerFaceColor',col_ang(an,:), ...
            'MarkerEdgeColor',col_ang(an,:),'Marker',markers{m})
    end
end

plot(X,Yhat,'k-','Linewidth',2)

title({'uStim electrode tuning direction',' & saccade precision'}, 'FontSize', fs);
ylabel({'\Delta distance from target (\circ)', '(uStim - no uStim)'}, 'FontSize',fs);
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