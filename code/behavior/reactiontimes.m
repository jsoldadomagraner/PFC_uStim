%% Behavior analysis: reaction times

%% Paths and specs

run ../addpaths

[datapath, statspath, figurepath] = addpaths;

monkeys    = {'Wa','Sa'};

savefig    = false;
savename   = 'RTs';

% analyse correct trials
code       = 150;
datafolder = 'correct';

% Plotting
[col_ang, markers, dispName] = plottingspecs;
   
%% Compute stats

statsfile = [statspath 'behavior/' 'RTstats.mat'];

if ~exist(statsfile,'file')

disp('Computing RT stats...')

RTs     = cell(length(monkeys),2,2);
RTs_ang = cell(length(monkeys),2,2);

for m=1:length(monkeys)

monkey = monkeys{m};

[filenames,numuStim,numang,stimorder] = datafiles(monkey);

fn = length(filenames);

saccangmean = nan(fn,numuStim,numang);
saccangsem  = nan(fn,numuStim,numang);
saccmean    = nan(fn,numuStim,1);
saccsem     = nan(fn,numuStim,1);
RTangmean   = nan(fn,numuStim,numang);
RTangsem    = nan(fn,numuStim,numang);
RTmean      = nan(fn,numuStim,1);
RTsem       = nan(fn,numuStim,1);

for i=1:fn
    
filename = filenames{i};

fprintf('loading file %s \n',filename)
load([datapath datafolder '/' filename],'behavior')

% change order of uStim conditions to have control condition first
behavior(:,:,:,:)  = behavior(stimorder,:,:,:);

trials   = size(behavior,4);

RTsCell  = nan(numuStim,numang,trials);

%% Compute RTs

goCue  = "fixoff"; % in memory guided tasks

am = 1; % uStim amplitude is fixed
for u = 1:numuStim
    for an = 1:numang
        for tr = 1:trials
            if ~isempty(behavior{u,am,an,tr})
                codes     = behavior{u,am,an,tr}.codes;
                if goCue=="fixmove"
                    RTsCell(u,an,tr) = min(codes(codes(:,1)==141,2)) - codes(codes(:,1)==4,2); % 141:saccade, 4:fixmove for visually guided
                elseif goCue == "fixoff"
                    fixoff = codes(codes(:,1)==3,2);
                    if sum(codes(:,1)==141)~=0
                        RTsCell(u,an,tr) = codes(codes(:,1)==141,2) - fixoff(1); % 141:saccade, 3:fixation off for memory guided
                    end
                end              
            end
        end
    end
end

% per angle mean RTs
RTangmean(i,:,:) = nanmean(RTsCell,3);
RTangsem(i,:,:)  = nanstd(RTsCell,0,3)./sqrt(sum(~isnan(RTsCell),3));

% mean RTs averaged across target angle conditions
RTmean(i,:,:) = nanmean(RTsCell,[2,3]);
RTsem(i,:,:)  = nanstd(RTsCell,0,[2,3])./sqrt(sum(~isnan(RTsCell),[2,3]));

end

% RTs averaged across target angles
RTmeanuStim   = RTmean(:,2:end,:);
RTmeancontrol = RTmean(:,1,:);
RTmeancontrol = repmat(RTmeancontrol,[1 numuStim-1 1]);

RTsemuStim   = RTsem(:,2:end,:);
RTsemcontrol = RTsem(:,1,:);
RTsemcontrol = repmat(RTsemcontrol,[1 numuStim-1 1]);

RTs{m,1,1} = RTmeancontrol*1000; %ms
RTs{m,1,2} = RTmeanuStim*1000;
RTs{m,2,1} = RTsemcontrol*1000;
RTs{m,2,2} = RTsemuStim*1000;

% per angle RTs
RTmeanuStim   = RTangmean(:,2:end,:);
RTmeancontrol = RTangmean(:,1,:);
RTmeancontrol = repmat(RTmeancontrol,[1 numuStim-1 1]);

RTsemuStim   = RTangsem(:,2:end,:);
RTsemcontrol = RTangsem(:,1,:);
RTsemcontrol = repmat(RTsemcontrol,[1 numuStim-1 1]);

RTs_ang{m,1,1} = RTmeancontrol*1000; %ms
RTs_ang{m,1,2} = RTmeanuStim*1000;
RTs_ang{m,2,1} = RTsemcontrol*1000;
RTs_ang{m,2,2} = RTsemuStim*1000;

end
    save(statsfile,'RTs','RTs_ang')
else
    load(statsfile)
end

%% Test whether uStim and no uStim RTs are statistically different

X = [RTs{1,1,1}(:) ; RTs{2,1,1}(:)]; % no uStim
Y = [RTs{1,1,2}(:) ; RTs{2,1,2}(:)]; % uStim

% Two-tailed paired t-test
% Test of the null hypothesis that the two samples come from distributions 
% with equal means, against the alternative that they do not.
[hy,p] = ttest(X,Y);
fprintf('t-test p=%.3f \n \n',p)

fprintf('RT no uStim = %.0f ± %.0f \n',nanmean(X),nanstd(X))
fprintf('RT uStim = %.0f ± %.0f \n',nanmean(Y),nanstd(Y))

fprintf('RT difference = %.2f \n \n',nanmean(X)-nanmean(Y))

% Linear regression
disp('linear regression')
% Compute p-values of the full model based on the F-statistic (H0: reduced 
% model, just intercept; HA: full model, intercept + slope).
R = [X, ones(size(X))];
[b,~,~,~,stats] = regress(Y,R);
fprintf('slope=%.2f, offset=%.2f, p=%.2e \n', b(1), b(2), stats(3))


%% Plot mean RTs in uStim vs. no uStim

h=figure; hold on
pos = get(h,'position');
set(h,'position',[pos(1:2) pos(3)*1 pos(4)*1])

% font size
fs = 20;

for m=1:length(monkeys)
    
RTmeancontrol = RTs{m,1,1};
RTmeanuStim   = RTs{m,1,2};

RTsemcontrol = RTs{m,2,1};
RTsemuStim   = RTs{m,2,2};

plot(RTmeancontrol(:),RTmeanuStim(:),'o','MarkerFaceColor','k', ...
    'MarkerEdgeColor','k','Marker',markers{m})

end

title("Mean reaction time", 'FontSize', fs);
ylabel('RT, uStim (ms)', 'FontSize',fs);
xlabel('RT, no uStim (ms)', 'FontSize', fs);

axis tight

% plot unity line
xLim = xlim;
yLim = ylim;
minLim = min([xLim yLim]);
maxLim = max([xLim yLim]);
plot([minLim maxLim],[minLim maxLim],'k-','LineWidth',2)
% make data breathe
gap = 3;
xlim([xLim(1)-gap, xLim(2)+gap])
ylim([yLim(1)-gap, yLim(2)+gap])
% make axis equal
xLim = xlim;
yLim = ylim;
minLim = min([xLim yLim]);
maxLim = max([xLim yLim]);
xlim([minLim maxLim])
ylim([minLim maxLim])

if savefig
    saveas(h,[figurepath 'behavior/' savename],'svg')
end


%% Plot delta RTs as a function of the tuning of the uStime

% In particular, we plot the change in RTs (uStim - no uStim) as a 
% function of the angular distance between the vector representing the 
% spatial tuning of the stimulated site and the vector representing the 
% target angle presented to the monkey. 
% Thus, we ask whether the difference between the target angle and the 
% tuning present at the stimulated location was impactful for behavior.

load([statspath 'FRs/' 'uStime_tuning.mat'],'targetdist')

X = [targetdist{1}(:) ; targetdist{2}(:)];
Y = [RTs_ang{1,1,2}(:) - RTs_ang{1,1,1}(:) ; RTs_ang{2,1,2}(:) - RTs_ang{2,1,1}(:)];

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

[~,numuStim, numang] = size(RTs_ang{m});

for m=1:length(monkeys)
    for an = 1:numang
        RTmeancontrol = RTs_ang{m,1,1}(:,:,an);
        RTmeanuStim   = RTs_ang{m,1,2}(:,:,an);
        
        deltaRT = (RTmeanuStim - RTmeancontrol);
        targdis = targetdist{m}(:,:,an);
        
        plot(targdis(:),deltaRT(:),'o','MarkerFaceColor',col_ang(an,:), ...
            'MarkerEdgeColor',col_ang(an,:),'Marker',markers{m})
    end
end

plot(X,Yhat,'k-','Linewidth',2)

title({'uStim electrode tuning direction', '& reaction time'}, 'FontSize', fs);
ylabel({'\Delta RT (ms)', '(uStim - no uStim)'}, 'FontSize',fs);
xlabel('Angle to uStime tuning vector (\circ)', 'FontSize', fs);

axis tight

% make data breathe
xLim = xlim;
yLim = ylim;
gapx = 10; gapy = 3;
xlim([xLim(1)-gapx, xLim(2)+gapx])
ylim([yLim(1)-gapy, yLim(2)+gapy])

dim = [.7 .72 .1 .1]; %[x y w h]
annotation('textbox',dim,'String',strP,'Linestyle','none','FontSize',12);

if savefig
    saveas(h,[figurepath 'behavior/' savename '_tuning'],'svg')
end