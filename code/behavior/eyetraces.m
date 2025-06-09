%% Eye traces

%% Paths and specs

run ../addpaths

[datapath, statspath, figurepath] = addpaths;

monkeys    = {'Wa','Sa'};

savefig    = false;
savename   = 'eyetraces_';

% analyse correct trials
code       = 150;
datafolder = 'correct';

% Plotting
[col_ang, markers, dispName] = plottingspecs;

% Plot no uStim condition trials (1), and uStim condition 2 trials
uStim = [1,2];

%% Plot eye traces
   
for m =1:length(monkeys)

monkey = monkeys{m};

[filenames,numuStim,numang,stimorder] = datafiles(monkey);

switch monkey
    case 'Wa'
        top = 1500; % movement artifact threshold
        i   = 1;    % which session to plot
    case 'Sa'
        top = 2500;
        i   = 1;
end

filename = filenames{i};

fprintf('loading file %s \n',filename)
load([datapath datafolder '/' filename],'behavior')

% change order of uStim conditions to have control condition first
behavior(:,:,:,:)  = behavior(stimorder,:,:,:);

trials = size(behavior,4);

% plotting
h=figure; hold on
pos = get(h,'position');
set(h,'position',[pos(1:2) pos(3)*1.5 pos(4)*1])

t  = tiledlayout(1,2);
t.TileSpacing = 'compact';
t.Padding     = 'compact';

% smoothing window size. 
% The standard deviation of the Gaussian filter for smoothdata is fixed to 
% 1/5th of the total window width.
win   = 10; % std = 10ms / 5 = 2ms

clear p
clear pl

pl = gobjects(1,numang);
ax = gobjects(1,length(uStim));

am = 1; % uStim amplitude is fixed

for u = uStim
    ax(u)=nexttile;hold on
    for an = 1:numang
        if ~isempty(behavior{u,am,an,1})
            
        % target angle locations
        trialinfo = behavior{u,am,an,1}.trialinfo;
        switch monkey
            case 'Wa'
                [targetX, targetY] = pol2cart(deg2rad(trialinfo.saccadeDir), ...
                    pix2deg(trialinfo.saccadeLength, trialinfo.screenDistance, trialinfo.pixPerCM));
            case 'Sa'
                [targetX, targetY] = pol2cart(deg2rad(trialinfo.angle), ...
                    pix2deg(trialinfo.distance, trialinfo.screenDistance, trialinfo.pixPerCM));
        end
        
        % plot eye traces
        for tr = 1:trials
            if ~isempty(behavior{u,am,an,tr})
                %find saccade time
                codes = behavior{u,am,an,tr}.codes;
                times = behavior{u,am,an,tr}.times;
                t1 = codes(codes(:,1)==3,2);    %fix off
                t1 = t1(1);
                t2 = codes(codes(:,1)==code,2); % correct/error trials
                if ~isempty(t2)
                bins = find(times>=t1 & times<=t2);
                eyes = behavior{u,am,an,tr}.eyedat(:,bins);
                trialinfo = behavior{u,am,an,tr}.trialinfo;
                if ~any(abs(eyes(1,:))>top) && ~any(abs(eyes(2,:))>top) % exclude trials with movement artifacts
                    % convert to degrees
                    eyes(1:2,:) = eye2deg(eyes(1:2,:),trialinfo);
                    eyes = smoothdata(eyes,2,'gaussian',win);
                    pl(an) = plot(eyes(1,:),eyes(2,:),'Color',col_ang(an,:));
                end
                end
            end
        end       
        % plot target location
        plot(targetX,targetY,'kx','MarkerSize',25,'LineWidth',3);
        end
    end

axis tight
xLim = xlim;
yLim = ylim;
xlim([xLim(1)-0.5, xLim(2)+0.5])
ylim([yLim(1)-0.5, yLim(2)+0.5])

axis off

end

linkaxes(ax,'xy')

if savefig
    saveas(h,[figurepath 'behavior/' savename monkey],'svg')
end

end