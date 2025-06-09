% Plot dPCA example trajectories

monkey   = 'Wa';
session  = 8;

savefig  = false;

% plot no uStim condition and uStim condition 3: {'0','12','41','91'};
myconditions = {'0','91'}; 

% plot target angles 45 and 315: [45 135 225 315]
angorder     = [1 4 2 3];
numangplot   = 2;

% if true, plot no uStim trajectories only, and for all angles
nouSallang = false;

if nouSallang
    myconditions  = {'0'}; 
    numangplot    = 4; 
    nouSallanglab = '_nouSallang'; 
end

[datapath, statspath, figurepath] = addpaths;

% analyse correct trials
datafolder = 'correct';

numdims    = 4; % mean optimal latent dimensionality

% numbering offset for electrodes on different banks (Wa=0, S=256)
enumoffset   = 0;
% control condition is 0 uStime
cntrol       = '0'; % for Wa '0', for Sa '0  0'

%% Load data

[filenames,~,~,stimorder] = datafiles(monkey);

filename = filenames{session};

fprintf('file %s \n', filenames{session})

% trial period to analyse (pre and post uStim)
trialperiod = [-300,800]; %ms

load([datapath datafolder '/' filename],'spiketrain','params')

switch monkey
    case 'Sa'
        % in monkey Sa, padd short trials with nans
        load([datapath datafolder '/' filename],'behavior')
        spiketrain = paddtrialswithnans(spiketrain,1,behavior);
end

% change order of uStim conditions to have control condition first
spiketrain = spiketrain(stimorder,:,:,:);
params.cond_uStimChan = params.cond_uStimChan(stimorder);

% remove shorted channels
chidx    = selectchannels(monkey,datapath,filename);

binWidth = params.binSize_spikerate*1000;

% compute binned spike counts (pre-uStim, post-uStim and all times)
[spikeTrainsBin01, spikeTrainsBin02, spikeTrainsBin0] = ...
    binspiketrains(spiketrain,chidx,binWidth,trialperiod);

firingRates01 = spikeTrainsBin01/(binWidth/1000);
firingRatesAverage01 = nanmean(firingRates01,5);

firingRates02 = spikeTrainsBin02/(binWidth/1000);
firingRatesAverage02 = nanmean(firingRates02,5);

firingRates0  = spikeTrainsBin0/(binWidth/1000);

firingRatesAverage = cat(4,firingRatesAverage01,firingRatesAverage02);
T = size(firingRatesAverage,4);

%% Fit dPCA model to mean firing rates

% dPCA marginalizations
% {'Stimulus', 'Angle', 'Condition-independent(Time)', 'S/A Interaction'};
combinedParams = {{1, [1 3]}, {2, [2 3]}, {3}, {[1 2], [1 2 3]}};
numComp = [numdims numdims numdims numdims];
 
[W,V,whichMarg] = dpca(firingRatesAverage, numComp, ...
    'combinedParams', combinedParams);

%% Single trial data projections in the dPCA subspace

% Find memory dimensions (marginalization 2)
marg = 2;
mg   = find(whichMarg==marg);
mg   = mg(1:numdims);
W    = W(:,mg);

Yplotstd{1} = firingRates01;
Yplotstd{2} = firingRates02;

for d=1:2
    dataDim = size(Yplotstd{d});
    X     = reshape(Yplotstd{d},dataDim(1),[])';
    Xcen  = bsxfun(@minus, X, nanmean(X));
    Z     = Xcen * W;
    Zfull = reshape(Z', [numdims dataDim(2:end)]);
    Xm{d} = Zfull;
end

% Data that includes all trial times, for control condition
Yplotstd0 = firingRates0;

dataDim = size(Yplotstd0);
X     = reshape(Yplotstd0,dataDim(1),[])';
Xcen  = bsxfun(@minus, X, nanmean(X));
Z     = Xcen * W;
Zfull = reshape(Z', [numdims dataDim(2:end)]);
Xm0   = Zfull;

% Zero FRs vector
Xv0 = zeros(dataDim(1),1)';
Xv0cen = bsxfun(@minus, Xv0, nanmean(X));
Zv0 = Xv0cen * W;
v0  = Zv0';


%% Plotting specs

smoothlndata = true;

Dims = [1 2]; % latents dimensions to plot

Z    = 1.96; % confidence interval

[colang, ~,~,stimcol] = plottingspecs;

colang = colang(angorder,:);

alpha  = 1;
alphal = 0.2;
mks    = 6;
lw1    = 2;
lw2    = 3;
mks1   = 15;
mks2   = 25;

% flip axes to match target angle order
flip   = [1 -1];

%% Format dPCA data

uStime = params.cond_uStimChan(stimorder);
angle  = params.cond_targetAngle;

Xm0   = permute(Xm0,[1 4 5 2 3]);
Xm{1} = permute(Xm{1},[1 4 5 2 3]); % {d}(ln,S,D,T01,E); d=pre/post
Xm{2} = permute(Xm{2},[1 4 5 2 3]); % (ln,S,D,T01,E)

if smoothlndata
    xmpre  = smoothdata(Xm{1},2,'gaussian',6);
    xmpost = smoothdata(Xm{2},2,'gaussian',6);
    xm0prepost = smoothdata(Xm0,2,'gaussian',6);
else
    xmpre  = Xm{1}; 
    xmpost = Xm{2};
    xm0prepost = Xm0; 
end
xm = cat(2,xmpre,xmpost); %(ln,ti,tr,s,a)
xm0 = xm0prepost;         %(ln,ti0,tr,s,a)

trialperiod    = [-300,750]; %ms
[~,~,~,tpre,~] = analysisperiod(trialperiod);
taupre = 1:tpre;

tau     = 1:size(xm,2);
taupost = tau(tpre+1:end);
meand   = [3,4];

for c3 = 1:length(angle)
    xmang{c3}  = squeeze(xm(:,:,:,:,c3));
    xmang0{c3} = squeeze(xm0(:,:,:,:,c3));
end
tr = size(xm,3);


%% Plot latent trajectories

if nouSallang 
    xx= 0.8; yy=1;
    a = 1; b = 1;
else
    xx=1.8; yy=1;
    a = 1 ; b = 2;
end

h1=figure;
pos = get(h1,'position');
set(h1,'position',[pos(1:2) pos(3)*xx pos(4)*yy])

tl = gobjects(a,b);
p  = tiledlayout(a,b);
p.Padding     = 'compact';

for j=1:length(uStime)

    ustime = str2num(uStime{j});
    idx    = ustime~=0;
    ustime(idx) = ustime(idx)-enumoffset;

    if ismember(num2str(ustime),myconditions)

    for r = 1:numangplot
        dim1 = Dims(1); dim2 = Dims(2);
        xm   = xmang{angorder(r)};
        xm0  = xmang0{angorder(r)};

        % plot latents
        if nouSallang
            tl(1)=nexttile(1); hold on;
        else
            tl(r)=nexttile(r); hold on;
        end
        
        % plot no uStim trajectory
        if strcmp(num2str(ustime),cntrol)
            tau  = taupre(end):taupost(end);
            x1   = xm0(dim1,tau,:,j)*flip(1);
            x2   = xm0(dim2,tau,:,j)*flip(2);
            x1m  = nanmean(x1,meand);
            x2m  = nanmean(x2,meand);
            x1std = Z*nanstd(x1,0,meand)/sqrt(tr);
            x2std = Z*nanstd(x2,0,meand)/sqrt(tr);
            % plot ellipse
            Dot = [1 size(x1m,2)];
            for d = 1:length(Dot)
                dot = Dot(d);
                pd = plot(x1m(dot),x2m(dot),'.');
                pd.MarkerSize = mks2;
                pd.Color      = [stimcol(1,:) alpha];
                h = plotEllipses([x1m(dot) x2m(dot)],[x1std(dot) x2std(dot)]);
                h.FaceColor = [stimcol(1,:) 0.3];
                h.EdgeColor = colang(r,:);
                h.LineWidth = 2;
            end
            %plot trajectory
            p1=plot(x1m,x2m,'.-');
            p1.LineWidth  = lw1;
            p1.LineStyle  = '-';
            p1.MarkerSize = mks1;
            p1.Color      = [stimcol(1,:) alpha];
            % plot arrow
            a = 5; b = 6; % arrow placement
            alphar = 1;   % Size of arrow head relative to the length of the vector
            betar  = 1;   % Width of the base of the arrow head relative to the length
            p4 = [x1m(a) x2m(a)]; % x(4): 4th x value point of a trajectory
            p5 = [x1m(b) x2m(b)]; % y(4): 4th y value point of a trajectory
            diff45 = p5 - p4;
            diff45 = diff45/norm(diff45); % normalize the vector length
            hu = [p5(1)-alphar*(diff45(1)+betar*(diff45(2)+eps)); p5(1); p5(1)-alphar*(diff45(1)-betar*(diff45(2)+eps))];
            hv = [p5(2)-alphar*(diff45(2)-betar*(diff45(1)+eps)); p5(2); p5(2)-alphar*(diff45(2)+betar*(diff45(1)+eps))];
            % Plot arrow head
            plot(hu(:),hv(:), 'color', stimcol(1,:),'LineWidth',lw1)
        end

        % plot uStim trajectory
        if ~strcmp(num2str(ustime),cntrol)
        for tt = 1:2
        if tt==1
            tau=[taupre(end) taupost(1)]; 
            lnty='--'; lw = lw2; mks = mks1;
            a = 1; b = 2; % arrow placement
        else
            tau=taupost; lnty = '-'; 
            lw = lw1; mks = mks1;
            a = 2; b = 3; % arrow placement
        end
        x1   = xm(dim1,tau,:,j)*flip(1);
        x2   = xm(dim2,tau,:,j)*flip(2);
        x1m  = nanmean(x1,meand);
        x2m  = nanmean(x2,meand);
        x1std = Z*nanstd(x1,0,meand)/sqrt(tr);
        x2std = Z*nanstd(x2,0,meand)/sqrt(tr);

        % plot ellipse
        Dot = [1 size(x1,2)];
        if tt==1; d=1;end
        if tt==2; d=2;end
        dot = Dot(d);
        pd = plot(x1m(dot),x2m(dot),'.');
        pd.MarkerSize = mks2;
        pd.Color      = [stimcol(2,:) alpha];
        h = plotEllipses([x1m(dot) x2m(dot)],[x1std(dot) x2std(dot)]);
        h.FaceColor = [stimcol(2,:) 0.3];
        h.EdgeColor = colang(r,:);
        h.LineWidth = 2;

        % plot trajectory
        p1=plot(x1m,x2m,'.-');
        p1.LineWidth  = lw;
        p1.LineStyle  = lnty;
        p1.MarkerSize = mks;
        p1.Color      = [stimcol(2,:) alpha];

        if tt==2
        pd = plot(x1m(dot),x2m(dot),'.');
        pd.MarkerSize = mks2;
        pd.Color      = [stimcol(2,:) alpha];
        end

        % plot arrow
        alphar = 1; %0.05  % Size of arrow head relative to the length of the vector
        betar  = 1; %0.4 % Width of the base of the arrow head relative to the length
        p4 = [x1m(a) x2m(a)]; % x(4): 4th x value point of a trajectory
        p5 = [x1m(b) x2m(b)]; % y(4): 4th y value point of a trajectory
        diff45 = p5 - p4;
        diff45 = diff45/norm(diff45); % normalize the vector length
        hu = [p5(1)-alphar*(diff45(1)+betar*(diff45(2)+eps)); p5(1); p5(1)-alphar*(diff45(1)-betar*(diff45(2)+eps))];
        hv = [p5(2)-alphar*(diff45(2)-betar*(diff45(1)+eps)); p5(2); p5(2)-alphar*(diff45(2)+betar*(diff45(1)+eps))];
        % Plot arrow head
        plot(hu(:),hv(:), 'color', stimcol(2,:), 'LineWidth', lw)
        end
        end
        
        if nouSallang
            tl(1).YTick=[]; tl(1).XTick=[];
        else
            tl(r).YTick=[]; tl(r).XTick=[];
        end

    % plot vector zero
    plot(v0(1),v0(2),'x','MarkerSize',25,'Color','k','LineWidth',4);

    end
    end
end 

linkaxes(tl,'xy')

if nouSallang; numangplot = 1;end

for r = 1:numangplot 
    nexttile(r)
    yLim = ylim; xLim = xlim;
    l1 = plot([0 0],yLim,'k','LineWidth',1);
    l2 = plot(xLim,[0 0],'k','LineWidth',1);
    l1.Color = [0 0 0 alphal];
    l2.Color = [0 0 0 alphal];
    set(gca,'Visible','off')
end

if savefig
    saveas(h1,[figurepath 'subspaces/' 'dPCAtrajectories' nouSallanglab],'svg')
end

%% Plot no-uStim data tpre to tpost for all angles

if ~nouSallang
    
Tau = [tpre,taupost(1)];
for co=1:length(uStime)
    ustime = str2num(uStime{j});
    idx    = ustime~=0;
    ustime(idx) = ustime(idx)-enumoffset;
    
    if ismember(num2str(ustime),myconditions{end})
        for a=1:length(angle)
            for d=1:length(Dims)
                xmean(a,d,:) = squeeze(nanmean(xmang{a}(Dims(d),Tau,:,co),3))*flip(d);
                xci(a,d,:)   = (Z*squeeze(nanstd(xmang{a}(Dims(d),Tau,:,co),0,3))/sqrt(tr));
            end
        end
    end
    
end

% plotting specs
[colang, ~,~, stimcol] = plottingspecs;
colplane = stimcol(3,:);

lw2    = 5;
lw3    = 4;
alpha  = 0.3;
alphap = 0.05;
alphal = 0.2;
offset = 2;

% Arrow head size
Alpha  = 2.2;
Beta   = 1;

h2 = figure; hold on
xx =0.8; yy=0.9;
pos = get(h2,'position');
set(h2,'position',[pos(1:2) pos(3)*xx pos(4)*yy])

% make dummy figure to know size of colored planes
for a=1:length(angle)
    
    Xpre  = [xmean(a,1,1), xmean(a,2,1)];
    Xpost = [xmean(a,1,2), xmean(a,2,2)];
    CIpre = [xci(a,1,1), xci(a,2,1)];
    
    e = plotEllipses(Xpre,CIpre);
    e.FaceColor = [1 1 1];
    e.EdgeColor = [1 1 1];
    e.LineWidth = lw3;
    
    x1 = Xpre(1); x2 = Xpost(1);
    y1 = Xpre(2); y2 = Xpost(2);
    l  = plot([x1 x2],[y1 y2],'--');
    l.LineWidth = lw2;
    l.Color     = [1 1 1];
    
    % plot arrow
    alphar = Alpha; % Size of arrow head relative to the length of the vector
    betar  = Beta;  % Width of the base of the arrow head relative to the length
    p4 = [x1 y1]; % x(4): 4th x value point of a trajectory
    p5 = [x2 y2]; % y(4): 4th y value point of a trajectory
    diff45 = p5 - p4;
    diff45 = diff45/norm(diff45); % normalize the vector length
    hu = [p5(1)-alphar*(diff45(1)+betar*(diff45(2)+eps)); p5(1); p5(1)-alphar*(diff45(1)-betar*(diff45(2)+eps))];
    hv = [p5(2)-alphar*(diff45(2)-betar*(diff45(1)+eps)); p5(2); p5(2)-alphar*(diff45(2)+betar*(diff45(1)+eps))];
    % Plot arrow head
    plot(hu(:),hv(:), 'color', [1 1 1], 'LineWidth', lw2)
    
end

ax=gca;
ax.XTick = []; ax.YTick = [];
axis tight
off = offset;
xlim([ax.XLim(1)-off ax.XLim(2)+off])
ylim([ax.YLim(1)-off ax.YLim(2)+off])

yLim = ylim; xLim = xlim;
l1 = plot([0 0],yLim,'k','LineWidth',1);
l2 = plot(xLim,[0 0],'k','LineWidth',1);
l1.Color = [0 0 0 alphal];
l2.Color = [0 0 0 alphal];


% Plot planes according to figure size
XL = get(gca, 'XLim');
YL = get(gca, 'YLim');
patch([XL(1), XL(2), XL(2), XL(1)], [YL(1), YL(1), YL(2), YL(2)], ...
    [0 0 0 0], 'FaceColor', colplane,'FaceAlpha',alphap,...
    'EdgeColor', colplane,'LineWidth',2);

% Plot data on top of planes
for a=1:length(angle)
    
    Xpre  = [xmean(a,1,1), xmean(a,2,1)];
    Xpost = [xmean(a,1,2), xmean(a,2,2)];
    CIpre = [xci(a,1,1), xci(a,2,1)];

    e = plotEllipses(Xpre,CIpre);
    e.FaceColor = [stimcol(2,:) alpha];
    e.EdgeColor = [colang(a,:) 1];
    e.LineWidth = lw3;

    x1 = Xpre(1); x2 = Xpost(1);
    y1 = Xpre(2); y2 = Xpost(2);
    l  = plot([x1 x2],[y1 y2],'--');
    l.LineWidth = lw2;
    l.Color     = stimcol(2,:);
    
    % plot arrow
    alphar = Alpha; % Size of arrow head relative to the length of the vector
    betar  = Beta;  % Width of the base of the arrow head relative to the length
    p4 = [x1 y1]; % x(4): 4th x value point of a trajectory
    p5 = [x2 y2]; % y(4): 4th y value point of a trajectory
    diff45 = p5 - p4;
    diff45 = diff45/norm(diff45); % normalize the vector length
    hu = [p5(1)-alphar*(diff45(1)+betar*(diff45(2)+eps)); p5(1); p5(1)-alphar*(diff45(1)-betar*(diff45(2)+eps))];
    hv = [p5(2)-alphar*(diff45(2)-betar*(diff45(1)+eps)); p5(2); p5(2)-alphar*(diff45(2)+betar*(diff45(1)+eps))];
    % Plot arrow head
    plot(hu(:),hv(:), 'color', stimcol(2,:), 'LineWidth', lw2)
    
end

if savefig
    saveas(h2,[figurepath 'subspaces/' 'dPCAtrajectories_uSallang'],'svg')
end

end