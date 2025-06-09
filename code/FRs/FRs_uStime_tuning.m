%% FRs analysis: uStim electrode tuning

%% Paths and specs

run ../addpaths

[datapath, statspath] = addpaths;

monkeys    = {'Wa','Sa'};

% analyse correct trials
datafolder = 'correct';

statspath  = [statspath 'FRs/'];

%% Distance between target location and uStime tuning vector

% trial period to analyse: uStim period
trialperiod = [-300,750]; %ms
[~,uStimp]  = analysisperiod(trialperiod);

% compute angular distance between the target location and the 
% vector representing the tuning of the uStim electrode.
targetdist  = cell(1,length(monkeys));

for m = 1:length(monkeys)

monkey   = monkeys{m};
[filenames,numuStim,numang,stimorder] = datafiles(monkey);


for f=1:length(filenames)
    
filename = filenames{f};

fprintf('loading file %s \n',filename)
load([datapath datafolder '/' filename],'spikerate', 'params')

if strcmp(monkey,'Sa')
    % in monkey Sa, padd short trials with nans
    % although not needed if the trial period selected is the uStim period
    load([datapath datafolder '/' filename],'behavior')
    binWidth  = params.binSize_spikerate*1000;
    spikerate = paddtrialswithnans(spikerate,2,behavior,binWidth);
end

% change order of uStim conditions to have control condition first
spikerate = spikerate(stimorder,:,:,:);
params.cond_uStimChan = params.cond_uStimChan(stimorder);

% remove shorted channels
[chidx, enumoffset] = selectchannels(monkey,datapath,filename);


%% Compute tuning of uStim electrodes

% uStim electrodes (first is 0 = no uStim)
Uchan = params.cond_uStimChan;

% target angles
Ang   = params.cond_targetAngle; % [45 135 225 315]
angle = 1:numang;

% targets in vector space
theta  = deg2rad(str2double(Ang));
Thetav = exp(1i.*theta);

% trials
trials = size(spikerate,4);

for uidx = 1:numuStim-1

   uStime         = str2num(Uchan{uidx+1}) - enumoffset; % str2double doesn't work for 2 uStime
   nonshortedelec = ismember(uStime,chidx); % use non-shorted electrodes
   numuStime      = sum(nonshortedelec);
   uStime         = uStime(nonshortedelec);
   
   if numuStime>0
   
   tuningvec  = zeros(numuStime,length(Thetav));
   
   for i=1:length(uStime)
       
   ci   = uStime(i);
   aidx = 1; %uStim amplitude
   
   %Compute tuning of uStime during uStim period, in no uStim condition
   FRs  = nan(numang,trials,length(uStimp));
    for gidx = angle
        frs  = nan(trials,length(uStimp));
        for r = 1:trials
            if ~isempty(spikerate{1,aidx,gidx,r})
                frs(r,:) = spikerate{1,aidx,gidx,r}(ci,uStimp);
            end
        end
        FRs(gidx,:,:)= frs;
    end     

    % per target angle FRs averaged across trials and time
    FRs = squeeze(nanmean(FRs,[2,3]));
     
    % FR responses to all angles in vector space 
    tuningvec(i,:) = (FRs.*Thetav')./(nansum(FRs,1)+eps);
    
   end
   
   % sum vector responses across angles to compute tuning vector
   % if num uSime > 1, sum tuning vectors across all uStim electrodes
   tuningvec      = sum(tuningvec,'all');
   tuningvec      = tuningvec'/abs(tuningvec);
   
    % compute angular distance between the target location and the 
    % uStime tuning vector.
    distance = zeros(1,length(angle));
    for a = angle
        % compute target angle location in vector space
        angv      = zeros(length(angle),1);
        angv(a)   = 1;
        angvec    = angv.*Thetav';
        angvec    = sum(angvec,'all');
        
        % compute angular distance between the two vectors
        x = [real(tuningvec);imag(tuningvec)];
        y = [real(angvec);imag(angvec)];
        distance(a) = acosd(x'*y);
    end

    targetdist{m}(f,uidx,:) = distance;
    
   else

    targetdist{m}(f,uidx,:) = nan(1,length(Thetav));

   end
end


end

end

save([statspath 'uStime_tuning.mat'],'targetdist')


