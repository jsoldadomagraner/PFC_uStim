function [spikeTrainsBin01, spikeTrainsBin02, spikeTrainsBin0] = ...
    binspiketrains(spiketrain,chidx,binWidth,trialperiod)

% trial period to analyse
trialtimes = -800:1:800; %ms

% trial period to analyse 
trstart   = trialperiod(1);
trend     = trialperiod(2);

% uStim condition (pre-uStim and post-uStim)
tperiod01 = find(trialtimes>=trstart & trialtimes<0);
tau1      = length(tperiod01);
tperiod02 = find(trialtimes>=200 & trialtimes<=trend);
tau2      = length(tperiod02);

% control condition
tperiod0  = find(trialtimes>=trstart & trialtimes<=trend);
tau       = length(tperiod0);

T01   = floor(tau1/binWidth);
T02   = floor(tau2/binWidth);
T0    = floor(tau/binWidth);

numch = length(chidx);

S = size(spiketrain,1);   % number of uStim patterns
D = size(spiketrain,3);   % number of angles
E = size(spiketrain,4)-1; % trials (-1 to have the same per condition)
N = numch;

spikeTrainsBin01 = nan(N,S,D,T01,E);
spikeTrainsBin02 = nan(N,S,D,T02,E);
spikeTrainsBin0  = nan(N,S,D,T0,E);

aidx  = 1; % uStima
for si=1:S % uStime
    for di=1:D % angle
        for ei=1:E % trials
        if ~isempty(spiketrain{si,aidx,di,ei})
            spiketrain{si,aidx,di,ei} = spiketrain{si,aidx,di,ei}(chidx,:,:);
            spikespre = spiketrain{si,aidx,di,ei}(:,tperiod01);
            for t = 1:T01
                iStart = binWidth * (t-1) + 1;
                iEnd   = binWidth * t;
                spikeTrainsBin01(:,si,di,t,ei) = sum(spikespre(:,iStart:iEnd),2);
            end
            spikespost = spiketrain{si,aidx,di,ei}(:,tperiod02);
            for t = 1:T02
                iStart = binWidth * (t-1) + 1;
                iEnd   = binWidth * t;
                spikeTrainsBin02(:,si,di,t,ei) = sum(spikespost(:,iStart:iEnd),2);
            end
            spikesall = spiketrain{si,aidx,di,ei}(:,tperiod0);
            for t = 1:T0
                iStart = binWidth * (t-1) + 1;
                iEnd   = binWidth * t;
                spikeTrainsBin0(:,si,di,t,ei) = sum(spikesall(:,iStart:iEnd),2);
            end
        end
        end
    end
end

end