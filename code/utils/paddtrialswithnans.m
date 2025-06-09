% dat = 1 for spiketrains
% dat = 2 for spikerate

function spikes = paddtrialswithnans(spikes,dat,behavior,binWidth)

[numuStim,numamp,numang,trials] = size(behavior);

goCue = 'fixoff';  % memory guided tasks
% goCue = fixmove'; % visually guided tasks in Sa

switch goCue
    case 'fixoff'
        gocode = 3;
    case 'fixmove'
        gocode = 4;
end

for u = 1:numuStim
    for am = 1:numamp
        for an = 1:numang
            for tr = 1:trials
                if ~isempty(spikes{u,am,an,tr}) && ~isempty(behavior{u,am,an,tr})
                    codes  = behavior{u,am,an,tr}.codes;
                    gotime = codes(codes(:,1)==gocode,3);
                    if dat==1 %spiketrain
                        spikes{u,am,an,tr}(:,gotime:end)=nan;
                    else %spikerate
                        T = codes(end,3);
                        binned   = 1:binWidth:T;
                        [~,bins] = find(binned>gotime);
                        spikes{u,am,an,tr}(:,bins(1):end)=nan;
                    end
                end
            end
        end
    end
end
