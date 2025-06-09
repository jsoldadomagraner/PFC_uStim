function dpca_plot_default_uStim(data, time, yspan, explVar, compNum, ...
    events, signif, marg, whichplot)

% Modify this function to adjust how components are plotted.
%
% Parameters are as follows:
%   data      - data matrix, size(data,1)=1 because it's only one component
%   time      - time axis
%   yspan     - y-axis spab
%   explVar   - variance of this component
%   compNum   - component number
%   events    - time events to be marked on the time axis
%   signif    - marks time-point where component is significant
%   marg      - marginalization number

% Soldado-Magraner J. made modifications for custom plotting, 2025.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% displaying legend
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(data, 'legend')
    
    % if there is only time and no other parameter - do nothing
    if length(time) == 2
        return

    % if there is one parameter
    elseif length(time) == 3
        numOfStimuli = time(2); % time is used to pass size(data) for legend
        colors = lines(numOfStimuli);
        hold on
        
        for f = 1:numOfStimuli
            plot([0.5 1], [f f], 'color', colors(f,:), 'LineWidth', 2)
            text(1.2, f, ['Stimulus ' num2str(f)])
        end
        axis([0 3 -1 1.5+numOfStimuli])
        set(gca, 'XTick', [])
        set(gca, 'YTick', [])
        set(gca,'Visible','off')
        return

    % two parameters: stimulus and decision (decision can only have two
    % values)
    elseif length(time) == 4 && time(3) == 2
        numOfStimuli = time(2); % time is used to pass size(data) for legend
        colors = lines(numOfStimuli);
        hold on
        
        for f = 1:numOfStimuli
            plot([0.5 1], [f f], 'color', colors(f,:), 'LineWidth', 2)
            text(1.2, f, ['Stimulus ' num2str(f)])
        end
        plot([0.5 1], [-2 -2], 'k', 'LineWidth', 2)
        plot([0.5 1], [-3 -3], 'k--', 'LineWidth', 2)
        text(1.2, -2, 'Decision 1')
        text(1.2, -3, 'Decision 2')
        
        axis([0 3 -4.5 1.5+numOfStimuli])
        set(gca, 'XTick', [])
        set(gca, 'YTick', [])
        set(gca,'Visible','off')
        return
        
    elseif length(time) == 4 && time(3) == 4
        numOfStimuli = 4; % time is used to pass size(data) for legend
        colors = lines(numOfStimuli);
        stimLegends = ["control","Pattern 1","Pattern 2","Pattern 3"];
        hold on
        
        for f = 1:numOfStimuli
            plot([0.5 1], [f f], 'color', colors(f,:), 'LineWidth', 2)
            text(1.2, f, stimLegends(f));
        end
        plot([0.5 1], [-1 -1], 'k', 'LineWidth', 2)
        plot([0.5 1], [-2 -2], 'k--', 'LineWidth', 2)
        plot([0.5 1], [-3 -3], 'k:', 'LineWidth', 2)
        plot([0.5 1], [-4 -4], 'k-.', 'LineWidth', 2)
        
        text(1.2, -1, 'Angle 45')
        text(1.2, -2, 'Angle 135')
        text(1.2, -3, 'Angle 225')
        text(1.2, -4, 'Angle 315')
        
        axis([0 3 -4.5 1.5+numOfStimuli])
        set(gca, 'XTick', [])
        set(gca, 'YTick', [])
        set(gca,'Visible','off')
        return
    % other cases - do nothing
    else
        return
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% setting up the subplot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(time)
    time = 1:size(data, ndims(data));
end
axis([time(1) time(end) yspan])
hold on

if ~isempty(explVar)
    title(['Component #' num2str(compNum)])
     fprintf(['explained var ' num2str(explVar,'%.1f') '\n'])
else
    title(['Component #' num2str(compNum)])
end

if ~isempty(events)
%     plot([events; events], yspan, 'Color', [0.6 0.6 0.6])
end

if ~isempty(signif)
    signif(signif==0) = nan;
    plot(time, signif + yspan(1) + (yspan(2)-yspan(1))*0.05, 'k', 'LineWidth', 3)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotting the component
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ndims(data) == 2
    % only time - plot it
    plot(time, squeeze(data(1, :)), 'k', 'LineWidth', 2)

elseif ndims(data) == 3
    % different stimuli in different colours
    numOfStimuli = size(data, 2);
    plot(time, squeeze(data(1,:,:)), 'LineWidth', 2)    

elseif ndims(data) == 4 && size(data,3)==4
% different stimuli in different colours and binary condition as
    % solid/dashed
    numOfStimuli = size(data, 2);
    % matlab default colors
    gry = [0.5 0.5 0.5];
    yel = [0.9290, 0.6940, 0.1250];
    pur = [0.4940, 0.1840, 0.5560];
    olv = [0.4660, 0.6740, 0.1880];
    colors = [gry;olv;pur;yel]; 
    
    % angle color
    angcol  = [[0/255 0/255 205/255];[255/255 174/255 185/255];
        [205/255 38/255 38/255];[0/255 191/255 255/255]];
    
    if whichplot==1; colors = angcol;end
    if whichplot==2; colors = repmat([0 0 0],4,1);end

    for f=1:numOfStimuli 
        if whichplot == 2 || whichplot == 3
        plot(time, squeeze(data(1, f, 1, :)), 'color', colors(f,:), 'LineWidth', 2)
        plot(time, squeeze(data(1, f, 2, :)), 'color', colors(f,:), 'LineWidth', 2)
        plot(time, squeeze(data(1, f, 3, :)), 'color', colors(f,:), 'LineWidth', 2)
        plot(time, squeeze(data(1, f, 4, :)), 'color', colors(f,:), 'LineWidth', 2)
        end
        if whichplot == 1
        plot(time, squeeze(data(1, f, 1, :)), 'color', colors(1,:), 'LineWidth', 2)
        plot(time, squeeze(data(1, f, 2, :)), 'color', colors(2,:), 'LineWidth', 2)
        plot(time, squeeze(data(1, f, 3, :)), 'color', colors(3,:), 'LineWidth', 2)
        plot(time, squeeze(data(1, f, 4, :)), 'color', colors(4,:), 'LineWidth', 2)        
        end
        
    end

else
    % in all other cases pool all conditions and plot them in different
    % colours
    data = squeeze(data);
    dims = size(data);
    data = permute(data, [numel(dims) 1:numel(dims)-1]);
    data = reshape(data, size(data,1), []);
    data = data';
    
    plot(time, data, 'LineWidth', 2)    
end
