function [col_ang, markers, anglab, stimcol] = plottingspecs

set(0,'DefaultAxesTitleFontWeight','normal');
set(0, 'DefaultAxesFontSize',18);

anglab  = {'45°','135°','225°','315°'};
col_ang = [[0/255 0/255 205/255];[255/255 174/255 185/255];
          [205/255 38/255 38/255];[0/255 191/255 255/255]];
      
markers  = {'square','diamond'};

% matlab default colors
gry = [0.5 0.5 0.5];
yel = [0.9290, 0.6940, 0.1250];
pur = [0.4940, 0.1840, 0.5560];
olv = [0.4660, 0.6740, 0.1880];

stimcol = [gry;yel;pur;olv];

end