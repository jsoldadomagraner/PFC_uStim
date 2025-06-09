function [analyp,uStimp,trialtimes,tpre,tpost] = analysisperiod(analyperiod)

% trial period. 0 marks uStim onset
binWidth   = 50; 
trialtimes = -800:binWidth:750;
% analysis period
analyp     = find(trialtimes>=analyperiod(1) & trialtimes<=analyperiod(2));
% uStim period (150 ms uStim + 50ms post-uStim to exclude artifacts)
uStimp     = find(trialtimes>=0 & trialtimes<200);
[~,uStimp,~] = intersect(analyp,uStimp);
% time steps before and after uStim
tpre  = uStimp(1)-1;
tpost = uStimp(end)+1;

end