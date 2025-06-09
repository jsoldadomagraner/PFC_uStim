
function [datapath, statspath, figurepath] =  addpaths

% Place where data is stored, and where stats and figures will be stored
% The data is publicly available on KiltHub (see README)
mypath     = '~/PFC_uStim_data/';

datapath   = [mypath 'data/'];
statspath  = [mypath 'stats/'];
figurepath = [mypath 'figures/'];

% path to this code package
codepath = '~/PFC_uStim/';
addpath(genpath(codepath))

end