
function [datapath, statspath, figurepath] =  addpaths

% Folder where data is stored, and where stats and figures will be stored
% This folder can be dowloaded from a public repository (see README)
mypath     = '~/PFC_uStim_data/';

datapath   = [mypath 'data/'];
statspath  = [mypath 'stats/'];
figurepath = [mypath 'figures/'];

% path to this code package
codepath = '~/PFC_uStim/';
addpath(genpath(codepath))

end
