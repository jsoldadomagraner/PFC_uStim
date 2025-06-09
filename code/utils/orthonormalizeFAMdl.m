% Given an FA model and (optional) infered posterior means for latents,
% this function will orthonormalize the FA model and transform the means
% into the coordinate system for the new model. 
%
% Usage: [faMdl, postMeans] = orthonormalizeFAMdl(faMdl, postMeans)
%
% Inputs: 
%
%   faMdl - A structure with the following fields:
%
%       d - A P by 1 vector of the model mean, where P is the number of
%       observed variables. 
%
%       c - A P by L matrix giving the estimated loading matrix, where L is the
%       number of latent variables. 
%
%       psi - A P by P diagonal matrix of the estimated private noise
%       variances. 
%
%   postMeans - A N by L matrix of the posterior mean of the latent
%   variables, where N is the number of samples. 
%
% Optional Inputs: All optional inputs should be entered in string-value
% pair format. 
%
%   NO_ASSERT - True if assertions should be skipped.  Good for efficiency
%   buy may not catch unexpected inputs.  Default: false. 
%
% Outputs: 
%
%   faMdl - the faMdl with the field orthC added, which contains the
%   orthonormalized C matrix. 
%
%   postMeans - the posterior means expressed in the coordinate system for
%   the new cOrth. 
%
% Author: 
%
%   wbishop@cs.cmu.edu
%
function [faMdl, postMeans] = orthonormalizeFAMdl(faMdl, postMeans, varargin)

NO_ASSERT = false;
warnOpts(assignOpts(varargin)); 

nLatents = size(faMdl.c,2); 
nVars = size(faMdl.c,1); 
minVl = min(nLatents, nVars); 

% Calculate new C matrix 
[U, Sigma, V] = svd(faMdl.c); 
% The SVD almost gives a unique decomposition of C, but we still need to
% account for the the sign on the columns of U.  
for lI = 1:minVl
    uColSign = sign(sum(U(:,lI))); 
   if uColSign < 0
        U(:,lI) = -1*U(:,lI); 
        Sigma(lI,lI) = -1*Sigma(lI,lI);
    end
end

faMdl.orthC = U(:,1:minVl);

% Calculate new posterior means
if nargin > 1
    if nVars >= nLatents
        Sigma = Sigma(1:nLatents, 1:nLatents);
        postMeans = (Sigma*V'*postMeans')';
    else
        postMeans = (Sigma*V'*postMeans')'; 
    end
end

