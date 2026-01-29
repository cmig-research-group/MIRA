function [MIRA_vectors, muvec_lomot, V_mot] = MIRA_createComponents(corrmat, meanfdvec, nbins, refBin)
% [MIRA_vectors, muvec_lomot, V_mot] = MIRA_createComponents(corrmat, meanfdvec, nbins, refBin)
% Function to create MIRA projection vectors, given a set of concatenated
% resting state fMRI correlation matrices, and a vector of mean framewise
% displacement
%
%% Inputs:
% corrmat:      [n x v]     matrix of flattened correlation matrices for
%                           n subjects and v correlations
% 
% meanfdvec:    [n x 1]     vector of mean framewise displacement for n
%                           subjects (after censoring) 
% 
% nbins:        [1 x 1]     number indicating how many percentile bins
%                           should meanfdvec be divided into (default: 25)
% 
% refBin:       [1 x 1]     number indicating which percentile bin should
%                           be used as reference group to compute the bias
%                           from (default: low motion group is the 4th
%                           percentile, which is the first bin when nbins
%                           is 25)
% 
%% Output:
% MIRA_vectors: [n x p]     matrix of observations projected on the p
%                           projection vectors (where p is determined by
%                           the number of percentile bins that the
%                           meanfdvec is divided into; default to 25 bins)
% 
% muvec_lomot:  [1 x v]     vector of average corrmat for the low motion
%                           group
% 
% V_mot:        [v x b]     matrix of V matrix from the SVD decomposition,
%                           where b is the number of percentile bins that
%                           motion was split in
%
%% Notes:
% To project new data "corrmat_new" on to the same vectors, do:
%   bias_corrmat_new = corrmat_new - muvec_lomot;
%   MIRA_vectors_new = bias_corrmat_new * V_mot;

%% Check inputs
% Check corrmat
if ~exist('corrmat', 'var') || isempty(corrmat)
    error('Please provide a matrix containing the concatenated rsfMRI correlation matrix');
else
    nsubjs = size(corrmat, 1);
end

% Check meanfdvec
if ~exist('meanfdvec', 'var') || isempty(meanfdvec)
    error('Please provide a vector containing the mean FD values for each subject');
else
    % Ensure that this is a single column
    meanfdvec = reshape(meanfdvec, [], 1);

    if length(meanfdvec) ~= nsubjs
        error('Mismatch between number of subjects in the corrmat and in meanfdvec');
    end
end

% Check nbins
if ~exist('nbins', 'var') || isempty(nbins)
    nbins = 25;
else
    if ~isnumeric(nbins)
        error(['nbins should be a number specifying the number of percentile ', ...
               'bins meanfdvec should be divided into']);
    end
end

%% Step 1: divide mean FD into percentiles
allBins          = linspace(0, 100, nbins+1);
prctilelist      = prctile(meanfdvec, allBins);
prctilelist(end) = Inf;
nprcbins         = length(prctilelist)-1;

% Check refBin
if ~exist('refBin', 'var') || isempty(refBin)
    refBin = find(allBins,1);
else
    if ~isnumeric(refBin)
        error(['refBin should be a number specifying which percentile bin ', ...
               'should be used as reference']);
    else
        if refBin > length(allBins)
            error('refBin exceeds number of computed bins');
        end
    end
end

% Any columns of corrmat that needs to be ignored?
% defvec = isfinite(sum(corrmat,2));

%% Step 2: average of the low motion group
ivec_lomot  = meanfdvec < prctilelist(refBin); % & defvec
muvec_lomot = mean(corrmat(ivec_lomot,:), 'omitmissing');

%% Step 3: bias from the low motion average
corrmat_bias = corrmat - muvec_lomot;

%% Step 4: for every percentile bin, calculate average bias
mumat = nan(nprcbins, size(corrmat_bias,2));
for prci = 1:nprcbins
    ivec_tmp      = meanfdvec >= prctilelist(prci) & meanfdvec < prctilelist(prci+1); % & defvec
    mumat(prci,:) = mean(corrmat_bias(ivec_tmp,:), 'omitmissing');
end

%% Step 5: SVD on bias matrix for every percentile bin
[~, ~, V_mot] = svd(mumat, 'econ');

%% Step 6: project the bias matrix 
MIRA_vectors = corrmat_bias * V_mot;