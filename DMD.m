function [U,Sigma,V,Phi,Lambda,b,r_optimal] = DMD(X,Xprime,varargin)

% function [U,Sigma,V,Phi,Lambda,b,r_optimal] = DMD(X,Xprime,varargin)

% Function to perform the dynamic mode decomposition (DMD) on
% spatio-temporal data. 
% Adapted from codes available freely at http://dmdbook.com/
% (see https://tinyurl.com/sksodx5 for more details about the method.)

% IN: 

%      X: Data matrix of dimensions (n,m), where m is the number of
%         snapshots in time and n is the number of pixels / values 
%         per snapshot. For example, if you have a scalar vorticity 
%         field on a 2D grid of size (NX,NY), then n=NX*NY and each
%         column of X would be the flattened array of the vorticity 
%         field at that time instance. For multiple fields, X would 
%         have columns that are the concatenated fields. Note that 
%         the time interval between any two snapshots need not be 
%         the same throughout X.
%
% Xprime: The same data matrix as X, but shifted by a *fixed* time 
%         interval dt. So, in a scenario where you fixed sampling 
%         time in the data with m+1 samples in total, X would have
%         samples 1 through m and Xprime would have samples 
%         2 through m+1. 
%      
%      r: (optional) Number of DMD modes to be calculated.
%         If this is not provided, then the number of modes to be
%         kept will be calculated using the Optimal Hard Threshold
%         for SVD modes
%         (see the code `optimal_SVHT_coef.m' for more details)
% 
% OUT:

%    [U,Sigma,V]: The output of svd(X) 

% [Phi,Lambda,b]: The DMD modes, eigenvalues and amplitudes
%                 (see https://tinyurl.com/sksodx5)
%
%      r_optimal: The number of SVD modes kept using the Optimal Hard
%                 Threshold. (see the code `optimal_SVHT_coef.m'
%                 for more details)         
%
% Usage:
%
%    This example shows how to compute DMD on 2D velocity field data.
%    Assume you have velocities ux and uy stored on grids of size
%    (NX,NY) at times 1 through m+1, equally spaced in time. The data
%    matrix would be of size (n,m+1) with n=2*n0, and n0=NX*NY, with
%    data(1:n0,k) being ux at time t_k, flattened as a vector, and
%    similarly data(n0+1:2*n0,k) the uy. Then, the DMD can be computed
%    using:
%   
%   X = data(:,1:end-1) % Size (n,m)   
%   Xprime = data(:,2:end) % Size (n,m)
%
%   % If all the optimal thresholded modes are desired
%   [U,Sigma,V,Phi,Lambda,b,r_optimal] = DMD(X,Xprime)
%
%   % If only r modes are desired

%   r = 21; % Set desired value
%   [U,Sigma,V,Phi,Lambda,b,r_optimal] = DMD(X,Xprime,r)


    % evaluating the SVD
    [n,m] = size(X);
    [U,Sigma,V] = svd(X,'econ');

    % Compute the hard threshold cutoff for SVD:
    % reference for code: https://purl.stanford.edu/vg705qn9070
    % reference for paper: https://arxiv.org/abs/1305.5870

    % Note that the code online has a typo, which is corrected in the
    % file in this directory.
    sigma = diag(Sigma);
    opt_hard_thresh = optimal_SVHT_coef(m/n,0) * median(sigma);
    r_optimal = size(nonzeros(sigma > opt_hard_thresh),1);
    % r_optimal = size(nonzeros(sigma > (optimal_SVHT_coef(m/n,0) *
    %                                   median(sigma))),1);
    disp(['Total modes = ' num2str(m)]);
    disp(['Optimal # of modes = ' num2str(r_optimal)]);

    if nargin==2
        disp('Number of modes not provided. Using the optimal value.');
        r = r_optimal;
    elseif nargin==3
        arg = varargin{1};
        if ~isa(arg,'double')
            error(['ValueError: Expected an integer as the' ...
                            ' optional third argument, '...
                            'but got a ' class(arg) '.']);
        elseif arg>m
            % If desired number of modes is greater than the total
            % number of modes, compute all the modes.
            r = m;
            disp(['Desired number of modes is greater than the ' ...
                'total number of modes. Computing all the modes...'])
        else
            r = int32(arg);
            disp(['Computing only ' num2str(r) ...
                ' modes as requested...'])
        end
    elseif nargin>3
        error('ValueError: Expected only 2 or 3 arguments.');
    end
    
    Ur = U(:,1:r);
    Sigmar = Sigma(1:r,1:r);
    Vr = V(:,1:r); % Step 1
    Atilde = Ur'*Xprime*Vr/Sigmar;
    [W,Lambda] = eig(Atilde); % Step 2
    % Step 3
    Phi = Xprime*(Vr/Sigmar)*W;
    alpha1 = Sigmar*Vr(1,:)';
    b = (W*Lambda)\alpha1; % Step 4
end