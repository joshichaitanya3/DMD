# dmd
MATLAB Function to perform the dynamic mode decomposition (DMD) on spatio-temporal data spaced evenly in time.

Adapted from codes available freely at http://dmdbook.com/

(see https://tinyurl.com/sksodx5 for more details about the method.)

Usage:

This example shows how to compute DMD on 2D velocity field data.
Assume you have velocities ux and uy stored on grids of size
(NX,NY) at times 1 through m+1, equally spaced in time. The data
matrix would be of size (n,m+1) with n=2*n0, and n0=NX*NY, with
data(1:n0,k) being ux at time t_k, flattened as a vector, and
similarly data(n0+1:2*n0,k) the uy. Then, the DMD can be computed
using:

X = data(:,1:end-1) % Size (n,m)   
Xprime = data(:,2:end) % Size (n,m)

% If all the optimal thresholded modes are desired

[U,Sigma,V,Phi,Lambda,b,r_optimal] = DMD(X,Xprime)

% If only r modes are desired

r = 21; % Set desired value

[U,Sigma,V,Phi,Lambda,b,r_optimal] = DMD(X,Xprime,r)
