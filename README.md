# Dynamic Mode Decomposition in MATLAB
MATLAB Function to perform the dynamic mode decomposition (DMD) on spatio-temporal data spaced evenly in time.

In simple terms, it decomposes the data into oscillating spatio-temporal patterns, with a fixed frequency and growth/decay rate.

## Sources
This script is based on the techniques and codes presented in the book 'Data-Driven Science and Engineering' by Steven L. Brunton and J. Nathan Kutz, as well as codes available on their [DMD book website](http://dmdbook.com/).

See Steve's video below for an excellent description of the method. 
<a href="http://www.youtube.com/watch?feature=player_embedded&v=sQvrK8AGCAo
" target="_blank"><img src="http://img.youtube.com/vi/sQvrK8AGCAo/0.jpg" 
alt="Link to Youtube video describing Dynamic Mode Decomposition" width="480" height="360" border="10" /></a>

The script for finding the optimal threshold for the modes is 
developed by D. L. Donoho and M. Gavish in ["The Optimal Hard Threshold for Singular
Values is 4/sqrt(3)"](http://arxiv.org/abs/1305.5870)

## Usage:

This example shows how to compute DMD on 2D velocity field data.
Assume you have velocities `ux`and `uy` stored on grids of size
`(NX,NY)` at times `1` through `m+1`, equally spaced in time. The data
matrix would be of size `(n,m+1)` with `n=2*n0`, and `n0=NX*NY`, with
`data(1:n0,k)` being `ux` at time `t_k`, flattened as a vector, and
similarly `data(n0+1:2*n0,k)` being `uy`. Then, the DMD can be computed
using:

```matlab
X = data(:,1:end-1) % Size (n,m)   
Xprime = data(:,2:end) % Size (n,m)

% If all the optimal thresholded modes are desired

[U,Sigma,V,Phi,Lambda,b,r_optimal] = DMD(X,Xprime)

% If only r modes are desired

r = 21; % Set desired value

[U,Sigma,V,Phi,Lambda,b,r_optimal] = DMD(X,Xprime,r)
```
