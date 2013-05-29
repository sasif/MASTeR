MASTeR
======

Motion-adapative spatio-temporal regularization for accelerated dynamic MRI


%-------------------------------------------+
% MATLAB code for Motion-Adaptive Spatio-Temporal Regularization (MASTeR)
% for Accelerated Dynamic MRI
%-------------------------------------------+
% 
% Author: M. Salman Asif @ Georgia Tech
% Email: sasif@gatech.edu
% Web: http://users.ece.gatech.edu/~sasif/dynamicMRI
% Created: October 2011 -- March 2012
% Released: October 2012
% 
%-------------------------------------------+
% Copyright (c) 2012.  M. Salman Asif 
%-------------------------------------------+

%% Start
% *** add all the subfolders in MATLAB path and compile mex files where necessary *** 

% Run a single experiment using 
% demo_MASTeR 
% 
% Run a simulation with different values of reduction factor using
% job_MASTeR 
% job_ktMEMC
% job_rwt_MASTeR

%% Reading original data
%
% k-space data and sensitivity maps are organized in data folder as 
% 
% DATA USED IN THE PAPER
%
% Short-axis (with prescan)
% P501767_* or  P49664.7_  
% (To simulate prescan behavior, use one of these two scans as a prescan,
% which means use its sensitivity maps to reconstruct the other)
%
% wrap_003_*    : Two-chamber
% (Without prescan; sensitivity maps were generated using a center portion 
% of the original k-space)

%% Paper: 
% 
% "Motion-Adaptive Spatio-Temporal Regularization (MASTeR) for Accelerated Dynamic MRI"
% by Salman Asif, Lei Hamilton, Marijn Brummer, and Justin Romberg
% Accepted in Magnetic Resonance in Medicine, August 2012
% Preprint available at http://users.ece.gatech.edu/~sasif/

%% Acknowledgements:
%
% k-space and sensitivity maps data provided by L. Hamilton and M. Brummer
% along with scripts to read those files. 
%
% NESTA package for L1-norm minimization 
% available at http://www-stat.stanford.edu/~candes/nesta/
%
% Optical flow code available at 
% http://people.csail.mit.edu/celiu/OpticalFlow/

%%
%-------------------------------------------+
% Copyright (2012): M. Salman Asif
%-------------------------------------------+
%
% MASTeR is distributed under the terms of the GNU General Public
% License 2.0.
% 
% Permission to use, copy, modify, and distribute this software for
% any purpose without fee is hereby granted, provided that this entire
% notice is included in all copies of any software which is or includes
% a copy or modification of this software and in all copies of the
% supporting documentation for such software.
% This software is being provided "as is", without any express or
% implied warranty.  In particular, the authors do not make any
% representation or warranty of any kind concerning the merchantability
% of this software or its fitness for any particular purpose."
%
%-------------------------------------------+


%% RELEASE INFORMATION:
%
% First release: October 2012
%
% This code is in development stage; thus any comments or 
% bug reports are very welcome.
% 
% Please feel free to send questions or comments to
% 
% Salman Asif -- sasif@gatech.edu
% Justin Romberg -- jrom@ece.gatech.edu
%

