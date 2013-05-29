% Setup paths

mname = mfilename;
mpath = mfilename('fullpath');
mdir = mpath(1:end-length(mname));
cd(mdir);

addpath(genpath('data'));
addpath(genpath('utils'));
addpath('operators');


%% May need to compile the following mex files.
% Wavelet transform (fwt, fwt2, fwt2_CWT, ifwt, ifwt2, ifwt2_CWT...)

if exist('fwt')~=3
    cd utils\utils_wavelet
    disp('Compiling CWT');
    Compile_CWT
    cd(mdir);
end
