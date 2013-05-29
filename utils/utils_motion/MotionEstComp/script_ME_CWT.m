% function [hor_ind, ver_ind] = MotionEst_CWT_fbf(Ip_m, Ie)
%
% Script for computing motion compensated residual between two
% images with preset parameters 
% Ip_m : reference frame
% Ie : target frame
%

J_m = floor(log2(min(ROW,COL)))-3;

% Complex wavelets
% [Faf, Fsf] = FSfarras; % 1st stage anal. & synth. filters
% % [Faf, Fsf] = AntonB;
% [af, sf] = dualfilt1;

symF = 1; symh = 1; symg = 2;
[Faf, Fsf, af, sf] = BiOrthDualFilt_mod;

opts = [];
opts.Ic_type = 1; % 1 - full interpolation, 2 - block replacement
opts.iker = 1; % 1 - box interpolation, 2 - linear interpolation for wavelet coefficients
opts.plots = 0;
opts.SD_combine = 1;
opts.wt_SD = 4; % 4 has better psnr, 2 looks better in some cases
opts.motion_wrap = 1;
opts.stp_scale = 1;
opts.itype = 'linear';
opts.itype_dir = 0;

% Motion compensation with CWT
Omega_nm = estimate_localFreq(Ip_m, 1, COL, ROW, J_m, Faf, af); % this represents 2^scale \Omega^{n,m)
opts.Omega_nm = Omega_nm;

%     wp = cplxdual2D_mod(Ip_m, J_m, Faf, af);
%     wk = cplxdual2D_mod(Ie, J_m, Faf, af);

Fsep =  CWT2_separate(Ip_m, Faf, af, J_m, symF, symh, symg);
W_temp = Fsep;% [Whh Whg Wgh Wgg]; % without pm
wsym_temp = CWTmat2struct(W_temp, J_m);
wp = separate2cw(wsym_temp,J_m);
% W_sym = CWTstruct2mat(wp);
% W_tot = CWTstruct2mat(wp);
% figure(1001); imagesc([W_sym- W_tot]);
% max(abs(W_sym(:)-W_tot(:)))
%  max(max(abs(W_tot-W_sym)))

Fsep =  CWT2_separate(Ie, Faf, af, J_m, symF, symh, symg);
W_temp = Fsep;
wsym_temp = CWTmat2struct(W_temp, J_m);
wk = separate2cw(wsym_temp,J_m);

[f01 f02] = motionEstCWT_mod(wp, wk, opts,  J_m, Faf, Fsf, af, sf, ROW, COL);

opts.itype = 'linear';
[hor_ind ver_ind] = GetWarpedLocs(f01,f02,opts, J_m, ROW, COL);
