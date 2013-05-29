%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Complex wavelet and motion estimation parameters %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%%%%%%%%%%%%%%%%%%%%
%% FORWARD motion %%
%%%%%%%%%%%%%%%%%%%%

MOTION_FIELD_FORWARD = [];

F_I0 = Ir_cube_me(:,:,1); 
FI_cube = zeros(ROW,COL,T_frames);

for frame = 1:T_frames
    Ip_m = Ir_cube_me(:,:,frame); % REFERENCE FRAME
    if frame == T_frames
        Ie = Ir_cube_me(:,:,1); % assuming periodicity of image seq.
    else
        Ie = Ir_cube_me(:,:,frame+1); % TARGET FRAME
    end
    Ip_m = Ip_m(AR_ver,AR_hor);
    [row col] = size(Ip_m);
    Ie = Ie(AR_ver,AR_hor); % assuming periodicity of image seq.
    
    Ip_m = abs(Ip_m);
    Ie = abs(Ie);
    
    % Motion compensation with CWT
    Omega_nm = estimate_localFreq(Ip_m, 1, col, row, J_m, Faf, af); % this represents 2^scale \Omega^{n,m)
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
    
    [f01 f02] = motionEstCWT_mod(wp, wk, opts,  J_m, Faf, Fsf, af, sf, row, col);
    
    opts.itype = 'linear';
    [hor_ind1 ver_ind1] = GetWarpedLocs(f01,f02,opts, J_m, row, col);
    [hor_ind ver_ind] = meshgrid(1:COL,1:ROW);
    hor_ind(AR_ver,AR_hor) = hor_ind1+AR_hor(1)-1;
    ver_ind(AR_ver,AR_hor) = ver_ind1+AR_ver(1)-1;
    
    % Icw = linear_forward(Ip_m, hor_ind, ver_ind);
    %     Icw = round(Icw);
    %     Icw_neg = Icw<0;
    %     Icw(Icw_neg) = 0;
    
    if rnd == 1
        MOTION_FIELD_FORWARD{frame}{1} = round(hor_ind);
        MOTION_FIELD_FORWARD{frame}{2} = round(ver_ind);
    else
        MOTION_FIELD_FORWARD{frame}{1} = hor_ind;
        MOTION_FIELD_FORWARD{frame}{2} = ver_ind;
    end
    Ip_m = Ir_cube_me(:,:,frame); % REFERENCE FRAME
    F_I0 = linear_mc(Ip_m,hor_ind, ver_ind, ROW, COL);
    if frame == T_frames
        FI_cube(:,:,1) = F_I0;
    else
        FI_cube(:,:,frame+1) = F_I0;
    end
end

%%%%%%%%%%%%%%%%%%%%
%% BACKWARD motion %%
%%%%%%%%%%%%%%%%%%%%

MOTION_FIELD_BACKWARD = [];
B_I1 = Ir_cube_me(:,:,T_frames);
BI_cube = zeros(ROW,COL,T_frames);
for frame = 1:T_frames
    Ip_m = Ir_cube_me(:,:,frame); % REFERENCE FRAME
    if frame == 1
        Ie = Ir_cube_me(:,:,T_frames); % assuming periodicity of image seq.
    else
        Ie = Ir_cube_me(:,:,frame-1); % TARGET FRAME
    end
    Ip_m = Ip_m(AR_ver,AR_hor);
    [row col] = size(Ip_m);
    Ie = Ie(AR_ver,AR_hor);
    
    Ip_m = abs(Ip_m);
    Ie = abs(Ie);
    
    Omega_nm = estimate_localFreq(Ip_m, 1, col, row, J_m, Faf, af); % this represents 2^scale \Omega^{n,m)
    opts.Omega_nm = Omega_nm;
    
    %      wp = cplxdual2D_mod(Ip_m, J_m, Faf, af);
    %      wk = cplxdual2D_mod(Ie, J_m, Faf, af);
    
    Fsep =  CWT2_separate(Ip_m, Faf, af, J_m, symF, symh, symg);
    W_temp = Fsep;% [Whh Whg Wgh Wgg]; % without pm
    wsym_temp = CWTmat2struct(W_temp,J_m);
    wp = separate2cw(wsym_temp,J_m);
    % W_sym = CWTstruct2mat(wp);
    % W_tot = CWTstruct2mat(wp);
    % figure(1001); imagesc([W_sym- W_tot]);
    % max(abs(W_sym(:)-W_tot(:)))
    %  max(max(abs(W_tot-W_sym)))
    
    Fsep =  CWT2_separate(Ie, Faf, af, J_m, symF, symh, symg);
    W_temp = Fsep;
    wsym_temp = CWTmat2struct(W_temp,J_m);
    wk = separate2cw(wsym_temp,J_m);
    
    [f01 f02] = motionEstCWT_mod(wp, wk, opts,  J_m, Faf, Fsf, af, sf, row, col);
    
    opts.itype = 'linear';
    [hor_ind1 ver_ind1] = GetWarpedLocs(f01,f02,opts, J_m, row, col);
    [hor_ind ver_ind] = meshgrid(1:COL,1:ROW);
    hor_ind(AR_ver,AR_hor) = hor_ind1+AR_hor(1)-1;
    ver_ind(AR_ver,AR_hor) = ver_ind1+AR_ver(1)-1;
    
    % Icw = linear_forward(Ip_m, hor_ind, ver_ind);
    %     Icw = round(Icw);
    %     Icw_neg = Icw<0;
    %     Icw(Icw_neg) = 0;
    
    if rnd == 1
        MOTION_FIELD_BACKWARD{frame}{1} = round(hor_ind);
        MOTION_FIELD_BACKWARD{frame}{2} = round(ver_ind);
    else
        MOTION_FIELD_BACKWARD{frame}{1} = hor_ind;
        MOTION_FIELD_BACKWARD{frame}{2} = ver_ind;
    end
    Ip_m = Ir_cube_me(:,:,frame); % REFERENCE FRAME
    B_I1 = linear_mc(Ip_m,hor_ind, ver_ind, ROW, COL);
    if frame == 1
        BI_cube(:,:,T_frames) = B_I1;
    else
        BI_cube(:,:,frame-1) = B_I1;
    end
end


%% Construct MC frames
% FI_cube = zeros(ROW,COL,T_frames);
% BI_cube = zeros(ROW,COL,T_frames);
% for frame = 1:T_frames
%     
%     ind_n = mod(frame,T_frames)+1;
%     ind_p = mod(frame-2,T_frames)+1;
%     
%     In = Ir_cube_me(:,:,ind_n);
%     Ip = Ir_cube_me(:,:,ind_p);
%     
%     hor_ind = MOTION_FIELD_FORWARD{ind_p}{1};
%     ver_ind = MOTION_FIELD_FORWARD{ind_p}{2};
%     FI_cube(:,:,frame) = linear_mc(Ip,hor_ind,ver_ind,ROW,COL);
%     
%     hor_ind = MOTION_FIELD_BACKWARD{ind_n}{1};
%     ver_ind = MOTION_FIELD_BACKWARD{ind_n}{2};
%     BI_cube(:,:,frame) = linear_mc(In,hor_ind,ver_ind,ROW,COL);
% end