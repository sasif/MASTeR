J = floor(log2(min(ROW,COL)))-3;
% Complex wavelets
[Faf, Fsf] = FSfarras; % 1st stage anal. & synth. filters
% [Faf, Fsf] = AntonB;
[af, sf] = dualfilt1;

% symF = 1; symh = 1; symg = 2;
% [Faf, Fsf, af, sf] = BiOrthDualFilt_mod;

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

%%
MOTION_FIELD_FORWARD = [];
for frame = 1:T_frames
    Ip_m = Ir_cube_me(:,:,frame); % Reference image
    Ip_m = Ir_cube_me(AR_ver,AR_hor,frame);
    [row col] = size(Ip_m);
    if frame == T_frames
        Ie = Ir_cube_me(:,:,1); % assuming periodicity of image seq.
        Ie = Ir_cube_me(AR_ver,AR_hor,1); % assuming periodicity of image seq.
    else
        Ie = Ir_cube_me(:,:,frame+1); % Target image
        Ie = Ir_cube_me(AR_ver,AR_hor,frame+1); % assuming periodicity of image seq.
    end
    Ip_m = abs(Ip_m);
    Ie = abs(Ie);
    
    % Motion compensation with CWT
    switch motion_type
        case 1
            Omega_nm = estimate_localFreq(Ip_m, 1, col, row, J, Faf, af); % this represents 2^scale \Omega^{n,m)
            opts.Omega_nm = Omega_nm;
            
            wp = cplxdual2D_mod(Ip_m, J, Faf, af);
            wk = cplxdual2D_mod(Ie, J, Faf, af);
            
            %             % wp = cplxdual2D_mod(Ip_m, J_m, Faf, af);
            %             Fsep =  CWT2_separate(Ip_m, Faf, af, J, symF, symh, symg);
            %             W_temp = Fsep;% [Whh Whg Wgh Wgg]; % without pm
            %             wsym_temp = CWTmat2struct(W_temp,J);
            %             wp = separate2cw(wsym_temp,J_m);
            %
            %             % wk = cplxdual2D_mod(Ie, J_m, Faf, af);
            %             Fsep =  CWT2_separate(Ie, Faf, af, J, symF, symh, symg);
            %             W_temp = Fsep;
            %             wsym_temp = CWTmat2struct(W_temp,J);
            %             wk = separate2cw(wsym_temp,J);
            
            [f01 f02] = motionEstCWT_mod(wp, wk, opts,  J, Faf, Fsf, af, sf, row, col);
            
            opts.itype = 'linear';
            [hor_ind1 ver_ind1] = GetWarpedLocs(f01,f02,opts, J, row, col);
            [hor_ind ver_ind] = meshgrid(1:COL,1:ROW);
            hor_ind(AR_ver,AR_hor) = hor_ind1+AR_hor(1)-1;
            ver_ind(AR_ver,AR_hor) = ver_ind1+AR_ver(1)-1;
            
            % Icw = linear_forward(Ip_m, hor_ind, ver_ind);
            %     Icw = round(Icw);
            %     Icw_neg = Icw<0;
            %     Icw(Icw_neg) = 0;
            
        case 2
            % Block motion estimation.
            blk_size = 4;
            p = 4;
            [mv, es, Ib] = motionEstES_mod(Ie,Ip_m, blk_size, p);
            f01_b = reshape(mv(2,:),COL/blk_size, ROW/blk_size)';
            f02_b = reshape(mv(1,:),COL/blk_size, ROW/blk_size)';
            opts.itype = 'box';
            [hor_ind ver_ind] = GetWarpedLocs(f01_b,f02_b,opts,[],ROW,COL);
            
        case 3
            %% Control grid interpolation
            mask = ones(ROW,COL);
            [col_indices,row_indices] = meshgrid(1:COL, 1:ROW);
            [mvrow, mvcol, eout] = get_motion_field(Ie,Ip_m,8,1,1,mask);
            [xx yy] = meshgrid(1:COL,1:ROW);
            int_locs = xx+mvcol;
            int_locs_neg = int_locs<1;
            int_locs_out = int_locs>COL;
            int_locs(int_locs_neg) = 1;
            int_locs(int_locs_out) = COL;
            hor_ind = int_locs;
            
            int_locs = yy+mvrow;
            int_locs_neg = int_locs<1;
            int_locs_out = int_locs>ROW;
            int_locs(int_locs_neg) = 1;
            int_locs(int_locs_out) = ROW;
            ver_ind = int_locs;
        otherwise
            disp('chan oye chal');
    end
    
    MOTION_FIELD_FORWARD{frame}{1} = hor_ind;
    MOTION_FIELD_FORWARD{frame}{2} = ver_ind;
end

%%
MOTION_FIELD_BACKWARD = [];
for frame = T_frames:-1:1
    Ip_m = Ir_cube_me(:,:,frame); % REFERENCE FRAME
    Ip_m = Ir_cube_me(AR_ver,AR_hor,frame);
    
    if frame == 1
        Ie = Ir_cube_me(:,:,T_frames); % assuming periodicity of image seq.
        Ie = Ir_cube_me(AR_ver,AR_hor,T_frames); % assuming periodicity of image seq.
    else
        Ie = Ir_cube_me(:,:,frame-1); % Target image
        Ie = Ir_cube_me(AR_ver,AR_hor,frame-1); % assuming periodicity of image seq.
    end
    Ip_m = abs(Ip_m);
    Ie = abs(Ie);
    
    % Motion compensation with CWT
    switch motion_type
        case 1
            Omega_nm = estimate_localFreq(Ip_m, 1, col, row, J, Faf, af); % this represents 2^scale \Omega^{n,m)
            opts.Omega_nm = Omega_nm;
            
            wp = cplxdual2D_mod(Ip_m, J, Faf, af);
            wk = cplxdual2D_mod(Ie, J, Faf, af);
            
            %             % wp = cplxdual2D_mod(Ip_m, J_m, Faf, af);
            %             Fsep =  CWT2_separate(Ip_m, Faf, af, J, symF, symh, symg);
            %             W_temp = Fsep;% [Whh Whg Wgh Wgg]; % without pm
            %             wsym_temp = CWTmat2struct(W_temp,J);
            %             wp = separate2cw(wsym_temp,J_m);
            %             % W_sym = CWTstruct2mat(wp);
            %             % figure(1001); imagesc([W_sym- W_tot]);
            %             % max(abs(W_sym(:)-W_tot(:)))
            %
            %             % wk = cplxdual2D_mod(Ie, J_m, Faf, af);
            %             Fsep =  CWT2_separate(Ie, Faf, af, J, symF, symh, symg);
            %             W_temp = Fsep;
            %             wsym_temp = CWTmat2struct(W_temp,J);
            %             wk = separate2cw(wsym_temp,J);
            
            [f01 f02] = motionEstCWT_mod(wp, wk, opts,  J, Faf, Fsf, af, sf, row, col);
            
            opts.itype = 'linear';
            [hor_ind1 ver_ind1] = GetWarpedLocs(f01,f02,opts, J, row, col);
            [hor_ind ver_ind] = meshgrid(1:COL,1:ROW);
            hor_ind(AR_ver,AR_hor) = hor_ind1+AR_hor(1)-1;
            ver_ind(AR_ver,AR_hor) = ver_ind1+AR_ver(1)-1;
            
            % Icw = linear_forward(Ip_m, hor_ind, ver_ind);
            %     Icw = round(Icw);
            %     Icw_neg = Icw<0;
            %     Icw(Icw_neg) = 0;
            
        case 2
            % Block motion estimation.
            blk_size = 4;
            p = 4;
            [mv, es, Ib] = motionEstES_mod(Ie,Ip_m, blk_size, p);
            f01_b = reshape(mv(2,:),COL/blk_size, ROW/blk_size)';
            f02_b = reshape(mv(1,:),COL/blk_size, ROW/blk_size)';
            opts.itype = 'box';
            [hor_ind ver_ind] = GetWarpedLocs(f01_b,f02_b,opts,[],ROW,COL);
            
        case 3
            %% Control grid interpolation
            mask = ones(ROW,COL);
            [col_indices,row_indices] = meshgrid(1:COL, 1:ROW);
            [mvrow, mvcol, eout] = get_motion_field(Ie,Ip_m,8,1,1,mask);
            [xx yy] = meshgrid(1:COL,1:ROW);
            int_locs = xx+mvcol;
            int_locs_neg = int_locs<1;
            int_locs_out = int_locs>COL;
            int_locs(int_locs_neg) = 1;
            int_locs(int_locs_out) = COL;
            hor_ind = int_locs;
            
            int_locs = yy+mvrow;
            int_locs_neg = int_locs<1;
            int_locs_out = int_locs>ROW;
            int_locs(int_locs_neg) = 1;
            int_locs(int_locs_out) = ROW;
            ver_ind = int_locs;
        otherwise
            disp('chan oye chal');
    end
    
    MOTION_FIELD_BACKWARD{frame}{1} = hor_ind;
    MOTION_FIELD_BACKWARD{frame}{2} = ver_ind;
end


%% Construct MC frames
FI_cube = zeros(ROW,COL,T_frames);
BI_cube = zeros(ROW,COL,T_frames);
for frame = 1:T_frames
    %     switch frame
    %         case 1
    %             In = Ir_cube_me(:,:,frame+1);
    %             Ip = Ir_cube_me(:,:,T_frames);
    %         case T_frames
    %             In = Ir_cube_me(:,:,1);
    %             Ip = Ir_cube_me(:,:,frame-1);
    %         otherwise
    %             In = Ir_cube_me(:,:,frame+1);
    %             Ip = Ir_cube_me(:,:,frame-1);
    %     end
    %     ind_p = frame; ind_n = frame;
    
    ind_n = mod(frame,T_frames)+1;
    ind_p = mod(frame-2,T_frames)+1;
    
    In = Ir_cube_me(:,:,ind_n);
    Ip = Ir_cube_me(:,:,ind_p);
    
    hor_ind = MOTION_FIELD_FORWARD{ind_p}{1};
    ver_ind = MOTION_FIELD_FORWARD{ind_p}{2};
    FI_cube(:,:,frame) = linear_mc(Ip,hor_ind,ver_ind,ROW,COL);
    
    hor_ind = MOTION_FIELD_BACKWARD{ind_n}{1};
    ver_ind = MOTION_FIELD_BACKWARD{ind_n}{2};
    BI_cube(:,:,frame) = linear_mc(In,hor_ind,ver_ind,ROW,COL);
end