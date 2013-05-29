% Block motion estimation.
% ROW and COL must be divisible by blk_size

%%%%%%%%%%%%%%%%%%%%
%% FORWARD motion %%
%%%%%%%%%%%%%%%%%%%%

MOTION_FIELD_FORWARD = [];

F_I0 = Ir_cube_me(:,:,1); % REFERENCE FRAME
FI_cube = zeros(ROW,COL,T_frames);

for frame = 1:T_frames
    Ip_m = Ir_cube_me(:,:,frame); % Reference image
    if frame == T_frames
        Ie = Ir_cube_me(:,:,1); % assuming periodicity of image seq.
    else
        Ie = Ir_cube_me(:,:,frame+1); % Target image
    end
    
    Ip_m = abs(Ip_m);
    Ie = abs(Ie);
    
    [mv, es, Ib] = motionEstES_mod(Ie,Ip_m, blk_size, p);
    f01_b = reshape(mv(2,:),COL/blk_size, ROW/blk_size)';
    f02_b = reshape(mv(1,:),COL/blk_size, ROW/blk_size)';
    opts.itype = 'box';
    [hor_ind ver_ind] = GetWarpedLocs(f01_b,f02_b,opts,[],ROW,COL);
    
    % Icw = linear_forward(Ip_m, hor_ind, ver_ind);
    %     Icw = round(Icw);
    %     Icw_neg = Icw<0;
    %     Icw(Icw_neg) = 0;
    
    MOTION_FIELD_FORWARD{frame}{1} = hor_ind;
    MOTION_FIELD_FORWARD{frame}{2} = ver_ind;
    
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
        Ie = Ir_cube_me(:,:,frame-1); % Target image
    end
    
    Ip_m = abs(Ip_m);
    Ie = abs(Ie);
    
    [mv, es, Ib] = motionEstES_mod(Ie,Ip_m, blk_size, p);
    f01_b = reshape(mv(2,:),COL/blk_size, ROW/blk_size)';
    f02_b = reshape(mv(1,:),COL/blk_size, ROW/blk_size)';
    opts.itype = 'box';
    [hor_ind ver_ind] = GetWarpedLocs(f01_b,f02_b,opts,[],ROW,COL);
    
    % Icw = linear_forward(Ip_m, hor_ind, ver_ind);
    %     Icw = round(Icw);
    %     Icw_neg = Icw<0;
    %     Icw(Icw_neg) = 0;
    
    MOTION_FIELD_BACKWARD{frame}{1} = hor_ind;
    MOTION_FIELD_BACKWARD{frame}{2} = ver_ind;
    
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