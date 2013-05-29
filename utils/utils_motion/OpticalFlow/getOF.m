% method = 'Sun';% 'Sun' or 'Liu'
% opts.method = 'classic++';
% % Method can be
% %   'classic+nl-fast' (default)  'classic+nl' 'classic+nl-full'
% %   'classic++'  'classic-c'  'classic-l'/'ba' 'hs'
% opts.para = [];


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
    
    Ip_m = abs(Ip_m); Ie = abs(Ie);
    
    [u,v] = handle_OF(Ie,Ip_m,method,opts);
    [xx yy] = meshgrid(1:COL,1:ROW);
    
    locs = xx+u;
    locs(locs<1) = 1; 
    locs(locs>COL) = COL;
    hor_ind = locs; 
    
    locs = yy+v;
    locs(locs<1) = 1; 
    locs(locs>ROW) = ROW;
    ver_ind = locs; 
    
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
    
    Ip_m = abs(Ip_m); Ie = abs(Ie);
    
    [u,v] = handle_OF(Ie,Ip_m,method,opts);
    [xx yy] = meshgrid(1:COL,1:ROW);
    
    locs = xx+u;
    locs(locs<1) = 1; 
    locs(locs>COL) = COL;
    hor_ind = locs; 
    
    locs = yy+v;
    locs(locs<1) = 1; 
    locs(locs>ROW) = ROW;
    ver_ind = locs; 
    
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

