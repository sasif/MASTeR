% clear;
% disp(['-----------------',datestr(now),'-----------------']);
% 
% %% DATA selection
% vid_seq_name = 'foreman';
% opts = [];
% opts.format = 'qcif';
% opts.frames = 10;
% opts.square = 1;
% Y_channel = read_Yseq(vid_seq_name,opts);
% T_frames = opts.frames;
% 
% I_cube = [];
% for frame = 1:T_frames
%     I_cube(:,:,frame) = Y_channel(:,:,frame)';
% end
% [ROW COL T_frames] = size(I_cube);
% 
% Ir_cube_me = I_cube;

vec = @(z) z(:);

%% ARPS + OBMC as in kt-FOCUSS
% ME/MC
% mbSize=2; %block size
% p=7; %searching window size/2

%% MASK and index pre-compilation
% For MC we start off from the top left of the image
% we will walk in steps of 1
% for every marcoblock that we look at we will read the motion vector
% and put that macroblock from refernce image in the compensated image
%
% Here we precompute the MASK required to normalize the overlapped MC
% blocks pixels
MASK=zeros(ROW,COL);
hor_OB = zeros((ROW-mbSize+1)*(COL-mbSize+1)*mbSize^2,1);
ver_OB = zeros((ROW-mbSize+1)*(COL-mbSize+1)*mbSize^2,1);
mbCount = 1;
t_ind = 0;
for i = 1:ROW-mbSize+1
    for j = 1:COL-mbSize+1
        hor_OB(t_ind+1:t_ind+mbSize^2) = vec(repmat(j:j+mbSize-1,mbSize,1));
        ver_OB(t_ind+1:t_ind+mbSize^2) = vec(repmat(i:i+mbSize-1,1,mbSize));
        t_ind = t_ind+mbSize^2;
        
        MASK(i:i+mbSize-1,j:j+mbSize-1) = MASK(i:i+mbSize-1,j:j+mbSize-1)+1;
        mbCount = mbCount + 1;
    end
end
MASK(MASK==0)=1;

%% ME/MC
% fprintf('Backward      ');
%% Backward ME
BI_cube = zeros(ROW, COL,T_frames);
MOTION_FIELD_BACKWARD = [];
MOTION_FIELD_BACKWARD = zeros(2,length(hor_OB),T_frames);
for frame=1:T_frames
    % fprintf('\b\b\b\b\b%0.2d...',frame);
    
    % Motion Estimation via Block Matching with sliding window
    % [motionVect, EScomputations, imgC]  = motionEstES_full(imgP, imgI, mbSize, p);
    % [motionVect, EScomputations, imgC]  = motionEst4SS_full(imgP, imgI, mbSize, p);
    % [motionVect, EScomputations, imgC]  = motionEstARPS_full(imgP, imgI, mbSize, p);
    % motionVect(1,:) - vertical direction (row co-ordinates)
    % motionVect(2,:) - horizontal direction (col co-ordinates)
    
    
    if frame == T_frames
        imgI = Ir_cube_me(:,:,1);
    else
        imgI = Ir_cube_me(:,:,frame+1);
    end
    imgP = Ir_cube_me(:,:,frame);
        
    motionVect = motionEst(imgP, imgI, mbSize, p);
    % MEMCb_ref(:,:,frame)=motionComp_full(imgI, motionVect, mbSize);
    
    hor_MC = hor_OB+(vec(repmat(motionVect(2,:),mbSize^2,1)));
    ver_MC = ver_OB+(vec(repmat(motionVect(1,:),mbSize^2,1)));
    if frame == T_frames
        MOTION_FIELD_BACKWARD(:,:,1) = [hor_MC'; ver_MC'];
    else
        MOTION_FIELD_BACKWARD(:,:,frame+1) = [hor_MC'; ver_MC'];
    end
    BI_cube(:,:,frame) = operator_OBMC(imgI, hor_OB, ver_OB, hor_MC, ver_MC, MASK);
end
% BI_cube(:,:,T_frames) = Ir_cube_me(:,:,T_frames);

% fprintf('\n');

%% Forward ME
% fprintf('Forward      ');
FI_cube = zeros(ROW,COL,T_frames);
FI_cube(:,:,1) = Ir_cube_me(:,:,1);
MOTION_FIELD_FORWARD = [];
MOTION_FIELD_FORWARD = zeros(2,length(hor_OB),T_frames);
for frame=1:T_frames
    % fprintf('\b\b\b\b\b%0.2d...',frame);
    
    
    imgI = Ir_cube_me(:,:,frame);
    if frame == T_frames
        imgP = Ir_cube_me(:,:,1);
    else
        imgP = Ir_cube_me(:,:,frame+1);
    end
    motionVect = motionEst(imgP, imgI, mbSize, p);
    % MEMCf_ref(:,:,frame)=motionComp_full(imgI, motionVect, mbSize);
    
    hor_MC = hor_OB+vec(repmat(motionVect(2,:),mbSize^2,1));
    ver_MC = ver_OB+vec(repmat(motionVect(1,:),mbSize^2,1));
    
    MOTION_FIELD_FORWARD(:,:,frame) = [hor_MC'; ver_MC'];
    if frame == T_frames
        FI_cube(:,:,1) = operator_OBMC(imgI, hor_OB, ver_OB, hor_MC, ver_MC, MASK);
    else
        FI_cube(:,:,frame+1) = operator_OBMC(imgI, hor_OB, ver_OB, hor_MC, ver_MC, MASK);
    end
end
% fprintf('\n');
