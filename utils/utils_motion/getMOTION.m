% function [FI_cube BI_cube] = getMOTION(Ir_cube_me,MOTION_TYPE, AR_hor, AR_ver);
%
% Estimate the forward and backward motion and store motion vectors in 
% MOTION_FIELD_FORWARD and _BACKWARD
%
% Possible ME schemes to choose from
%   CWT     - complex wavelets phase based motion estimation
%   OBMC    - overlapped block matching 
%   OF      - optical flow schemes
%   BM      - block matching
% 
%
% global ROW COL T_frames MOTION_FIELD_FORWARD MOTION_FIELD_BACKWARD

switch MOTION_TYPE
    case 'CWT'
        motion_type = 1; rnd = 0;
        getMotionField_CWT
        % calculate_MotionField_bi;
    case 'OBMC'
        mbSize=2;   %block size
        p=7;        %searching window size/2
        getOBMC;    
    case 'OF'
        % Two possible schems I found the best are 
        method = 'Liu'; % set optical flow parameters (see Coarse2FineTwoFrames.m for the definition of the parameters)
        % method = 'Sun'; opts.method = 'classic+nl';
        rnd = 0;
        disp(sprintf('ME %s, method = %s, itr = %d',mtype,method,itr));
        alpha = 1;
        ratio = 0.5;
        minWidth = 40;
        nOuterFPIterations = 3;
        nInnerFPIterations = 1;
        nSORIterations = 20;
        para = [alpha,ratio,minWidth,nOuterFPIterations,nInnerFPIterations,nSORIterations];
        
        opts.para = para;
        getOF;
    case 'BM'
        blk_size = 4; %block size
        p = 7;      %searching window size/2
        getBM;
end