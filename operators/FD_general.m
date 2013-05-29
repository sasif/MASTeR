function f_diff = FD_general(f_vec)
%
% penalize all the high frequencies outside loFreq_frac range.
%
global P N ROW COL T_frames C_coils 
global OMEGA SENSITIVITY_MASK MOTION_FIELD STATIC_MASK

static_index = find(STATIC_MASK);
dynamic_index = find(~STATIC_MASK);
N_S = length(static_index);
N_D = ROW*COL-N_S;

fvec_s = f_vec(1:N_S);
fvec_d = f_vec(N_S+1:end);
fs = zeros(N_S,1);
fd = zeros(N_D*T_frames,1);
for frame = 1:T_frames
    switch frame
        case 1
            %             f_frameC = [fvec_s(N_S_down+1:N_S,:); fvec_d((frame-1)*N_D+1:frame*N_D,:); fvec_s(1:N_S_down,:)];
            %             f_frameN = [fvec_s(N_S_down+1:N_S,:); fvec_d((frame)*N_D+1:(frame+1)*N_D,:); fvec_s(1:N_S_down,:)];
            
            f_frame = zeros(ROW, COL);
            f_frame(static_index) = fvec_s;
            f_frame(dynamic_index) = fvec_d((frame-1)*N_D+1:frame*N_D);
            f_frameC = f_frame;
            
            f_frame = zeros(ROW, COL);
            f_frame(static_index) = fvec_s;
            f_frame(dynamic_index) = fvec_d((frame)*N_D+1:(frame+1)*N_D);
            f_frameN = f_frame;
            
            diff_frame = f_frameC-f_frameN;
        case T_frames
            %             f_frameP = [fvec_s(N_S_down+1:N_S,:); fvec_d((frame-2)*N_D+1:(frame-1)*N_D,:); fvec_s(1:N_S_down,:)];
            %             f_frameC = [fvec_s(N_S_down+1:N_S,:); fvec_d((frame-1)*N_D+1:(frame)*N_D,:); fvec_s(1:N_S_down,:)];
            f_frame = zeros(ROW, COL);
            f_frame(static_index) = fvec_s;
            f_frame(dynamic_index) = fvec_d((frame-2)*N_D+1:(frame-1)*N_D);
            f_frameP = f_frame;
            
            f_frame = zeros(ROW, COL);
            f_frame(static_index) = fvec_s;
            f_frame(dynamic_index) = fvec_d((frame-1)*N_D+1:(frame)*N_D);
            f_frameC = f_frame;
            
            diff_frame = f_frameC-f_frameP;
        otherwise
            %previous frame
            f_frameP = [fvec_s(N_S_down+1:N_S,:); fvec_d((frame-2)*N_D+1:(frame-1)*N_D,:); fvec_s(1:N_S_down,:)];
            %current frame
            f_frameC = [fvec_s(N_S_down+1:N_S,:); fvec_d((frame-1)*N_D+1:(frame)*N_D,:); fvec_s(1:N_S_down,:)];
            %next frame
            f_frameN = [fvec_s(N_S_down+1:N_S,:); fvec_d((frame)*N_D+1:(frame+1)*N_D,:); fvec_s(1:N_S_down,:)];
             
            f_frame = zeros(ROW, COL);
            f_frame(static_index) = fvec_s;
            f_frame(dynamic_index) = fvec_d((frame-2)*N_D+1:(frame-1)*N_D);
            f_frameP = f_frame;
            
            f_frame = zeros(ROW, COL);
            f_frame(static_index) = fvec_s;
            f_frame(dynamic_index) = fvec_d((frame-1)*N_D+1:(frame)*N_D);
            f_frameC = f_frame;
            
            f_frame = zeros(ROW, COL);
            f_frame(static_index) = fvec_s;
            f_frame(dynamic_index) = fvec_d((frame)*N_D+1:(frame+1)*N_D);
            f_frameN = f_frame;
           
            diff_frame = 2*f_frameC-f_frameP-f_frameN;
    end
    
    fs_temp = diff_frame(static_index); %[diff_frame(N-N_S_down+1:N,:); diff_frame(1:N_S-N_S_down,:)];
    fd_temp = diff_frame(dynamic_index); %diff_frame(N_S-N_S_down+1:N-N_S_down,:);
    
    fs = fs+fs_temp;
    fd((frame-1)*N_D+1:frame*N_D) = fd_temp;
end
f_diff = [fs; fd];