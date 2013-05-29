function f_trunc = TemporalLoPass_general(f_vec, kt_MASK)
%
% penalize all the high frequencies outside loFreq_frac range.
%
global COL ROW T_frames
global STATIC_MASK

static_index = find(STATIC_MASK);
dynamic_index = find(~STATIC_MASK);
N_S = length(static_index);
N_D = ROW*COL-N_S;

fvec_s = f_vec(1:N_S);
fvec_d = f_vec(N_S+1:end);
fs = zeros(N_S,1);
fd = zeros(N_D*T_frames,1);

I_cube = zeros(ROW,COL, T_frames);
for frame = 1:T_frames
    f_frame = zeros(ROW, COL);
    f_frame(static_index) = fvec_s;
    f_frame(dynamic_index) = fvec_d((frame-1)*N_D+1:frame*N_D);
    
    I_cube(:,:,frame) = f_frame;
end
FI_cube = fft(I_cube,[],3);

FI_trunc = (kt_MASK.^2).*FI_cube;
I_trunc = ifft(FI_trunc,[],3);
for frame = 1:T_frames
    f_temp = I_trunc(:,:,frame);
    fs_temp = f_temp(static_index); %[f_temp(ROW-N_S_down+1:ROW,:); f_temp(1:N_S-N_S_down,:)];
    fd_temp = f_temp(dynamic_index); %f_temp(N_S-N_S_down+1:ROW-N_S_down,:);
    
    fs = fs+fs_temp;
    fd((frame-1)*N_D+1:frame*N_D) = fd_temp;
end
f_trunc = [fs; fd];