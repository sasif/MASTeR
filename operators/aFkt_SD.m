function f_vec = aFkt_SD(fk_vec)
% adjoint of Fkt_SD

global ROW COL T_frames
global STATIC_MASK
% N is ROW
% P is COL

F_cube = ifft(reshape(fk_vec,ROW,COL,T_frames),[],3)*sqrt(T_frames);

static_index = find(STATIC_MASK);
dynamic_index = find(~STATIC_MASK);

N_S = length(static_index);
N_D = ROW*COL-N_S;

fvec_s = zeros(N_S,1);
fvec_d = zeros(N_D*T_frames,1);

for frame = 1:T_frames
    f_ct = F_cube(:,:,frame);
    fs = f_ct(static_index); 
    fd = f_ct(dynamic_index);
    
    fvec_s = fvec_s+fs;
    fvec_d((frame-1)*N_D+1:frame*N_D) = fd;
end
f_vec = [fvec_s; fvec_d];
