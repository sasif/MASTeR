function fk_vec = Fkt_SD(f_vec)
% Takes FFT of image cube along temporal direction.

global ROW COL T_frames
global STATIC_MASK
% N is ROW
% P is COL

static_index = find(STATIC_MASK);
dynamic_index = find(~STATIC_MASK);
N_S = length(static_index);
N_D = ROW*COL-N_S;

fvec_s = f_vec(1:N_S);
fvec_d = f_vec(N_S+1:end);

F_cube = zeros(ROW,COL,T_frames);

for frame = 1:T_frames
    f_frame = zeros(ROW, COL);
    f_frame(static_index) = fvec_s;
    f_frame(dynamic_index) = fvec_d((frame-1)*N_D+1:frame*N_D);
    F_cube(:,:,frame) = f_frame;
end
fk_vec = fft(F_cube,[],3)/sqrt(T_frames);
fk_vec = fk_vec(:);
