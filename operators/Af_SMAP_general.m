function y_vec = Af_SMAP_general(f_vec, SAMPLING_MASK)
% This function takes a vector which contains T columns, each of length N,
% corresponding to T frames, and applies MRI sensing operator.
% Sensing operator applies a sensitivity mask on each frame, and takes
% samples at OMEGA frequencies from its centeralized Fourier transform.
%
% y_channel_frame = R_{OMEGA} FFT SENSITIVITY_MAPS f_frame

global ROW COL T_frames C_coils SENSITIVITY_MAPS 

F_seq = reshape(f_vec, ROW, COL, 1, T_frames);
F_seq = repmat(F_seq, [1, 1, C_coils, 1]);
SF = SENSITIVITY_MAPS.*F_seq;
FSF = fft(SF, [],1)/sqrt(ROW);
%Y_seq = repmat(permute(MASK3,[1 2 4 3]), [1,1,C_coils,1]).*FSF;
Y_seq = SAMPLING_MASK.*FSF;
y_vec = Y_seq(:);
