function f_vec = At_SMAP_general(y_vec, SAMPLING_MASK)
% Adjoint for Af_SMAP_general

global ROW COL T_frames C_coils SENSITIVITY_MAPS 

Y_seq = reshape(y_vec, ROW,COL,C_coils,T_frames);
% Y_seq = repmat(permute(MASK3,[1 2 4 3]), [1,1,C_coils,1]).*Y_seq;
Y_seq = SAMPLING_MASK.*Y_seq;
iFSF = ifft(Y_seq, [],1)*sqrt(ROW);
iSF = conj(SENSITIVITY_MAPS).*iFSF;
F_seq = sum(iSF,3);
f_vec = F_seq(:);
