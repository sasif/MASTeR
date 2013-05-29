function Ic_cube = RecUnderlyingSeq(full_k_space,SENSITIVITY_MAPS)

[ROW COL C_coils T_frames] = size(full_k_space);

I_coils = ifft2(full_k_space)*sqrt(ROW*COL);
I_coils = conj(SENSITIVITY_MAPS).*I_coils;
I_coils = sum(I_coils,3);
I_coils = squeeze(I_coils);

S2_diag = sum(abs(SENSITIVITY_MAPS).^2,3);
S2_diag = squeeze(S2_diag);

Ic_cube = I_coils./S2_diag;