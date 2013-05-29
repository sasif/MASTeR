function I0 = adjoint_OBMC(I1, hor_OB, ver_OB, hor_MC, ver_MC, MASK);

% Performs motion compensation after OBMC (exhaustive BM with sliding block
% window)
% 
% Inputs: 
% I0 : reference image;
% hor_OB, ver_OB : hor. and ver. indices for the sliding block in the target
% image
% hor_MC, ver_MC : hor. and ver. indices in the ref. image
% MASK : Normalizing OB mask for each pixel 

[ROW, COL] = size(I1);
I1 = I1./MASK;

ind_OB = ver_OB+ROW*(hor_OB-1);
ind_MC = ver_MC+ROW*(hor_MC-1);

I0 = accumarray(ind_MC,I1(ind_OB), [ROW*COL 1]);
I0 = reshape(I0, ROW,COL);

