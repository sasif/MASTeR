function I1 = operator_OBMC(I0, hor_OB, ver_OB, hor_MC, ver_MC, MASK);
% Performs motion compensation after OBMC (exhaustive BM with sliding block
% window)
% 
% Inputs: 
% I0 : reference image;
% hor_OB, ver_OB : hor. and ver. indices for the sliding block in the target
% image
% hor_MC, ver_MC : hor. and ver. indices in the ref. image
% MASK : Normalizing OB mask for each pixel 

[ROW, COL] = size(I0);

ind_OB = ver_OB+ROW*(hor_OB-1);
ind_MC = ver_MC+ROW*(hor_MC-1);

I1 = accumarray(ind_OB,I0(ind_MC), [ROW*COL 1]);
I1 = reshape(I1, ROW,COL);
I1 = I1./double(MASK);
