% ifwt2_CWT.m
%
% Inverse fast 2D wavelet transform using Complex Wavelets
% Here it takes inverse of one pair of output from the CWT2 -- inv(Wxx)

% In iCWT2.m we calculate repective inverses of all Whh, Wgg, Whg and Wgh 
% seperately and then add/subtract accordingly to get the image back 

% Usage: x = ifwt2_CWT(w, g0_r, g1_r, g0_c, g1_c, J, symr, symc, rosym_flip, cosym_flip);
% w - nxn input image of wavelet coeffients
% g0_r - lowpass reconstruction filter on rows
% g1_r - highpass reconstruction filter on rows
%      g0 and g1 should be zero padded appropriately so that they have the
%      same length, and this length is even.
% g0_c - lowpass reconstruction filter on columns
% g1_c - highpass reconstruction filter on columns

% J - number of levels in the filter bank.
%     Right now, it must be chose so that n*2^(-J) >= length(g0)
% symc, symr - How the input was extended when the forward transform was
% taken on rows and columns
%     sym=0: periodization
%     sym=1: type-I symmetric extension ( [... x(2) x(1) x(2) x(3) ...])
%            The wavelet filters must have type-I even symmetry 
%            (e.g. daub79)
%     sym=2: type-II symmetric extension ( [... x(2) x(1) x(1) x(2) x(3) ...])
%            The lowpass filter must have type-II even symmetry, 
%            The highpass filter must have type-II odd symmetry.
%            (e.g. daub1018)
% rosym_flip, cosym_flip - 0 : do nothing
%                          1 : flip the symmetric extensions for row (r) or
%                          column (c) from Es(1,2) on scaling to Es(2,1)
%                          and Es(2,1) or wavelet to Es(1,2).
% Written by: Salman Asif
% Created: February 2008
