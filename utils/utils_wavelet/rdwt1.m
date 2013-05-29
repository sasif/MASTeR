% rdwt1.m
%
% Performs multiple levels of the redundant discrete wavelet transform.
% Lines the coefficients up in the "tree structure".
% Usage : w = rdwt1(x, h0, h1, nlev)
% x - signal, 1xN
% h0,h1 - filter coefficients
% nlev - number of levels in the decomposition
% w - redundant wavelet transform, nlev+1xN.  w(1,:) are the scaling
% coefficients, w(l+1,:) are the wavelet coefficients for level l.
%
% Written by : Justin Romberg
% Created : 8/10/2001

function w = rdwt1(x, h0, h1, nlev)

N = length(x);
w = zeros(nlev+1, N);

h0i = h0;
h1i = h1;
w(1,:) = x;
for ll = 1:nlev
  m0i = length(h0i);
  m1i = length(h1i);
  zi = floor((m0i-m1i)/4)+1;
  si = w(ll,:);
  w(ll+1,:) = cshift(cconv(si, h0i, floor((m1i+1)/2)-1), zi, 'l');
  w(ll,:) = cconv(si, h1i, floor((m1i+1)/2)-1);
  % make the horizontal positions correpond with shifts
  w(ll,:) = cshift(fliplr(w(ll,:)), 1, 'r');
  % upsample filters
  h0i = reshape([h0i; zeros(1,m0i)], [1 2*m0i]);
  h1i = reshape([h1i; zeros(1,m1i)], [1 2*m1i]);
end
w(nlev+1,:) = cshift(fliplr(w(nlev+1,:)), 1, 'r');
w = flipud(w);

