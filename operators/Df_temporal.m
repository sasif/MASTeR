% Df_temporal.m
%
% Differencing-in-time forward operator
% Input is vectorized

function d = Df_temporal(x,order)

global ROW COL T_frames

vec = @(z) z(:);
cube = @(z) reshape(z, ROW, COL, T_frames);

C = cube(x);
Cd = diff(C, 1, 3);
d = vec(Cd);
