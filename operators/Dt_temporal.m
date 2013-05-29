% Dt_temporal.m
%
% Differencing-in-time adjoint operator
% Output is vectorized

function v = Dt_temporal(d)

global ROW COL T_frames

vec = @(z) z(:);

Cd = reshape(d, ROW, COL, T_frames-1);
C = zeros(ROW, COL, T_frames);
C(:,:,1) = -Cd(:,:,1);
C(:,:,2:T_frames-1) = diff(-Cd(:,:,1:T_frames-1),1, 3);
C(:,:,T_frames) = Cd(:,:,T_frames-1);
v = vec(C);


