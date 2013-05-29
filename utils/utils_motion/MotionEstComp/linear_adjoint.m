function I0 = linear_adjoint(I1,hor_ind, ver_ind, ROW, COL);

% global ROW COL

I0_vec = zeros(ROW*COL,1);

% Interpolated locations in vectorized form
ver_floor = min(floor(ver_ind(:)),ROW-1);
ver_ceil = ver_floor+1;
hor_floor = min(floor(hor_ind(:)),COL-1);
hor_ceil = hor_floor+1;

% Consider a box of pixels surrounding the interpolating pixel
int_ind1 = ver_floor + ROW*(hor_floor-1); % top left corner
int_ind2 = ver_ceil + ROW*(hor_floor-1); % bottom left corner
int_ind3 = ver_floor + ROW*(hor_ceil-1); % top right corner
int_ind4 = ver_ceil + ROW*(hor_ceil-1); % bottom right corner

% 
I1_temp = I1(:);
I1_left_t = I1_temp.*(hor_ceil-hor_ind(:));
I1_right_t = I1_temp.*(hor_ind(:)-hor_floor);

%
I1_ind1_t = I1_left_t.*(ver_ceil-ver_ind(:));
I1_ind2_t = I1_left_t.*(ver_ind(:)-ver_floor);
I1_ind3_t = I1_right_t.*(ver_ceil-ver_ind(:));
I1_ind4_t = I1_right_t.*(ver_ind(:)-ver_floor);

% Accumulation
I0_vec1 = accumarray(int_ind1,I1_ind1_t, [ROW*COL 1]);
I0_vec2 = accumarray(int_ind2,I1_ind2_t, [ROW*COL 1]);
I0_vec3 = accumarray(int_ind3,I1_ind3_t, [ROW*COL 1]);
I0_vec4 = accumarray(int_ind4,I1_ind4_t, [ROW*COL 1]);

I0_temp = I0_vec1+I0_vec2+I0_vec3+I0_vec4;
I0 = reshape(I0_temp,ROW,COL);
