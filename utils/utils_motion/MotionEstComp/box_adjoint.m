function I0 = box_adjoint(I1,hor_ind, ver_ind);

global ROW COL

I1_vec = I1(:);
prev_ind = ver_ind(:)+(hor_ind(:)-1)*ROW;

I0_vec = accumarray(prev_ind,I1_vec, [ROW*COL 1]);
I0 = reshape(I0_vec,ROW,COL);

% % accumarray performs the following task
% for ii = 1:ROW*COL;
%     I0_vec(ii) = sum(I1_vec(prev_ind==ii));
% end
% if ~isempty(find(I0_temp-I0_vec))
%     stp = 1;
% end