function I1 = box_forward(I0,hor_ind, ver_ind);

global OMEGA N M J Faf Fsf af sf jj L QMF ROW COL

I0_vec = I0(:);
prev_ind = ver_ind(:)+(hor_ind(:)-1)*ROW;
I1_vec = I0_vec(prev_ind);
I1 = reshape(I1_vec,ROW,COL);

end