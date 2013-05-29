% %% Compute ALPHA
% % Use multi-hypothesis to decide weighting of forward and backward
% % predictions.
% % For the case of CS use measurements instead of images (why? JL lemma :-)
% % [cite Video CS with multihypothesis by Trame and Fowler])

% % Ir_cube_me = abs(I_cube);
% ALPHA = zeros(T_frames,2);
% for frame = 1:T_frames
%     x_mat = abs([vec(FI_cube(:,:,frame)) vec(BI_cube(:,:,frame))]);
%     alpha = x_mat\vec(Ir_cube_me(:,:,frame));%inv(x_mat'*x_mat)*(x_mat'*vec(Ir_cube_me(:,:,frame)));
%     ALPHA(frame,:) = alpha';
% end
% ALPHA_orig = ALPHA;

%% ALPHA computation with CS measurements
% global OMEGA SENSITIVITY_MAPS
if 1
    vec = @(z) z(:); cube = @(z) reshape(z,ROW,COL,T_frames);
    SAMPLING_MASK = repmat(permute(MASK3,[1 2 4 3]), [1,1,C_coils,1]);
    Af_seq = @(z) Af_SMAP_general(z, SAMPLING_MASK);
    % Af_seq = @(z) vec(MASK3.*fft(cube(abs(z)),[],2));
    y_f = Af_seq(FI_cube(:));
    y_b = Af_seq(BI_cube(:));
    y_orig = Af_seq(I_cube(:));
    y_orig = 0*y_orig;
    y_orig(SAMPLING_MASK) = y;
    
    ALPHA = zeros(T_frames,2);
    M_ind = 0;
    for frame = 1:T_frames
        M = ROW*COL;
        y_mat = [y_f(M_ind+1:M_ind+M) y_b(M_ind+1:M_ind+M)];
        alpha = y_mat\y_orig(M_ind+1:M_ind+M);%inv(y_mat'*y_mat)*(y_mat'*y_orig(M_ind+1:M_ind+M));
        ALPHA(frame,:) = alpha';
        M_ind = M_ind+M;
    end
    % ALPHA = abs(ALPHA);
    ALPHA = abs(ALPHA)./repmat(sum(abs(ALPHA),2),1,2);
    clear y_f y_b y_orig SAMPLING_MASK;
end
