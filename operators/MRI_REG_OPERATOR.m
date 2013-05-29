switch KT_TYPE
    case 'BRICK'
        % Fixed flat weights on the KT coeffs.
        mask = ones(T_frames,1);
        mask(1:ceil(T_frames*frac_LFt/2))=0;
        mask(T_frames:-1:T_frames-floor(T_frames*frac_LFt/2)+1)=0;
        
        %         kt_MASK = repmat(mask, [1,ROW,COL]);
        %         kt_MASK = shiftdim(kt_MASK,1);
        %         AtA = @(z) At(Af(z));
        %         AtA_h = @(z) AtA(z) + wtB_LFt*TemporalLoPass_general(z,kt_MASK);
        
        cube = @(z) reshape(z,ROW,COL,T_frames);
        
        wt_LFt = wtB_LFt;
        q_LFt = frac_LFt;
        
        
        kt_MASK = shiftdim(mask,-2);
        A_kt_wt = @(z) vec(bsxfun(@times, fft(cube(SD2seq(z)),[],3)/sqrt(T_frames), kt_MASK));
        At_kt_wt = @(z) seq2SD(vec(ifft(bsxfun(@times, cube(z), kt_MASK),[],3)*sqrt(T_frames)));
        % AtA_h = @(z) At(Af(z)) + wt_LFt*At_kt_wt(A_kt_wt(z));
        
    case 'SMOOTH'
        % Fixed smooth-decaying weights on the KT coeffs.
        mask = ones(T_frames,1);
        mask(1:ceil(T_frames/2)) = [0:ceil(T_frames/2)-1].^exp_LFt;
        mask(T_frames:-1:ceil(T_frames/2)+1) = [1:floor(T_frames/2)].^exp_LFt;
   
        %         kt_MASK = repmat(mask, [1,ROW,COL]);
        %         kt_MASK = shiftdim(kt_MASK,1);
        %         AtA = @(z) At(Af(z));
        %         AtA_h = @(z) AtA(z) + wtS_LFt*TemporalLoPass_general_A(z,kt_MASK);
        
        cube = @(z) reshape(z,ROW,COL,T_frames);
        
        wt_LFt = wtS_LFt;
        q_LFt = exp_LFt;
        
        
        kt_MASK = shiftdim(mask,-2);
        A_kt_wt = @(z) vec(bsxfun(@times, fft(cube(SD2seq(z)),[],3)/sqrt(T_frames), kt_MASK));
        At_kt_wt = @(z) seq2SD(vec(ifft(bsxfun(@times, cube(z), kt_MASK),[],3)*sqrt(T_frames)));
        % AtA_h = @(z) At(Af(z)) + wt_LFt*At_kt_wt(A_kt_wt(z));
        
    case 'RWT'
        % Select weights according to the KT coeffs. 
        
end
