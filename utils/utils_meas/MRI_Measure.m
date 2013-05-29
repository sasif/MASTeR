% SM - Sampling mode
% options: PINOT, RAND, HYBRID, ...
% PINOT is deterministice but does not work as well as RAND or HYBRID

OMEGA = [];
matrixdisplay = zeros(ROW,10*T_frames);
if ~exist('emptyROWs')
    emptyROWs = 0;
end
switch SM
    case 'RAND'
        RAND_LO_FRAC = 0;
        for frame=1:T_frames
            if RAND_LO_FRAC
                lo_frac = RAND_LO_FRAC;
                K_row = round((ROW-emptyROWs)/reduction_factor);
                K_lo = ceil(K_row*lo_frac/2);
                K_hi = floor(K_row*lo_frac/2); 
                center_index = setxor(K_lo+1:ROW-K_hi,ROW/2-emptyROWs/2+1:ROW/2+emptyROWs/2);
                Q_c = center_index(randperm(length(center_index)));
                matrixdisplay([1:K_lo ROW:-1:ROW-K_hi+1 Q_c(1:K_row-K_lo-K_hi)]',frame*10)=1;

%                 hi_index = [1:round(ROW/2)-floor(ROW/lo_frac/2) round(ROW/2)+floor(ROW/lo_frac/2)+1:ROW]';
%                 lo_index = [round(ROW/2)-floor(ROW/lo_frac/2)+1:round(ROW/2)-1 round(ROW/2)+1:round(N/2)+floor(ROW/lo_frac/2)]';
%                 Q_lo = lo_index(randperm(length(lo_index)));
%                 Q_hi = hi_index(randperm(length(hi_index)));
%                 matrixdisplay([round(ROW/2); Q_lo(1:K_lo); Q_hi(1:K_hi)],frame*10)=1;
            else
                K_row = round((ROW-emptyROWs)/reduction_factor);
                
%                 Q_row = randperm(ROW-emptyROWs)+emptyROWs/2;
%                 matrixdisplay(Q_row(1:K_row),frame*10)=1;
%                 matrixdisplay(Q_row(K_row),frame*10) = matrixdisplay(ROW/2+1,frame*10);
%                 matrixdisplay(ROW/2+1,frame*10) = 1;
                Q_row = randperm(ROW-1)+1;
                matrixdisplay(Q_row(1:K_row-1),frame*10)=1;
                matrixdisplay(1,frame*10)= 1;
            end
        end
        % matrixdisplay = fftshift(matrixdisplay,1);
        subsampled_kspace_vec = [];
        SAMP_MASK = logical(matrixdisplay(:,10:10:T_frames*10));
        MASK3 = permute(repmat(SAMP_MASK,[1,1,COL]),[1 3 2]);
        for frame=1:T_frames
            rows_in_this_frame=find(matrixdisplay(:,frame*10));
            for antena=1:C_coils
                kspacefromonecoilinoneframe=full_k_space(rows_in_this_frame,:,antena,frame);
                % kspacefromonecoilinoneframe=ifft(fftshift(kspacefromonecoilinoneframe(:,:),2),[],2);% IFFT IN HORIZONTAL DIRECTION AND FFTSHIFT (LIKE ONLY SENSE)
                % kspacefromonecoilinoneframe=fftshift(kspacefromonecoilinoneframe,2);
                
                kspacefromonecoilinoneframe=ifft(kspacefromonecoilinoneframe,[],2)*sqrt(COL);
                
                %subsampled_k_spaces=[subsampled_k_spaces ; kspacefromonecoilinoneframe];    % I PUT ALL THE K SPACES, ONE UNDER THE OTHER
                subsampled_kspace_vec = [subsampled_kspace_vec; kspacefromonecoilinoneframe(:)];
                OMEGA{frame}{antena} = rows_in_this_frame;
            end
        end
        %     case 'RAND_C'
        %         subsampled_kspace_vec = [];
        %         K_row = round(N/reduction_factor);
        %         for frame=1:T_frames
        %             for antena=1:C_coils
        %                 Q_row = randperm(N);
        %                 rows_in_this_coil = Q_row(1:K_row);
        %
        %                 kspacefromonecoilinoneframe=full_k_space(rows_in_this_coil,:,antena,frame);
        %                 % kspacefromonecoilinoneframe=ifft(fftshift(kspacefromonecoilinoneframe(:,:),2),[],2);% IFFT IN HORIZONTAL DIRECTION AND FFTSHIFT (LIKE ONLY SENSE)
        %                 % kspacefromonecoilinoneframe=fftshift(kspacefromonecoilinoneframe,2);
        %
        %                 kspacefromonecoilinoneframe=ifft(kspacefromonecoilinoneframe,[],2)/sqrt(COL);
        %
        %                 %subsampled_k_spaces=[subsampled_k_spaces ; kspacefromonecoilinoneframe];    % I PUT ALL THE K SPACES, ONE UNDER THE OTHER
        %                 subsampled_kspace_vec = [subsampled_kspace_vec; kspacefromonecoilinoneframe(:)];
        %                 OMEGA{frame}{antena} = rows_in_this_coil;
        %             end
        %         end
    case 'HYBRID'
        [MASK3 SAMP_MASK]= Hybrid_DownsamplingMASK(ROW, COL, T_frames, reduction_factor);
        MASK3=logical(ifftshift(MASK3,1));
        SAMP_MASK=logical(ifftshift(SAMP_MASK,1));
        subsampled_kspace_vec = [];
        for frame=1:T_frames
            rows_in_this_frame=find(SAMP_MASK(:,frame));
            for antena=1:C_coils
                kspacefromonecoilinoneframe=full_k_space(rows_in_this_frame,:,antena,frame);
                % kspacefromonecoilinoneframe=ifft(fftshift(kspacefromonecoilinoneframe(:,:),2),[],2);% IFFT IN HORIZONTAL DIRECTION AND FFTSHIFT (LIKE ONLY SENSE)
                % kspacefromonecoilinoneframe=fftshift(kspacefromonecoilinoneframe,2);
                
                kspacefromonecoilinoneframe=ifft(kspacefromonecoilinoneframe,[],2)*sqrt(COL);
                %subsampled_k_spaces=[subsampled_k_spaces ; kspacefromonecoilinoneframe];    % I PUT ALL THE K SPACES, ONE UNDER THE OTHER
                subsampled_kspace_vec = [subsampled_kspace_vec; kspacefromonecoilinoneframe(:)];
                OMEGA{frame}{antena} = rows_in_this_frame;
            end
        end
    case 'PINOT'
        %------------------------------------------
        % code for PINOT provided by L. Hamilton
        %------------------------------------------
        reduction_factor_spacerip = round(reduction_factor/(ROW*COL/(N_S/T_frames+N_D)));
        reduction_factor = reduction_factor_spacerip*(ROW*COL/(N_S/T_frames+N_D));
        %%%%% WHICH ROWS ARE WE GOING TO TAKE FROM THE ORIGINAL K_SPACE FOR SENSE STEP ( 1 in 2, 1 in 3 or 1 in 4)
        [rows_to_be_sampled] = samplerowsenglish(full_k_space, reduction_factor_spacerip);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%% SPACE_RIP + NOQUIST RECONSTRUCTION   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        %%%% FROM THE LINES SELECTED ABOVE IN rows_to_be_sampled, I TAKE SOME OF THEM (NOQUIST STEP)
        [RowsPerFrame, matrixdisplay, rows_to_be_samplednew]=samplefinalrowsenglish(static, staticup, dynamic, reduction_factor_spacerip, rows_to_be_sampled,T_frames);
        % THE INTERESTING OUTPUT IS matrixdisplay. IT HAS SO MANY COLUMNS AS FRAMES YOU WANT TO USE, WHERE I PUT A "1" IN THE LINES I WANT TO TAKE
        % IN THIS FRAME
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        clear rows_to_be_samplednew
        clear RowsPerFrame
        % matrixdisplay gives the indices of rows in the Fourier plane which will be sampled
        % and used for the reconstruction
        
        %%%%% SUBSAMPLING
        % PUT THE KSPACES ONE UNDER THE OTHER
        % subsampled_k_spaces=[];                                                     % IN THIS MATRIX I'M GOING TO PUT ALL THE K SPACES
        subsampled_kspace_vec = [];
        matrixdisplay = fftshift(matrixdisplay,1);
        SAMP_MASK = logical(matrixdisplay(:,10:10:T_frames*10));
        MASK3 = permute(repmat(SAMP_MASK,[1,1,COL]),[1 3 2]);
        
        for frame=1:T_frames
            rows_in_this_frame=find(matrixdisplay(:,frame*10));
            for antena=1:C_coils
                kspacefromonecoilinoneframe=full_k_space(rows_in_this_frame,:,antena,frame);
                % kspacefromonecoilinoneframe=ifft(fftshift(kspacefromonecoilinoneframe(:,:),2),[],2);% IFFT IN HORIZONTAL DIRECTION AND FFTSHIFT (LIKE ONLY SENSE)
                % kspacefromonecoilinoneframe=fftshift(kspacefromonecoilinoneframe,2);
                
                kspacefromonecoilinoneframe=ifft(kspacefromonecoilinoneframe,[],2)*sqrt(COL);
                
                %subsampled_k_spaces=[subsampled_k_spaces ; kspacefromonecoilinoneframe];    % I PUT ALL THE K SPACES, ONE UNDER THE OTHER
                subsampled_kspace_vec = [subsampled_kspace_vec; kspacefromonecoilinoneframe(:)];
                OMEGA{frame}{antena} = rows_in_this_frame;
            end
        end
        
    otherwise
        error('Not a valid option');
end

