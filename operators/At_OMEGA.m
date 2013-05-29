function f_vec = At_OMEGA(y_vec)
% general - assign any part of the image static
%
% Inputs:
%     y_vec - A vector of length T*C*size(OMEGA) corresponding to all the
%     measurements;
%     OMEGA{frame} - Sampling freq. locations for all frames
%     SENSITIVITY_MAPS{frame}{channel} - Sensitivity mask for each coil at every frames
%     T - total number of frames
%     N - lenght of each column in the image frame
%     N_S - length of static part (this is assumed to be the top portion of
%     the vector fvec.
%     N_D - length of dynamic part
%     C - number of coils (channels)
%
% Outputs:
%     f_vec - A vector of length N_S+TN_D, which represents concatenation of a column of length N
%     at a particular location in T consecutive frames
%

global ROW COL T_frames C_coils OMEGA SENSITIVITY_MAPS STATIC_MASK

% ROW = opts.ROW;
% COL = opts.COL;
% T_frames = opts.T_frames;
% C_coils = opts.C_coils;
% OMEGA = opts.OMEGA;
% SENSITIVITY_MAPS = opts.SENSITIVITY_MAPS;
% STATIC_MASK = opts.STATIC_MASK;

static_index = find(STATIC_MASK);
dynamic_index = ~STATIC_MASK;
N_S = length(static_index);
N_D = ROW*COL-N_S;

fvec_s = zeros(N_S,1);
fvec_d = zeros(N_D*T_frames,1);
y_ind = 0;
for frame = 1:T_frames
    fs_temp = zeros(N_S,1);
    fd_temp = zeros(N_D,1);
    
    for coil = 1:C_coils
        Omega_f = OMEGA{frame}{coil};
        M = length(Omega_f);

        y_c = reshape(y_vec(y_ind+1:y_ind+M*COL),M,COL); 
        y_ind = y_ind+M*COL; %(coil+(frame-1)*C_coils-1)*M+1:(coil+(frame-1)*C_coils)*M);
        y_c_mod = zeros(ROW,COL);
        y_c_mod(Omega_f,:) = y_c; % *sqrt(ROW/length(Omega_f));

%         iFsf = fftshift(ifft(fftshift(y_c_mod,1),[],1),1)*sqrt(ROW);
        iFsf = ifft(y_c_mod,[],1)*sqrt(ROW);
        
        S_mask = SENSITIVITY_MAPS(:,:,coil,frame);
        f_ct = conj(S_mask(:,1:COL)).*iFsf;
        
        %         fs = f_ct(SENSITIVITY_INDEX(1:N_S),:);
        %         fd = f_ct(SENSITIVITY_INDEX(N_S+1:end),:);
        
        % fs = f_ct(static_index); %[f_ct(N-N_S_down+1:N,:); f_ct(1:N_S-N_S_down,:)];
        % fd = f_ct(dynamic_index);
        
        
        fs_temp = fs_temp+f_ct(static_index);
        fd_temp = fd_temp+f_ct(dynamic_index);
    end
    fvec_s = fvec_s+fs_temp;
    fvec_d((frame-1)*N_D+1:frame*N_D) = fd_temp;
end
f_vec = [fvec_s; fvec_d];
