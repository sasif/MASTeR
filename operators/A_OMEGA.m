function y_vec = A_OMEGA(f_vec)
% general - assign any part of the image static
% 
% This function takes a vector which contains T columns, each of length N,
% corresponding to T frames, and applies MRI sensing operator.
% Sensing operator applies a sensitivity mask on each frame, and takes
% samples at OMEGA frequencies from its centeralized Fourier transform.
%
% y_channel_frame = R_{OMEGA} FFT SENSITIVITY_MAPS f_frame
%
% Inputs:
%     f_vec - A vector of length N_S+TN_D, which represents concatenation of a column of length N
%     at a particular location in T consecutive frames
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
%     y_vec - A vector of length T*C*size(OMEGA) corresponding to all the
%     measurements;
%

global ROW COL T_frames C_coils OMEGA SENSITIVITY_MAPS STATIC_MASK

% ROW = opts.ROW;
% COL = opts.COL;
% T_frames = opts.T_frames;
% C_coils = opts.C_coils;
% OMEGA = opts.OMEGA;
% SENSITIVITY_MAPS = opts.SENSITIVITY_MAPS;
% STATIC_MASK = opts.STATIC_MASK;

y_vec = zeros(ROW*COL*T_frames*C_coils,1);
static_index = find(STATIC_MASK);
dynamic_index = ~STATIC_MASK;
N_S = length(static_index);
N_D = ROW*COL-N_S;

fvec_s = f_vec(1:N_S);
fvec_d = f_vec(N_S+1:end);
y_ind = 0;
for frame = 1:T_frames
    f_frame = zeros(ROW, COL);
    f_frame(static_index) = fvec_s;
    f_frame(dynamic_index) = fvec_d((frame-1)*N_D+1:frame*N_D);
    % f_frame = [fvec_s(N_S_down+1:N_S,:); fvec_d((frame-1)*N_D+1:frame*N_D,:); fvec_s(1:N_S_down,:)];
    for coil = 1:C_coils
        Omega_f = OMEGA{frame}{coil};
    
        S_mask = SENSITIVITY_MAPS(:,:,coil,frame);
        % S_mask = ones(N_S+N_D,1);
        Sf = S_mask(:,1:COL).*f_frame;
        
        % Sf_mod = [Sf(N_S_down+1:end,:); Sf(1:N_S_down,:)];
        Sf_mod = Sf;
        
%         FSf = fftshift(fft(fftshift(Sf_mod,1),[],1),1)/sqrt(ROW); 
        FSf = fft(Sf_mod,[],1)/sqrt(ROW);
        
        % If k_space data is not fftshifted, we don't need it either
%         FSf = fft(Sf_mod)/sqrt(N_S+N_D); 

        % I think there should be another fftshift in the correct
        % implementation of system matrix but PINOT implementation does not
        % have that and therefore the reconstructed image will be modulated form 
        % of their original image. LOOK INTO THIS FURTHER 
        
 
        y_c = FSf(Omega_f,:); % *sqrt(ROW/length(Omega_f));
        M = length(Omega_f);
        y_vec(y_ind+1:y_ind+M*COL) = y_c(:);
        y_ind = y_ind+M*COL; %(coil+(frame-1)*C_coils-1)*M+1:(coil+(frame-1)*C_coils)*M);
    end
end
% y_ind = y_ind-M*COL;
y_vec = y_vec(1:y_ind);

% OLD version
% y_vec = [];
% fvec_s = f_vec(1:N_S,:);
% fvec_d = f_vec(N_S+1:end,:);
% for frame = 1:T_frames
%     f_frame = [fvec_s(N_S_down+1:N_S,:); fvec_d((frame-1)*N_D+1:frame*N_D,:); fvec_s(1:N_S_down,:)];
%     Omega_f = OMEGA{frame};
%     for coil = 1:C_coils
%         S_mask = SENSITIVITY_MAPS(:,:,coil,frame);
%         % S_mask = ones(N_S+N_D,1);
%         Sf = S_mask(:,1:COL).*f_frame;
%         
%         % Sf_mod = [Sf(N_S_down+1:end,:); Sf(1:N_S_down,:)];
%         Sf_mod = Sf;
%         
%         FSf = fftshift(fft(fftshift(Sf_mod,1),[],1),1)/sqrt(N_S+N_D); 
%         % If k_space data is not fftshifted, we don't need it either
% %         FSf = fft(Sf_mod)/sqrt(N_S+N_D); 
% 
%         % I think there should be another fftshift in the correct
%         % implementation of system matrix but PINOT implementation does not
%         % have that and therefore the reconstructed image will be modulated form 
%         % of their original image. LOOK INTO THIS FURTHER 
%         
%  
%         y_c = FSf(Omega_f,:);
%         y_vec = [y_vec; y_c];
%     end
% end