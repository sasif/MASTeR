function MRI_rwt_MASTeR_function(rf, opts)

% Main function for MASTeR with k-t FOCUSS
%
%-------------------------------------------+
% Author: M. Salman Asif @ Georgia Tech
% Email: sasif@gatech.edu
%-------------------------------------------+

global ROW COL T_frames C_coils OMEGA SENSITIVITY_MAPS SAMPLING_MASK
global STATIC_MASK PERIODIC MOTION_FIELD_FORWARD MOTION_FIELD_BACKWARD

DATASET = opts.DATASET; %Change to choose dataset [ 1 2 3 4]
MAX_frames = opts.MAX_frames;
rwt_MAXITER = opts.rwt_MAXITER;
PRESCAN = opts.PRESCAN;

MC_FT = opts.MC_FT; % apply FFT along KT direction (1) or not (0);
MOTION_DIR = opts.MOTION_DIR; 

SM = opts.SM;
reduction_factor = rf;
mname = opts.mname;
sim_dir = opts.sim_dir;

cg_tol = 1e-6;
cg_verbose = 0;  % Display progress after xx iterations with plots of all frames
cg_maxiter = opts.cg_maxiter;
cg_opt = opts.cg_opt; 
% If cg_opt = 1, select the solution at CG iteration that gives smalles RMS value w.r.t the original signal
% otherwise, use the convergence criterion

SAVE_RESULTS = opts.SAVE_RESULTS;
DIFF_AMP = 5;
% Want to save the reconstructed images, set it true and setup location below

STATIC_MASK = []; % Mask to represent static regions in the video (all nonzero places are static)


% Motion parameters
MOTION = true;
MOTION_TYPE = opts.MOTION_TYPE;
MOTION_Refine_itr = 2*MOTION; % Number of times want to refine reconstruction by re-estimating motion and solving CG again
PERIODIC = true; % Assumes that video sequence is periodic and links first and last frame with motion too
oracleMotion = opts.oracleMotion;
oracleWeights = opts.oracleWeights; 

switch DATASET
    % ROI (red boxes in the report)
    case 1 % 'P517127'
        Filename = 'P517127';
    case 2 % 'P501767'
        Filename = 'P501767';
    case 3 % 'wrap003'
        Filename = 'wrap003';
    case 4 % 555
        Filename = '555';
    case 6 %
        Filename = 'P48128';
end

MAX_SIM = 1;

SER_stack = [];

disp(['----------------',mname,' ',Filename,' @ rf=',num2str(rf),' ',SM,' ',datestr(now),'----------------------------------'])

stack_ind = 1;

rseed = opts.rseed;
rand('state',rseed);
randn('state',rseed);

for sim_iter = 1:MAX_SIM
    
    %% DATA Read
    MRI_ReadData;
    
    T_frames = min(MAX_frames,T_frames);
    I_coils = [1:coils];
    C_coils = length(I_coils);
    
    full_k_space = full_k_space(:,:,I_coils,1:T_frames);
    SENSITIVITY_MAPS = SENSITIVITY_MAPS(:,:,I_coils,1:T_frames);
    STATIC_MASK_1 = STATIC_MASK;
    
    vec = @(z) z(:);
    cube = @(z) reshape(z,ROW,COL,T_frames);
    
    RCT = ROW*COL*T_frames;
    
    STATIC_MASK = zeros(ROW,COL);
    static_index = find(STATIC_MASK);
    dynamic_index = find(~STATIC_MASK);
    N_S = length(static_index);
    N_D = ROW*COL-N_S;
    
    I_cube = RecUnderlyingSeq(full_k_space,SENSITIVITY_MAPS);
    
    % I_cube = squeeze(sqrt(sum(abs(ifft2(full_k_space)).^2,3)))*sqrt(ROW*COL);
    % I_cube = I_cube/max(abs(I_cube(:)));
    
    %% Sensing
    rand_state = rand('state');
    randn_state = randn('state');
    MRI_Measure;
    
    % y = repmat(permute(MASK3,[1 2 4 3]), [1,1,C_coils,1]).*ifft(full_k_space,[],2)*sqrt(COL);
    % y = y(:);
    % y = subsampled_kspace_vec;
    % clear full_k_space subsampled_kspace_vec;
    
    %% ReWeighted REGULARIZED LS
    % A_h = @(z) A_OMEGA(z);
    % At_h = @(z) At_OMEGA(z);
    % AtA_h = @(z) At_h(A_h(z));
    
    aFkt = @(z) vec(ifft(reshape(z,ROW,COL,T_frames),[],3)*sqrt(T_frames));
    Fkt = @(z) vec(fft(reshape(z,ROW,COL,T_frames),[],3)/sqrt(T_frames));
    
    % aFkt = @(z) aFkt_SD(z);
    % Fkt = @(z) Fkt_SD(z);
    
    SAMPLING_MASK = repmat(permute(MASK3,[1 2 4 3]), [1,1,C_coils,1]);
    Af = @(z) Af_SMAP_general(z, SAMPLING_MASK);
    At = @(z) At_SMAP_general(z, SAMPLING_MASK);
    y = SAMPLING_MASK.*ifft(full_k_space,[],2)*sqrt(COL);
    y = y(:);
    
    clear full_k_space sensitivity_maps surfacecoilimages subsampled_kspace_vec;
    for inner_group = 1:rwt_MAXITER
        Init_wt = 'LR'; % Low resolution (LR) or back projection (BP)
        kt_sparse = 'vec'; % Vector (vec) or group sparse (gs)
        if inner_group == 1
            switch Init_wt
                case 'LR'
                    % Initialize weights from sum-of-squares of low-resolution coil images.
                    ykt = reshape(y,ROW,COL,C_coils,T_frames);
                    %                 ykt = ifft(kt.*MASK3,[],2)*sqrt(COL);
                    W0=zeros(size(ykt));
                    W0(1:4,:,:,:)=ykt(1:4,:,:,:);
                    W0(end-3:end,:,:,:)=ykt(end-3:end,:,:,:);
                    W0=ifft(W0,[],1)*sqrt(ROW);
                    %         W0 = conj(SENSITIVITY_MAPS).*W0;
                    %         W0 = squeeze((sum(abs(W0),3)));
                    
                    W0 = squeeze(sqrt(sum(abs(W0).^2,3)));
                    %                 W0 = abs(W0);
                    W0=fft(W0,[],3)/sqrt(T_frames);
                    clear ykt;
                case 'BP'
                    % Initialize weights from backprojection image.
                    W0 = fft(reshape(At(y),ROW,COL,T_frames),[],3)/sqrt(T_frames);
            end
            
            switch kt_sparse
                case 'vec'
                    % FOCUSS
                    W=abs(W0).^0.5;
                    % W=W/max(W(:));
                    clear W0;
                case 'gs'
                    % Group sparse
                    W = repmat(sum(abs(W0(:,:,2:end)).^2,3),[1,1,T_frames]).^0.5;
                    % W=W/max(W(:));
            end
        else
            switch kt_sparse
                case 'vec'
                    % FOCUSS
                    W=abs(reshape(FI_r,ROW,COL,T_frames)).^0.5;
                    %W=W/max(W(:));
                    % W(:,:,1) = ones(ROW,COL);
                case 'gs'
                    % Group sparse
                    FI_mat = reshape(FI_r,ROW,COL,T_frames);
                    W = repmat(sum(abs(FI_mat(:,:,2:T_frames)).^2,3),[1,1,T_frames]).^0.5;
                    %W=W/max(W(:));
            end
        end
        W = W(:); Wi = 1;        
        
        lambda = 1e-3*max(W)^2;     
        
        % Change of variables... 
        chg = 1;
        if chg
            % change of variables with Fourier alone...
            % Wi = 1./W; W = 1; 
            
            % change of variables with W and Fourier
            AtA_reg = @(z) W.*(Fkt(At(Af(aFkt(W.*z)))))+lambda*Wi.^2.*z;
            Aty = vec(W.*Fkt(At(y)));
            X0_orig = vec(fft(I_cube,[],3)/sqrt(T_frames))./W(:);
        else
            Wi = 1./W;
            AtA_reg = @(z) At(Af(z))+lambda*(aFkt(Wi.^2.*(Fkt(z))));
            Aty = vec(At(y));
            X0_orig = vec(I_cube);
        end        
        
        if inner_group == 1
            X0 = zeros((T_frames*ROW*COL),1);
        end
        X0_orig = vec(fft(I_cube,[],3)/sqrt(T_frames))./W(:);
        
        fprintf('Solving kt-FOCUSS @ rwt iter. %d, using CG with lambda = %3.4g ... \n', inner_group, lambda);
        [x_cg, best_x, res, iter, residuals, rms_table] = cgsolve2_MRI(X0, AtA_reg, Aty, cg_tol, cg_maxiter, X0_orig, cg_verbose);
        [val ind] = min(rms_table);            
        fprintf('bestrms = %3.4g, bestrms_iter = %d. cg_opt = %d, cg_maxiter = %d. \n', val,ind, cg_opt, cg_maxiter);
        
        % figure(151);
        % subplot(121); plot(residuals);
        % subplot(122); plot(rms_table);
           
        if cg_opt == 0
            best_x = x_cg;
        end
        X0 = best_x;
        if chg
            FI_r = W(:).*best_x;
            I_r = (ifft(reshape(FI_r,ROW,COL,T_frames),[],3))*sqrt(T_frames);
        else
            I_r = reshape(best_x,ROW,COL,T_frames);
            FI_r = fft(I_r,[],3)/sqrt(T_frames);
        end
        
        Ir_cube_init = I_r;
        
        SER_frame = [];
        SER_ROI_frame = [];
        vec = @(x) x(:);
        
        I_max = 1; % max(abs(I_cube(:)));
        Ir_max = 1; % max(abs(Ir_cube_init(:)));
        
        for i=1:T_frames
            SER_frame(i) = norm(vec(I_cube(:,:,i))/I_max)^2/norm(vec(I_cube(:,:,i))/I_max-vec(Ir_cube_init(:,:,i))/Ir_max)^2;
            SER_ROI_frame(i) = norm(vec(I_cube(ROI_ver,ROI_hor,i))/I_max)^2/norm(vec(I_cube(ROI_ver,ROI_hor,i))/I_max-vec(Ir_cube_init(ROI_ver,ROI_hor,i))/Ir_max)^2;
        end
        
        vec = @(x) x(:);
        SER_init = norm(vec(abs(I_cube/I_max)))^2/norm(vec(abs(Ir_cube_init/Ir_max))-vec(abs(I_cube/I_max)))^2;
        SER_ROI_init = norm(vec(abs(I_cube(ROI_ver,ROI_hor,:)/I_max)))^2/norm(vec(abs(Ir_cube_init(ROI_ver,ROI_hor,:)/Ir_max)-abs(I_cube(ROI_ver,ROI_hor,:)/I_max)))^2;
        fprintf('%s @ Dr = %g using %s sampling with rwt CG iter = %d: SER = %g, SER_ROI = %g \n', Filename, reduction_factor, SM, inner_group, SER_init, SER_ROI_init);
        
        % Ir_cube_history(:,:,:,1) = Ir_cube_init;
        m_iter = 0;
        SER_stack{stack_ind, 1} = reduction_factor;
        SER_stack{stack_ind,4*(m_iter)+2} = SER_init;
        SER_stack{stack_ind,4*(m_iter)+3} = SER_ROI_init;
        SER_stack{stack_ind,4*(m_iter)+4} = SER_frame;
        SER_stack{stack_ind,4*(m_iter)+5} = SER_ROI_frame;
    end
    clear W FI_r ;    
    
    %% MOTION INCLUSION
    % ykt = kt(:,:,1:T_frames).*MASK3;
    % recon=KTFOCUSS(ykt);
    % eval(sprintf('load temp_data%d_rf%d.mat;',DATASET, rf));
    
    yo = y;
    for m_iter = 1:MOTION_Refine_itr
        % reference frame for ME/MC
        %ref=kt(:,:,T_frames);
        %ref=ifft2(ref)*ROW;
        fprintf('Motion iteration:%d',m_iter);
        
        %% ME/MC
        if m_iter ==1
            Ir_cube_mc = Ir_cube_init;
        end
        
        switch oracleMotion
            case 0
                Ir_cube_me = Ir_cube_mc;
            case 1
                Ir_cube_me = I_cube;
        end
        
        fprintf(', MOTION_TYPE:%s, OracleMotion:%d, OracleWeights:%d \n',MOTION_TYPE,oracleMotion,oracleWeights);
        
        getMOTION;
        % (Ir_cube_me,MOTION_TYPE,AR_hor,AR_ver);

        Ir_cube_me = [];       
        
        % Motion operators
        desired_frames = 1:T_frames;
        opts_motion = [];
        opts_motion.mtype = MOTION_TYPE;
        opts_motion.desired_frames = desired_frames;
        if strcmp(MOTION_TYPE,'OBMC')
            opts_motion.OBMC.hor_OB = hor_OB;
            opts_motion.OBMC.ver_OB = ver_OB;
            opts_motion.OBMC.MASK = MASK;
        end
        motion_length = ROW*COL*length(desired_frames(:));
        
        switch MOTION_DIR
            case 'FOR'
                % Only forward direction
                m_red = 1;
                A_motion = @(z) A_forward_linearMC(SD2seq(z),opts_motion)-SD2seq(z,desired_frames);
                At_motion = @(z) seq2SD(At_forward_linearMC(z,opts_motion))-seq2SD(z,desired_frames);
            case 'BACK'
                % Only backward direction
                m_red = 1;
                A_motion = @(z) A_backward_linearMC(SD2seq(z),opts_motion)-SD2seq(z,desired_frames);
                At_motion = @(z) seq2SD(At_backward_linearMC(z,opts_motion))-seq2SD(z,desired_frames);
            case 'Bi'
                m_red = 2;
                A_motion = @(z) [A_forward_linearMC(SD2seq(z),opts_motion)-SD2seq(z,desired_frames); A_backward_linearMC(SD2seq(z),opts_motion)-SD2seq(z,desired_frames)]/sqrt(2);
                At_motion = @(z) (seq2SD(At_forward_linearMC(z(1:motion_length),opts_motion))-seq2SD(z(1:motion_length),desired_frames)+seq2SD(At_backward_linearMC(z(motion_length+1:end),opts_motion))-seq2SD(z(motion_length+1:end),desired_frames))/sqrt(2);           
        end
        
        %% Setup reg. on MC differences
                
        % Set Fkt and aFkt to identity if want sparsity in spatial domain 
        % along the KT direction, otherwise apply FFT along KT direction
        if MC_FT
            % Apply temporal Fourier on the MC residuals
            aFkt = @(z) vec(ifft(reshape(z,ROW,COL,m_red*T_frames),[],3)*sqrt(m_red*T_frames));
            Fkt = @(z) vec(fft(reshape(z,ROW,COL,m_red*T_frames),[],3)/sqrt(m_red*T_frames));
            fprintf('Apply Fourier transform on MC residuals\n');
        else
            % No temporal Fourier transform
            fprintf('No FFT on MC residuals\n');
            Fkt = @(z) z; aFkt = @(z) z;
        end
        
        for inner_group = 1:rwt_MAXITER
            % FOCUSS
            if oracleWeights
                W = abs(Fkt(A_motion(I_cube(:)))).^0.5;
            else
                W = abs(Fkt(A_motion(Ir_cube_mc(:)))).^0.5;
            end
            W = W(:);            
            
            if DATASET == 2
                % [lambda=1e-3*max(W)^2 does not work well with dataset 2, use 1e-4 for dataset 2]
                W(abs(W)<max(W)/1e5) = max(W)/1e5;
                lambda = 1e-4;
            else
                W(abs(W)<max(W)/1e5) = max(W)/1e5;
                lambda = 1e-3*max(W)^2;
            end
            
            Wi = 1./W(:);
            AtA_reg = @(z) At(Af(z))+lambda*At_motion(aFkt(Wi.^2.*(Fkt(A_motion(z)))));
            Aty = vec(At(y));
            X0_orig = vec(I_cube);
            
            X0 = Ir_cube_mc(:);            
            if inner_group == 1 && m_iter == 1;
                X0 = 0*X0;
            end
            
            fprintf('Solving rwt Motion adaptive regularization @ rwt iter. %d, using CG with lambda = %3.4g ... \n', inner_group, lambda);
            [x_cg, best_x, res, iter, residuals, rms_table] = cgsolve2_MRI(X0, AtA_reg, Aty, cg_tol, cg_maxiter, X0_orig, cg_verbose);
            
            [val ind] = min(rms_table);
            fprintf('bestrms = %3.4g, bestrms_iter = %d. cg_opt = %d, cg_maxiter = %d. \n', val,ind, cg_opt, cg_maxiter);
            
            % figure(151); 
            % subplot(121); plot(residuals);
            % subplot(122); plot(rms_table);
            
            if cg_opt == 0
                best_x = x_cg;
            end
            
            % Add back the prediction
            Ir_cube_mc = reshape(best_x,ROW,COL,T_frames);
            
            I_max = 1;%max(abs(I_cube(:)));
            Ir_max = 1;%max(abs(Ir_cube_mc(:)));
            
            vec = @(x) x(:);
            
            SER_frame = [];
            SER_frame_ROI = [];
            
            for i=1:T_frames
                SER_frame(i) = norm(vec(I_cube(:,:,i))/I_max)^2/norm(vec(I_cube(:,:,i))/I_max-vec(Ir_cube_mc(:,:,i))/Ir_max)^2;
                SER_ROI_frame(i) = norm(vec(I_cube(ROI_ver,ROI_hor,i))/I_max)^2/norm(vec(I_cube(ROI_ver,ROI_hor,i))/I_max-vec(Ir_cube_mc(ROI_ver,ROI_hor,i))/Ir_max)^2;
            end
            
            SER_mc = norm(vec(abs(I_cube/I_max)))^2/norm(vec(abs(Ir_cube_mc/Ir_max))-vec(abs(I_cube/I_max)))^2;
            SER_ROI_mc = norm(vec(abs(I_cube(ROI_ver,ROI_hor,:)/I_max)))^2/norm(vec(abs(Ir_cube_mc(ROI_ver,ROI_hor,:)/Ir_max)-abs(I_cube(ROI_ver,ROI_hor,:)/I_max)))^2;
            disp(sprintf('%s @ Dr = %g using %s sampling using kt ME/MC iter = %d, m_iter = %d: SER = %g, SER_ROI = %g', Filename, reduction_factor, SM, inner_group, m_iter, SER_mc, SER_ROI_mc));
            % Ir_cube_history(:,:,:,m_iter+1) = Ir_cube_mc;
        end
        SER_stack{stack_ind, 1} = reduction_factor;
        SER_stack{stack_ind,4*(m_iter)+2} = SER_mc;
        SER_stack{stack_ind,4*(m_iter)+3} = SER_ROI_mc;
        SER_stack{stack_ind,4*(m_iter)+4} = SER_frame;
        SER_stack{stack_ind,4*(m_iter)+5} = SER_ROI_frame;
    end
    stack_ind = stack_ind+1;
    %%
    
    if SAVE_RESULTS
        I_max = max(abs(I_cube(:)));
        Ir_max = I_max; %max(abs(Ir_cube_init(:)));
        Im_max = I_max; %max(abs(Ir_cube_mc(:)));
        
        I_cube = abs(I_cube)/I_max;
        Ir_cube_init = abs(Ir_cube_init)/Ir_max;
        Ir_cube_mc = abs(Ir_cube_mc)/Im_max;
        
        Images_ALL = [I_cube I_cube; Ir_cube_init Ir_cube_mc; abs(I_cube-Ir_cube_init)*DIFF_AMP abs(I_cube-Ir_cube_mc)*DIFF_AMP];
        
        % Main directory to save the data folders in
        directory_name=sprintf([Filename,'_SM-',SM,'_rf%1.3g_rwtITER%d'],reduction_factor, rwt_MAXITER);
        outdir = [sim_dir,'/',directory_name];
        status=mkdir(outdir);
        
        STATIC_MASK = STATIC_MASK(:,:,1);
        filename_save=[outdir,'/sim_parameters.mat DATASET reduction_factor SM ROW COL T_frames I_coils ROI_ver ROI_hor cg_maxiter cg_tol MOTION_Refine_itr DIFF_AMP rwt_MAXITER SER_stack rand_state randn_state PRESCAN'];
        eval(['save ', filename_save]);
        % filename_save=[outdir,'/Ir_cube_history.mat Ir_cube_history'];
        % eval(['save ', filename_save]);
        
        Images_ALL = permute(Images_ALL,[1 2 4 3]);
        ChaineName_Image=[outdir, '/', mname,'_',Filename,'_SM-',SM,'_rf',sprintf('%1.3g',reduction_factor),'.gif'];
        imwrite(abs(Images_ALL)*255,ChaineName_Image,'gif');
        Images_ALL = squeeze(Images_ALL);
        for ii = 1:T_frames;
            ChaineName_Image=[outdir, '/', mname,'_',Filename,'_SM-',SM,'_rf',sprintf('%1.3g',reduction_factor),'_frame',sprintf('%0.5d',ii),'.png'];
            imwrite(abs(Images_ALL(:,:,ii)),ChaineName_Image,'png');
        end
    end
end
