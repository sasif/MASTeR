function MASTeR_function(rf, opts)

% Main function for MASTeR 
%
%-------------------------------------------+
% Author: M. Salman Asif @ Georgia Tech
% Email: sasif@gatech.edu
%-------------------------------------------+

global ROW COL T_frames C_coils OMEGA SENSITIVITY_MAPS SAMPLING_MASK
global STATIC_MASK PERIODIC MOTION_FIELD_FORWARD MOTION_FIELD_BACKWARD

DATASET = opts.DATASET; %Change to choose dataset [ 1 2 3 4]
MAX_frames = opts.MAX_frames;
SD = opts.SD; % 0 (no SD regions) or 1 (SD regions as specified in STATIC_MASK)
PRESCAN = opts.PRESCAN;
SM = opts.SM;
reduction_factor = rf;

% Spatial regularization
SP_REG = opts.SP_REG;
L = opts.L;
SOLVER = opts.SOLVER;
CONSTR_SOLVER = opts.CONSTR_SOLVER; 
wave_red = opts.wave_red;

% KT regularization option
KT_REG = opts.KT_REG;   %
KT_NORM = opts.KT_NORM;
KT_TYPE = opts.KT_TYPE; % 'SMOOTH' or 'BRICK'
KEEP_REG = opts.KEEP_REG;

% Motion adaptation
MC_NORM = opts.MC_NORM;
MOTION_TYPE = opts.MOTION_TYPE;
MOTION_DIR  = opts.MOTION_DIR; % 'FOR','BACK', 'Bi', 'wtBi';
PERIODIC = opts.PERIODIC; % Assumes that video sequence is periodic and links first and last frame with motion too
oracleMotion = opts.oracleMotion;

% Apply weighting
RWT = opts.RWT;

% Save results
mname = opts.mname;
sim_dir = opts.sim_dir;

% Want to save the reconstructed images, set it true and setup location
% below
SAVE_RESULTS = opts.SAVE_RESULTS;
DIFF_AMP = 5;

STATIC_MASK = []; % Mask to represent static regions in the video (all nonzero places are static)

% KT Regularization parameters
wtB_LFt = 0.1; % Penalize high frequency in temporal domain, act as a lowpass filter
frac_LFt = 0.4; % Fraction of frequency allowed to be treated as lowpass
%
wtS_LFt = 3; % or 0.01
%
% USE L1 ON TEMPORAL FOURIER FOR INITIAL REGULARIZATION: expt_LFt = 0;
% USE L2 ON WEIGHTED TEMPORAL DFT FOR INITIAL REGULARIZATION: exp_LFt = 1;
switch KT_NORM
    case 'L2'
        exp_LFt = 1; % Smooth weight on fft in temporal: k^exp_LFt, and k = [0:T/2 T/2:-1:1].
        % BRICK_WT (wt_LFt = 0.1, frac_LFt ~ 0.4) or SMOOTH_WT (wt_LFt = 0.01, frac_LFt ~ 0.5)
    case 'L1'
        exp_LFt = 0;
end

% Motion parameters
MOTION = true;
MOTION_Refine_itr = 3*MOTION; % Number of times want to refine reconstruction by re-estimating motion and solving CG again

switch MC_NORM
    case 'L2'
        wt_MC0 = 0.2;   %1e-2*reduction_factor;%reduction_factor/20; 
        % Tradeoff parameter between the motion estimate and the measurement fidelity
        % e.g., if wt_MC = 0 it does not include any motion informatio
        % in the reconstruction
    case 'L1'
        wt_MC0 = 1;
end

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

switch SM
    case 'PINOT'
        MAX_SIM = 1;
    otherwise
        MAX_SIM = 1;
end

SER_stack = [];


stack_ind = 1;
rseed = opts.rseed;
rand('state',rseed);
randn('state',rseed);

% if SAVE_RESULTS
%     diary(sprintf([sim_dir,'/',mname,'_',Filename,'_Sampling-',SM,'_rf%1.3g.txt'],rf))
% end
fprintf(['----------------',mname,':  ',Filename,' @ rf=',num2str(rf),' ',SM,' ',datestr(now),'----------------------------------\n'])

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
    
    if ~SD
        STATIC_MASK = zeros(ROW,COL);
    end
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
    y = subsampled_kspace_vec;
    full_k_space = [];
    subsampled_kspace_vec = [];
    
    % Sensing operator
    A_meas = @(z) A_OMEGA(z);
    At_meas = @(z) At_OMEGA(z);
    
    %% Set wavelet parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Wavelets parameters %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % if isempty(L); L = 4; end
    % remember remember
    % L = log2(n)-J; % Coarsest level in the wavelets (Wavelab)
    % J = log2(n)-L; % Number of levels in wavelet transform (...)
    
    J = floor(log2(min(ROW,COL)))-4; % Number of scales in wavelet transform
    if wave_red == 1;
        % Complex wavelets
        % (Fhh +/- Fgg, Fhg +/- Fgh)
        red = 4;
        
        SYM = 3;
        [Faf, Fsf, af, sf] = BiOrthDualFilt_mod;
        
        % SYM = 0;
        % [Faf, Fsf] = FSfarras; % 1st stage anal. & synth. filters
        % % [Faf, Fsf] = AntonB;
        % [af, sf] = dualfilt1;
        
        C2D = 0;    % 1 -- complex coefficients
                    % 0 -- real and imaginary separate
        
        psiT = @(z) CWT2D_op(SD2seq(z), Faf, af, Fsf, sf, J, SYM, C2D, ROW, COL);
        psi = @(z) seq2SD(adj_CWT2D_op(z, Faf, af, Fsf, sf, J, SYM, C2D, ROW, COL));
        % invpsiT = @(z) iCWT2D_op(z, Faf, af, Fsf, sf, J, SYM, C2D, ROW, COL);
        
        fprintf('Overcomplete wavelets with SYM=%d\n',SYM);
    else
        % % Orthogonal wavelet parameters
        %     QMF = MakeONFilter('Daubechies',4);
        %     L = log2(ROW)-3;
        %     psiT = @(z) FWT2D_op(z, L, QMF, ROW, COL);
        %     psi = @(z) IWT2D_op(z, L, QMF, ROW, COL);
        red = 1;
        [h0, h1, g0, g1] = dauborth(4); sym = 0;
        % [h0, h1, g0, g1] = daub79; sym = 1;
        % [h0, h1, g0, g1] = daub1018; sym = 2;
        psiT = @(z) FWT2D_op(SD2seq(z), h0, h1, J, sym, ROW, COL);
        psi = @(z) seq2SD(adj_FWT2D_op(z, h0, h1, J, sym, ROW, COL));
        % invpsiT = @(z) IWT2D_op(z, g0, g1, J, sym, ROW, COL);
        
        % [af sf] = farras;
        % psiT = @(z) DWT2D_op(SD2seq(z), J, af, ROW, COL);
        % psi = @(z) seq2SD(iDWT2D_op(z, J, sf, ROW, COL));
        fprintf('Complete wavelets with SYM=%d\n',sym);
    end
    
    % %% KT
    % iFkt = @(z) vec(ifft(reshape(z,ROW,COL,T_frames),[],3)*sqrt(T_frames));
    % Fkt = @(z) vec(fft(reshape(z,ROW,COL,T_frames),[],3)/sqrt(T_frames));
    
    % SAMPLING_MASK = repmat(permute(MASK3,[1 2 4 3]), [1,1,C_coils,1]);
    % Af = @(z) Af_SMAP_general(z, SAMPLING_MASK);
    % At = @(z) At_SMAP_general(z, SAMPLING_MASK);
    % AtA_h = @(z) At(Af(z));
    % y = SAMPLING_MASK.*ifft(full_k_space,[],2)*sqrt(COL);
    % y = y(:);
    
    %% RECOVERY ROUTINES
    maxiter = MOTION_Refine_itr+1;
    len_y = length(y); MC_REG = 0; 
    
    for m_iter = 1:maxiter
        % eval(sprintf('load Data%d_rf%d_spreg Ir_cube;',DATASET,rf));
        if m_iter > 1 || oracleMotion==1;
            %% Motion adaptive iterations
            switch oracleMotion
                case 0
                    Ir_cube_me = Ir_cube;
                case 1
                    Ir_cube_me = I_cube;
                    wt_MC0 = 2*wt_MC0;
            end
            
            % [FI_cube BI_cube] = getMOTION(Ir_cube_me,MOTION_TYPE,AR_hor,AR_ver);
            fprintf('Estimating motion using %s method... \n',MOTION_TYPE);
            getMOTION;
            Ir_cube_me = [];
            
            % Motion
            if PERIODIC
                desired_frames = 1:T_frames;
            else
                desired_frames = 2:T_frames-1;
            end
            
            if m_iter == 2
                ALPHA = 0.5;
                ma = ones(ROW*COL*length(desired_frames),1)*0.5;
                mb = ma;
            else
                computeALPHA_MRI;
                ma = vec(repmat(ALPHA(desired_frames,1)',ROW*COL,1));
                mb = vec(repmat(ALPHA(desired_frames,2)',ROW*COL,1));
            end
            
            opts_motion = [];
            opts_motion.mtype = MOTION_TYPE;
            opts_motion.desired_frames = desired_frames;
            if strcmp(MOTION_TYPE,'OBMC')
                opts_motion.OBMC.hor_OB = hor_OB;
                opts_motion.OBMC.ver_OB = ver_OB;
                opts_motion.OBMC.MASK = MASK;
            end
            
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
                    % Bi-directional
                    m_red = 2;
                    motion_length = ROW*COL*length(desired_frames(:));
                    A_motion = @(z) [A_forward_linearMC(SD2seq(z),opts_motion)-SD2seq(z,desired_frames); A_backward_linearMC(SD2seq(z),opts_motion)-SD2seq(z,desired_frames)]/sqrt(2);
                    At_motion = @(z) (seq2SD(At_forward_linearMC(z(1:motion_length),opts_motion))-seq2SD(z(1:motion_length),desired_frames)+seq2SD(At_backward_linearMC(z(motion_length+1:end),opts_motion))-seq2SD(z(motion_length+1:end),desired_frames))/sqrt(2);
                    % A_motion = @(z) [A_motion_SD_forward(z); A_motion_SD_backward(z)]/sqrt(2);
                    % At_motion = @(z) [At_motion_SD_forward(z(1:RCT))+At_motion_SD_backward(z(RCT+1:2*RCT))]/sqrt(2);
                case 'wtBi'
                    % Bi-directional weighted MC
                    m_red = 1;
                    A_motion = @(z) (ma.*A_forward_linearMC(SD2seq(z),opts_motion)+mb.*A_backward_linearMC(SD2seq(z),opts_motion))-SD2seq(z,desired_frames);
                    At_motion = @(z) seq2SD(At_forward_linearMC(ma.*z,opts_motion)+At_backward_linearMC(mb.*z,opts_motion))-seq2SD(z,desired_frames);
            end
            MC_REG = 1;
            KT_REG = 0; % CHANGE THIS IF YOU WANT TO IMPOSE KT_REG IN MOTION STEPS TOO
        end
        %% Setup problem
        
        % define operators for L1 and L2 regularizations
        defineOperators
        
        if L1reg == 0
            % Use conjugate gradient to solve
            SOLVER = 'CG';
            lambda = 0; delta = 0;
            % CG parameters
            cg_tol = 1e-6;
            cg_verbose = 0;  % Display progress after xx iterations with plots of all frames
            cg_maxiter = 100;
            if m_iter == 1
                x_cg = zeros((T_frames*N_D+N_S),1);
            else
                x_cg = seq2SD(Ir_cube(:));
            end
            Aty = At_OMEGA(y);
            AtA_h = @(z) At_h(A_h(z));%+wt_LFt*TemporalLoPass_general(z,kt_MASK);
            fprintf('Solving CG with tol = %3.4g, maxiter=%d...\n',cg_tol,cg_maxiter);
            [x_cg, best_I, resid, niter, residuals, rms_table] = cgsolve2_MRI(x_cg, AtA_h, Aty, cg_tol, cg_maxiter, seq2SD(I_cube(:)), cg_verbose);
            [val ind] = min(rms_table);
            fprintf('bestrms = %3.4g, bestrms_iter = %d. cg_maxiter = %d. \n', val,ind, cg_maxiter);
            
            figure(151);
            subplot(121); plot(residuals);
            subplot(122); plot(rms_table);
            
            Ir_cube = cube(SD2seq(x_cg));
            outData = [];
        else
            opts = [];
            opts.Verbose = 50;
            % ALWAYS USE CONTINUATION; I DIDN'T SEE ANY SPEEDUP WITHOUT CONTINUATION... 
            % (Learned from experiments on 06/29/2012)--SAL
            opts.MaxIntIter = 5; 
            opts.maxiter = 300;
            opts.TolVar = 1e-6;
            opts.stoptest = 1;
            
            if m_iter > 1
                opts.xplug = Ir_cube(:);
            end
            delta = 0.1;      % 0.1
            
            Aty = At_OMEGA(y);
            muf = 1e-6; % 1e-8*max(abs(Aty));    % 1e-6
            
            lambda = 1e-3*max(abs(At_meas(y)))/(2*m_iter*wt_MC0);
            
            % DATASET 2, SP_REG, MC_ONLY
            if DATASET == 2
                lambda = 1e-3*max(abs(At_meas(y)));
            end
            
            La = 10;
            
            switch SOLVER
                case 'NESTA_L1'
                    % U_h = @(z) [alpha*Fkt(z); beta*psiT(z)];
                    % Ut_h = @(z) [alpha*iFkt(z(1:ROW*COL*T_frames))+beta*psi(z(1+ROW*COL*T_frames:end))]; ;
                    
                    % CG for AAtinv
                    % delta = 0;
                    % A_function = @(x) A_h(At_h(x));
                    % cg_tol = 1e-6; cg_maxit = 20;
                    % CGwrapper(); % (first, zero-out the CGwrapper counters)
                    % opts.AAtinv = @(b) CGwrapper(A_function,b,cg_tol,cg_maxit);
                    
                    opts.U = U_h;
                    opts.Ut = Ut_h;
                    opts.typemin = 'L1';
                    opts.errFcn = @(x) norm( abs(x) - abs(I_cube(:))) / norm(I_cube(:));
                    opts.outFcn = @(x) [norm( abs(x) - abs(I_cube(:)), 'inf' ), norm( abs(x) - abs(I_cube(:))) / norm(I_cube(:))];
                    
                    switch CONSTR_SOLVER
                        case 'constr'
                            % Works better with SP as well as KT-regularization
                            fprintf('Solving constrained NESTA using muf=%3.4g, delta=%3.4g, maxiter=%d...\n',muf,delta,opts.maxiter);
                            [xk_NESTA,niter,resid,outData] = NESTA(A_h,At_h,[y; zeros(len_A-len_y,1)],muf,delta,opts);
                        case 'unconstr'
                            % Works better with SP-regularization
                            fprintf('Solving unconstrained NESTA using muf=%3.4g, lambda=%3.4g, La=%1.2g, maxiter=%d...\n',muf,lambda,La,opts.maxiter);
                            [xk_NESTA,niter,resid,outData] = NESTA_UP(A_h,At_h,[y; zeros(len_A-len_y,1)],lambda, La, muf, opts);
                    end

                    Ir_cube = cube(SD2seq(xk_NESTA));
                    
                case 'NESTA_TV'
                    opts.typemin = 'tv';
                    
                    % [xk_NESTA,niter,resid,outData] = NESTA(A_h,At_h,[y; zeros(m_red*RCT,1)],muf,delta,opts);
                    [xk_NESTA,niter,resid,outData,out] = NESTA_UP(A_h,At_h,[y; zeros(len_A-len_y,1)],lambda, La, muf, opts);
                    
                    Ir_cube = cube(SD2seq(xk_NESTA));
                    
                case 'spg'
                    AF_h = @(z) A_h(psi(z));
                    FAt_h = @(z) psiT(At_h(z));
                    %             A_h = @(z) [Af(psi(z)); wt_reg*A_reg(psi(z))];
                    %             At_h = @(z) [psiT(At(z(1:len_y))+wt_reg*At_reg(z(len_y+1:end)))];
                    %             A_h = @(z) Af(psi(z));
                    %             At_h = @(z) psiT(At(z)); m_red = 0;
                    
                    opA = @(x,mode) opA_spg(x,AF_h,FAt_h,mode);
                    
                    opts = spgSetParms('verbosity',0);
                    opts.iterations = 500;
                    %         delta = norm(A_h(psiT(I_cube(:)))-b);
                    %         [x_bpdn,R,G,INFO_bpdn] = spg_bpdn(opA, [y; zeros(RCT,1)], delta, opts);
                    %         Ir_bpdn = cube(psi(x_bpdn));
                    %         clear x_bpdn R G;
                    
                    lambda = sum(abs(psiT(I_cube(:))));
                    % HOW ABOUT
                    % delta = max(abs(At(y)))*RCT
                    % or delta = sum(abs(psiT(Ir_cube_prev(:))));
                    [x_lasso,R,G,INFO_lasso] = spg_lasso(opA, [y; zeros(len_A-len_y,1)], lambda, opts);
                    clear R G;
                    Ir_cube = cube(SD2seq(psi(x_lasso)));
                case 'SALSA'
                    % SALSA constrained solver
                    %%%% algorithm parameters
                    outeriters = 300;
                    sigma = 0;
                    tol = 1e-3;
                    mu1 = 1e-3*max(abs(At_meas(y)))/max(m_iter,1);
                    mu2 = mu1;
                    epsilon = 0.1;
                    A_function = @(x) At_h(A_h(x));
                    cg_tol = 1e-6; cg_maxit = 20;
                    CGwrapper(); % (first, zero-out the CGwrapper counters)
                    invLS = @(b,mu) CGwrapperCSALSA(A_function,mu,b,cg_tol,cg_maxit);
                    PT = @(z) U_h(z);
                    P = @(z) Ut_h(z);
                    [xk_CSALSA, numA, numAt, objective, distance1, distance2, criterion, t_iter, mses] = ...
                        CSALSA_v2([y; zeros(len_A-len_y,1)], A_h, mu1, mu2, sigma,...
                        'AT', At_h, ...
                        'P', P, ...
                        'PT', PT, ...
                        'StopCriterion', 2, ...
                        'True_x', I_cube(:), ...
                        'ToleranceA', tol,...
                        'MAXITERA', outeriters, ...
                        'LS', invLS, ...
                        'VERBOSE', 1, ...
                        'EPSILON', epsilon, ...
                        'CONTINUATIONFACTOR', 5);
                    
                    % SALSA unconstrained solver
                    %%%% algorithm parameters
                    lambda = 0.02;
                    mu = lambda*5;
                    inneriters = 1;
                    outeriters = 10000;
                    tol = 1e-6;
                    A_function = @(x) mu*x+At_h(A_h(x));
                    cg_tol = 1e-6; cg_maxit = 40;
                    CGwrapper(); % (first, zero-out the CGwrapper counters)
                    invLS = @(b) CGwrapper(A_function,b,cg_tol,cg_maxit);
                    PT = @(z) psiT(z);
                    P = @(z) invpsiT(z);
                    
                    [xk_SALSA, numA, numAt, objective, distance, t_itr, mses] = ...
                        SALSA_v2([y; zeros(len_A-len_y,1)], A_h, lambda,...
                        'MU', mu, ...
                        'AT', At_h, ...
                        'StopCriterion', 1, ...
                        'True_x', I_cube(:), ...
                        'ToleranceA', tol,...
                        'MAXITERA', outeriters, ...
                        'P', P, ...
                        'PT', PT, ...
                        'LS', invLS, ...
                        'VERBOSE', 1);
                    
                otherwise
                    disp('');
            end
        end
        
        I_max = 1;%max(abs(I_cube(:)));
        Ir_max = 1;%max(abs(Ir_cube(:)));
        
        vec = @(x) x(:);
        SER = norm(vec(abs(I_cube/I_max)))^2/norm(vec(abs(Ir_cube/Ir_max))-vec(abs(I_cube/I_max)))^2;
        SER_ROI = norm(vec(abs(I_cube(ROI_ver,ROI_hor,:,:)/I_max)))^2/norm(vec(abs(Ir_cube(ROI_ver,ROI_hor,:,:)/Ir_max)-abs(I_cube(ROI_ver,ROI_hor,:,:)/I_max)))^2;
        
        SER_frame = [];
        SER_ROI_frame = [];
        
        for i=1:T_frames
            SER_frame(i) = norm(vec(I_cube(:,:,i))/I_max)^2/norm(vec(abs(I_cube(:,:,i)))/I_max-vec(abs(Ir_cube(:,:,i)))/Ir_max)^2;
            SER_ROI_frame(i) = norm(vec(I_cube(ROI_ver,ROI_hor,i))/I_max)^2/norm(vec(abs(I_cube(ROI_ver,ROI_hor,i)))/I_max-vec(abs(Ir_cube(ROI_ver,ROI_hor,i)))/Ir_max)^2;
        end
        
        % Ir_cube_history(:,:,:,m_iter) = Ir_cube;
        SER_stack{stack_ind, 1} = reduction_factor;
        SER_stack{stack_ind,4*(m_iter-1)+2} = SER;
        SER_stack{stack_ind,4*(m_iter-1)+3} = SER_ROI;
        SER_stack{stack_ind,4*(m_iter-1)+4} = SER_frame;
        SER_stack{stack_ind,4*(m_iter-1)+5} = SER_ROI_frame;
        
        if m_iter == 1
            Ir_cube_init = Ir_cube;
        else
            Ir_cube_mc = Ir_cube;
        end
        disp(sprintf('%s @ Dr = %g, %s sampling with SOLVER = %s, lambda = %g, SD%d, SP%d, KT%d_%s, KR%d, MC-%s_%s, m_iter = %d. SER = %g, SER_ROI = %g',...
            Filename, reduction_factor, SM, SOLVER, lambda, SD, SP_REG, KT_REG, KT_NORM, KEEP_REG, MOTION_DIR, MC_NORM, m_iter, SER, SER_ROI));
        
        if oracleMotion == 1
            Ir_cube_mc = Ir_cube;
            break;
        end
    end
    stack_ind = stack_ind+1;
    
    %%
    
    if SAVE_RESULTS
        I_max = max(abs(I_cube(:)));
        Ir_max = I_max; %max(abs(Ir_cube_init(:)));
        % Ir_max = I_max; %max(abs(Ir_cube_reg(:)));
        Im_max = I_max; %max(abs(Ir_cube_mc(:)));
        
        I_cube = abs(I_cube)/I_max;
        % Ir_cube_spatial = abs(Ir_cube_spatial)/Is_max;
        Ir_cube_init = abs(Ir_cube_init)/Ir_max;
        Ir_cube_mc = abs(Ir_cube_mc)/Im_max;
        
        % Images_ALL = [I_cube I_cube I_cube; Ir_cube_spatial Ir_cube_reg Ir_cube_mc; abs(I_cube-Ir_cube_spatial)*DIFF_AMP abs(I_cube-Ir_cube_reg)*DIFF_AMP abs(I_cube-Ir_cube_mc)*DIFF_AMP];
        Images_ALL = [I_cube I_cube;  Ir_cube_init Ir_cube_mc;  abs(I_cube-Ir_cube_init)*DIFF_AMP abs(I_cube-Ir_cube_mc)*DIFF_AMP];
        % Images_ALL = [I_cube;  Ir_cube_spatial;  abs(I_cube-Ir_cube_spatial)*DIFF_AMP];
        % Main directory to save the data folders in
        directory_name=sprintf([Filename,'_SM-',SM,'_rf%1.3g'],reduction_factor);
        outdir = [sim_dir,'/',directory_name];
        status=mkdir(outdir);
        
        filename_save=[outdir,'/sim_parameters.mat DATASET reduction_factor SM ROW COL T_frames I_coils ROI_ver ROI_hor MOTION_TYPE SOLVER delta wt_MC0 MOTION_Refine_itr PERIODIC DIFF_AMP ALPHA SER_stack rand_state randn_state L J Faf af sf Fsf wave_red PRESCAN KT_NORM MC_NORM I_max niter resid outData'];
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
