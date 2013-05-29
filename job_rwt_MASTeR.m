% This code implements kt FOCUSS with MASTeR
% in which we replaced reference frame-based motion-compensated residuals 
% terms of kt FOCUSS with ME/MC with the motion-adaptive transforms from MASTeR
%
% For further details see the Discussion section in our paper 
% 
%-------------------------------------------+
% Author: M. Salman Asif @ Georgia Tech
% Email: sasif@gatech.edu
%-------------------------------------------+

%% Path/machine setup
%
% numCores = str2double(getenv('NUMBER_OF_PROCESSORS'));
% numCores = feature('numCores');
% mpSize = 4;
% if matlabpool('size') ~= mpSize
%     if matlabpool('size')~=0
%         matlabpool close;
%     else
%         matlabpool('open', mpSize);
%     end
% end
 
mname = mfilename;
mpath = mfilename('fullpath');
mdir = mpath(1:end-length(mname));
cd(mdir); 

setup_path;

%% Setup simulation parameters

global ROW COL T_frames C_coils OMEGA SENSITIVITY_MAPS SAMPLING_MASK 
global STATIC_MASK PERIODIC MOTION_FIELD_FORWARD MOTION_FIELD_BACKWARD

DATASET = 2; PRESCAN = true;
% DATASET = 3; PRESCAN = false;
MAX_frames = 16;

rwt_MAXITER = 3; % maximum iteration for reweighting
cg_maxiter = 200; % maximum conjugate gradient iterations at every reweighting step
cg_opt = 0;     % 1 -- select the solution of CG that yields minimum MSE
                % 0 -- fixed threshold for the termination criterion

MOTION_TYPE = 'OBMC'; % {'CWT','OF','OBMC'}; % scheme for motion estimation
oracleMotion = 0;  % 0 -- estimate motion from reconstructed frames
                    % 1 -- estimate motion from original frames (oracle)
oracleWeights = 0;  % weights selection for the iterative reweighting

MC_FT = 0; % apply FFT along KT direction (1) or not (0);
MOTION_DIR = 'Bi'; % type for motion compensation, Bi, wtBi, FOR, BACK

% sampling scheme 
SampType = {'HYBRID'}; % HYBRID: few low-frequency + randomly sampled high frequency horizontal k-space lines
                       % RAND: randomly subsample horizontal k-space lines

% reduction factor 
rf_table = 8; % [2:2:12];

rseed = 0;
SAVE_RESULTS = false;

opts = [];
opts.DATASET = DATASET; 
opts.MAX_frames = MAX_frames;
opts.rwt_MAXITER = rwt_MAXITER; 
opts.MOTION_TYPE = MOTION_TYPE;
opts.oracleMotion = oracleMotion;
opts.oracleWeights = oracleWeights;
opts.MC_FT = MC_FT;
opts.MOTION_DIR = MOTION_DIR;
opts.PRESCAN = PRESCAN;
opts.cg_maxiter = cg_maxiter;
opts.cg_opt = cg_opt; 
opts.rseed = rseed;
opts.SAVE_RESULTS = SAVE_RESULTS;

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

sim_dir = sprintf(['results/',mname,'_pre%d_om%d_ow%d_MOTION-%s_MC_FT%d_cg_opt%d_cg_maxiter%d'],PRESCAN,oracleMotion,oracleWeights,MOTION_TYPE,MC_FT,cg_opt,cg_maxiter);

opts.mname = mname;
opts.sim_dir = sim_dir;

for ps = 1:length(SampType)

    SM = char(SampType(ps));
    opts.SM = SM;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if SAVE_RESULTS
        status = mkdir(sim_dir);
        diary(sprintf([sim_dir,'/',mname,'_',Filename,'_Sampling-',SM,'.txt']))
    end
    
    fprintf('data-%d, pre-%d, cg_maxiter-%d, rwt_iter-%d, MOTION-%s, om-%d, ow-%d\n',DATASET,PRESCAN,cg_maxiter, rwt_MAXITER, MOTION_TYPE, oracleMotion,oracleWeights);

    for pr = 1:length(rf_table)
        rf = rf_table(pr);
        MRI_rwt_MASTeR_function(rf,opts);
        % MRI_rwt_MASTeR_function;
    end
    
    % eval(sprintf(['save ', sim_dir,'/',mname,'_',Filename,'_Sampling-',SM,'_SD%d_rtype-',REG_TYPE,'_mtype-',MOTION_TYPE,'_SOLVER-%s.mat SER_stack'],SD,SOLVER))
    diary off
end
if matlabpool('size');
    matlabpool close 
end