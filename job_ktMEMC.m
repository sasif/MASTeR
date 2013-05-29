% This code implements kt FOCUSS with ME/MC 
% 
% Some parts of the code were adopted from kt FOCUSS code once available at 
% http://bisp.kaist.ac.kr/ktFOCUSS.htm
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

setup_path;

%% Setup simulation parameters

global ROW COL T_frames C_coils OMEGA SENSITIVITY_MAPS SAMPLING_MASK 
global STATIC_MASK PERIODIC MOTION_FIELD_FORWARD MOTION_FIELD_BACKWARD

% Change to choose dataset [ 1 2 3 4]
DATASET = 2; PRESCAN = true; rf_table = [4 8];
% DATASET = 3; PRESCAN = false; rf_table = [6 10];
MAX_frames = 16;

rwt_MAXITER = 3; % maximum iteration for reweighting
cg_maxiter = 200; % maximum conjugate gradient iterations at every reweighting step
cg_opt = 1;     % 1 -- select the solution of CG that yields minimum MSE
                % 0 -- fixed threshold for the termination criterion
Iref = 4; % Select the reference frame for ME/MC residual reconstruction from the 
            % 1 -- Last frame of the original image sequence
            % 2 -- Average of reconstructed sequence
            % 3 -- Last frame in the reconstructed sequence
            % 4 -- Average of diastole phase in the reconstructed sequence
                 % DIASTOLE_FRAMES assigned in MRI_ReadData
                 
MOTION_TYPE = 'OBMC'; % {'CWT','OF','OBMC'}; % scheme for motion estimation
oracleMotion = 0;   % 0 -- estimate motion from reconstructed frames
                    % 1 -- estimate motion from original frames (oracle)
oracleWeights = 0;  % weights selection for the iterative reweighting

% sampling scheme 
SampType = {'HYBRID'}; % HYBRID: few low-frequency + randomly sampled high frequency horizontal k-space lines
                       % RAND: randomly subsample horizontal k-space lines

rseed = 0;
SAVE_RESULTS = false;

% options in a structure form
opts = [];
opts.DATASET = DATASET; 
opts.MAX_frames = MAX_frames;
opts.rwt_MAXITER = rwt_MAXITER; 
opts.Iref = Iref;
opts.MOTION_TYPE = MOTION_TYPE;
opts.oracleMotion = oracleMotion;
opts.oracleWeights = oracleWeights;
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

sim_dir = sprintf(['results/',mname,'_pre%d_Iref%d_om%d_ow%d_MOTION-%s_cg_opt%d_cg_maxiter%d'],PRESCAN,Iref,oracleMotion,oracleWeights,MOTION_TYPE,cg_opt,cg_maxiter);

opts.mname = mname;
opts.sim_dir = sim_dir;

for ps = 1:length(SampType)

    SM = char(SampType(ps));
    opts.SM = SM;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if SAVE_RESULTS
        status = mkdir(sim_dir);
        diary(sprintf([sim_dir,'/',mname,'_',Filename,'_Sampling-',SM,'_Iref%d.txt'],Iref))
    end
    
    fprintf('data-%d, pre-%d, cg_maxiter-%d, cg_opt-%d, rwt_iter-%d, Iref-%d, MOTION-%s, om-%d, ow-%d\n',DATASET,PRESCAN,cg_maxiter, cg_opt, rwt_MAXITER, Iref, MOTION_TYPE, oracleMotion,oracleWeights);

    for pr = 1:length(rf_table)
        rf = rf_table(pr);
        % kt FOCUSS with ME/MC function
        MRI_rwt_ktMEMC_function(rf,opts);
        % MRI_rwt_ktMEMC_function;
    end
    
    % eval(sprintf(['save ', sim_dir,'/',mname,'_',Filename,'_Sampling-',SM,'_SD%d_rtype-',REG_TYPE,'_mtype-',MOTION_TYPE,'_SOLVER-%s.mat SER_stack'],SD,SOLVER))
    diary off
end
if matlabpool('size');
    matlabpool close 
end