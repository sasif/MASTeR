% function job_MASTeR
% (data, pre, sd, spreg, ktreg, ktnorm, keepReg, mcnorm, mcdir, per, om, rwt)
%
%-------------------------------------------+
% Author: M. Salman Asif @ Georgia Tech
% Email: sasif@gatech.edu
%-------------------------------------------+

%% Description of the parameters
% 
% data - dataset
% pre - with prescan (1) or without prescan (0)
% sd - static/dynamic partition (0/1)
% spreg - apply spatial regularization (1) or not (2)
% ktreg - apply kt reg (temporal DFT) for the initial iteration (1) or not (2)
% ktnorm - norm for KT reg. L1 (1) or L2 (2) (irrelevant if ktreg = 0)
% keepReg - use previous reg. with MC (1) or not (0)
% mcnorm - norm for MC reg L1 (1) or L2 (2)
% mcdir - type for motion compensation, Bi (1), wtBi (2), FOR (3), BACK (4)
% per - connect the first and the last frame (1) or not (0)
% om - flag to use oracle motion (1) or not (0)
% rwt - reweighting during motion adaptation steps.. ??

%% Path/machine setup
%
% numCores = str2double(getenv('NUMBER_OF_PROCESSORS'));
% numCores = feature('numCores');
% mpSize = min(numCores-1,4);
% if matlabpool('size')~=0
%     matlabpool close force;
% end
% matlabpool('open', mpSize);

setup_path

%% Setup simulation parameters

global ROW COL T_frames C_coils OMEGA SENSITIVITY_MAPS SAMPLING_MASK 
global STATIC_MASK PERIODIC MOTION_FIELD_FORWARD MOTION_FIELD_BACKWARD

SAVE_RESULTS = true;

NORM_LIST = {'L1','L2'};

vec = @(x) x(:);

% Choose an appropriate solver and regularization scheme
%
% We used the following in the paper 
spreg = 1;  keepReg = 0; ktreg = 0; 
% spatial regularization for the initialization, 
% motion-adaptive regularization only for the subsequent iterations
% no temporal DFT 
%
% Other choices:
% spreg = 1;  keepReg = 1; ktreg = 0;  % Keep spatial reg. during motion-adaptation step 
% spreg = 0;  keepReg = 0; ktreg = 1; % Use kt reg. for the initialization

% related parameters
ktnorm = 1; 
mcnorm = 1; 
mcdir = 0; per = 1; om = 0; sd = 0;
rwt = 0;

% solver selection for L1-norm minimization
% constr_solver = 'constr'; 
constr_solver = 'unconstr'; 

SOLVER = 'NESTA_L1'; 
% SOLVER = 'NESTA_TV';
% SOLVER = 'spg'; red_wave = 0;
% SOLVER = 'SALSA';
CONSTR_SOLVER = constr_solver;

% Sampling type 
SampType = {'HYBRID'}; % HYBRID: few low-frequency + randomly sampled high frequency horizontal k-space lines
                       % RAND: randomly subsample horizontal k-space lines
rseed = 0;

% Data setup
% data set; prescan? ; reduction factor
% data = 3; pre = 0; rf_table = [2:2:12]; % dataset 'wrap003'
data = 2; pre = 1; rf_table = [2:2:12]; % dataset 'P501767'

DATASET = data; %1:6
PRESCAN = pre; %0 or 1
PERIODIC = per; % 0 or 1

MAX_frames = 16;
SD = sd;

% spatial regularization 
SP_REG = spreg;     % 0 or 1 (use spatial regularization or not)
L = 4;              % wavelet level
wave_red = 1;       % 1 - use redundant wavelets (CWT) or 0 - use regular wavelets

% KT regularization 
KT_TYPE = 'SMOOTH';     % {'SMOOTH', 'BRICK'};
KT_REG = ktreg;     % 0 or 1 (use KT reg. or not)
KEEP_REG = keepReg;       % 0 or 1 (use KT reg. with MC iterations or not)
KT_NORM = NORM_LIST{ktnorm};

% ME/MC parameters
MOTION_TYPE = 'CWT';    % {'CWT','OF','OBMC','BM'};
MDIR_LIST = {'Bi','wtBi','FOR','BACK'};
MC_NORM = NORM_LIST{mcnorm};
MOTION_DIR = MDIR_LIST{mcdir}; % 1 - Bi, 2 - wtBi, 3 - forward, 4 - backward
oracleMotion = om; % Oracle motion 0 or 1

% Apply weighting on the regularizers
RWT = rwt;

fprintf('data-%d, pre-%d, sd-%d, spreg-%d, ktreg-%d, ktnorm-%d, keepReg-%d, mcnorm-%d, mdir-%d, per-%d, om-%d, rwt-%d \n',data,pre,SD,spreg,ktreg,ktnorm,keepReg,mcnorm,mcdir,per,om,rwt);

%% Identify the solver
if spreg
    fprintf('Spatial regularization\n')
else
    fprintf('No spatial regularization\n')
end
if SD
    fprintf('Static-dynamic partition applied\n');
else
    fprintf('No static-dynamic partition\n');
end
if ktreg
    fprintf('KT regularization with %s-norm',KT_NORM)
    if KEEP_REG
        fprintf(' , also use with MC iterations.\n')
    else
        fprintf(' for initialization only. \n')
    end
else
    fprintf('No KT regularization.\n')
end
if rwt 
    fprintf('Weighting applied on the regularizers \n');
end

opts = [];
opts.DATASET = DATASET; %Change to choose dataset [ 1 2 3 4]
opts.MAX_frames = MAX_frames;
opts.SD = SD; % 0 (no SD regions) or 1 (SD regions as specified in STATIC_MASK)
opts.PRESCAN = PRESCAN;

opts.SP_REG = SP_REG;   % use SP reg or not
opts.KT_REG = KT_REG;   % use KT reg or not
opts.KEEP_REG = KEEP_REG;     % use with MC or not
opts.KT_NORM = KT_NORM; % 'L1' or 'L2'
opts.KT_TYPE = KT_TYPE; % 'SMOOTH' or 'BRICK'

opts.MOTION_TYPE = MOTION_TYPE;
opts.MOTION_DIR = MOTION_DIR; % 'FOR','BACK', 'Bi', 'wtBi';
opts.MC_NORM = MC_NORM;
opts.oracleMotion = oracleMotion; 

opts.RWT = RWT;

opts.L = L; opts.wave_red = wave_red;
opts.SOLVER = SOLVER; 
opts.CONSTR_SOLVER = CONSTR_SOLVER; 

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

opts.PERIODIC = PERIODIC;
% mname = sprintf([mname,'_pre%d_sd%d_spreg%d_kt%d_%s_%d_mc%s_%s_%s_per%d_om%d'],PRESCAN,SD,SP_REG,KT_REG,KT_NORM,KEEP_REG,MC_NORM,MOTION_DIR,MOTION_TYPE,PERIODIC,oracleMotion);
sim_dir = sprintf(['results/',mname,'_pre%d_sd%d_spreg%d_kt%d_%s_%d_mc%s_%s_%s_per%d_om%d_rwt%d_%s%s'],PRESCAN,SD,SP_REG,KT_REG,KT_NORM,KEEP_REG,MC_NORM,MOTION_DIR,MOTION_TYPE,PERIODIC,oracleMotion, RWT, SOLVER, CONSTR_SOLVER);

disp([sim_dir,' ',Filename]);
opts.mname = mname;
opts.sim_dir = sim_dir;
opts.rseed = rseed;

for ps = 1:length(SampType)
    
    SM = char(SampType(ps));
    opts.SM = SM;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if SAVE_RESULTS
        status = mkdir(sim_dir);
        diary(sprintf([sim_dir,'/',mname,'_',Filename,'_SM-',SM,'.txt']))
    end
    for pr = 1:length(rf_table)
        rf = rf_table(pr);
        MASTeR_function(rf,opts);
    end
    % eval(sprintf(['save ',
    % sim_dir,'/',mname,'_',Filename,'_Sampling-',SM,'_SD%d_rtype-',REG_TYPE,'_mtype-',MOTION_TYPE,'_SOLVER-%s.mat SER_stack'],SD,SOLVER))
end
diary off
% matlabpool close