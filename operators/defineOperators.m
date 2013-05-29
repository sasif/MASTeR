% L2 reg operators
A_h = []; At_h = [];

% L1 reg operator
U = @(z) []; Ut = @(z) []; 

% Check if any L1 regularization is present
L1reg = 0;

% Sensing operator
% A_meas = @(z) A_OMEGA(z);
% At_meas = @(z) At_OMEGA(z);

len_A = len_y;
len_U = 0;

if SP_REG 
    U_sp = @(z) psiT(z);
    Ut_sp = @(z) psi(z);
    len_U = red*RCT;
    L1reg = 1;
    fprintf(': SP_REG :');
else
    U_sp = @(z) [];
    Ut_sp = @(z) 0;
    fprintf(': No SP_REG :');
end

if KT_REG    
    MRI_REG_OPERATOR; 
    switch KT_NORM
        case 'L2'
            A_kt = @(z) [A_meas(z); (wt_LFt)*A_kt_wt(z)];
            At_kt = @(z) [At_meas(z(1:len_A))+(wt_LFt)*At_kt_wt(z(len_A+1:end))];
            len_A = len_A + RCT;
            U_kt = @(z) U_sp(z);
            Ut_kt = @(z) Ut_sp(z);
        case 'L1'          
            U_kt = @(z) [U_sp(z); wt_LFt*A_kt_wt(z)];
            Ut_kt = @(z) [Ut_sp(z(1:len_U))+wt_LFt*At_kt_wt(z(len_U+1:end))];
            len_U = len_U + RCT;
            A_kt = @(z) A_meas(z);
            At_kt = @(z) At_meas(z);
            L1reg = 1;            
    end
    fprintf(': KT_NORM = %s, wt_KT = %3.4g :',KT_NORM, wt_LFt);
else
    U_kt = @(z) U_sp(z);
    Ut_kt = @(z) Ut_sp(z);
    A_kt = @(z) A_meas(z); 
    At_kt = @(z) At_meas(z);
    fprintf(': No KT_REG :');
end

if MC_REG    
    wt_mc = (m_iter)*wt_MC0;    
    if KEEP_REG == 0
        U_kt = @(z) [];
        Ut_kt = @(z) 0;
        L1reg = 0;
        len_U = 0;
        fprintf(': Only MC REG. :');        
    end
    len_motion = length(desired_frames)*ROW*COL;
    switch MC_NORM
        case 'L2'
            A_h = @(z) [A_kt(z); (wt_mc)*A_motion(z)];
            At_h = @(z) [At_kt(z(1:len_A))+(wt_mc)*At_motion(z(len_A+1:end))];
            len_A = len_A + m_red*len_motion;
            U_h = @(z) U_kt(z);
            Ut_h = @(z) Ut_kt(z);            
        case 'L1'
            U_h = @(z) [U_kt(z); wt_mc*A_motion(z)];
            Ut_h = @(z) [Ut_kt(z(1:len_U))+wt_mc*At_motion(z(len_U+1:end))];
            len_U = len_U + m_red*len_motion;
            A_h = @(z) A_kt(z);
            At_h = @(z) At_kt(z);
            L1reg = 1;
        case 'psiL1'
            U_h = @(z) [U_kt(z); wt_mc*psiT(A_motion(z))];
            Ut_h = @(z) [Ut_kt(z(1:len_U))+wt_mc*At_motion(psi(z(len_U+1:end)))];
            len_U = len_U + red*m_red*len_motion;
            A_h = @(z) A_kt(z);
            At_h = @(z) At_kt(z);
            L1reg = 1;
    end
    fprintf(': MC_NORM = %s, wt_MC = %d : \n',MC_NORM, wt_mc);
else
    U_h = @(z) U_kt(z);
    Ut_h = @(z) Ut_kt(z);
    A_h = @(z) A_kt(z); 
    At_h = @(z) At_kt(z);
    fprintf(': No MC_REG : \n');
end

if RWT && m_iter > 1
    W = abs(U_h(Ir_cube(:)));
    W(W<1e-4) = 1e-4;
    Wi = 1./W;
    U_h = @(z) Wi.*U_h(z); 
    Ut_h = @(z) Ut_h(Wi.*z);
    fprintf('Weighting applied to the regularizer -- min(W)=%3.4g, max(W)=%3.4g. \n',min(W), max(W));
end