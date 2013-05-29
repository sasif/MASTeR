function Omega_nm = estimate_localFreq(I0, sh_ref, COL, ROW, J, Faf, af);


symF = 1; symh = 1; symg = 2;
    
%% horizontal-direction
I0_sh = I0;
I0_sh(:,1:COL-sh_ref) = I0(:,sh_ref+1:COL);

% w0 = cplxdual2D_mod(I0, J, Faf, af);
% w1 = cplxdual2D_mod(I0_sh, J, Faf, af);

Fsep =  CWT2_separate(I0, Faf, af, J, symF, symh, symg);
W_temp = Fsep;% [Whh Whg Wgh Wgg]; % without pm
wsym_temp = CWTmat2struct(W_temp,J);
w0 = separate2cw(wsym_temp,J);
Fsep =  CWT2_separate(I0_sh, Faf, af, J, symF, symh, symg);
W_temp = Fsep;% [Whh Whg Wgh Wgg]; % without pm
wsym_temp = CWTmat2struct(W_temp,J);
w1 = separate2cw(wsym_temp,J);

for scale = J:-1:1
    %Real part
    for d1 = 1:2
        for d2 = 1:3
            % real part
            w_real = w0{scale}{1}{d1}{d2};
            % imaginary part
            w_imag = w0{scale}{2}{d1}{d2};
            w0_scale = w_real+1i*w_imag;
            
            % real part
            w_real = w1{scale}{1}{d1}{d2};
            % imaginary part
            w_imag = w1{scale}{2}{d1}{d2};
            w1_scale = w_real+1i*w_imag;
            
            Omega_nm{scale}{d1}{d2}{1} = angle((w1_scale+eps)./(w0_scale+eps))/sh_ref;
        end
    end
end

%% vertical-direction
I0_sh = I0;
I0_sh(1:ROW-sh_ref,:) = I0(sh_ref+1:ROW,:);

% w0 = cplxdual2D_mod(I0, J, Faf, af);
% w1 = cplxdual2D_mod(I0_sh, J, Faf, af);

Fsep =  CWT2_separate(I0, Faf, af, J, symF, symh, symg);
W_temp = Fsep;% [Whh Whg Wgh Wgg]; % without pm
wsym_temp = CWTmat2struct(W_temp,J);
w0 = separate2cw(wsym_temp,J);
Fsep =  CWT2_separate(I0_sh, Faf, af, J, symF, symh, symg);
W_temp = Fsep;% [Whh Whg Wgh Wgg]; % without pm
wsym_temp = CWTmat2struct(W_temp,J);
w1 = separate2cw(wsym_temp,J);

for scale = J:-1:1
    %Real part
    for d1 = 1:2
        for d2 = 1:3
            % real part
            w_real = w0{scale}{1}{d1}{d2};
            % imaginary part
            w_imag = w0{scale}{2}{d1}{d2};
            w0_scale = w_real+1i*w_imag;
            
            % real part
            w_real = w1{scale}{1}{d1}{d2};
            % imaginary part
            w_imag = w1{scale}{2}{d1}{d2};
            w1_scale = w_real+1i*w_imag;
            
            Omega_nm{scale}{d1}{d2}{2} = angle((w1_scale+eps)./(w0_scale+eps))/sh_ref;
        end
    end
end
