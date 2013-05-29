function [f01, f02] = motionEstCWT_mod(w0, w1, opts,  J, Faf, Fsf, af, sf, ROW, COL);
% w0 - reference frame CWT coeffs.
% w1 - motion frame CWT coeffs.
% opts - options structure

% Combine the surfaces after interpolation!!!

if isfield(opts,'stp_scale')
    stp_scale = opts.stp_scale;
else
    stp_scale = max(J-2,2); % use large stp_scale for large images
end

if isfield(opts,'iker');
    iker = opts.iker;
else
    iker = 1; % interpolation kernel
end

if isfield(opts,'plots');
    plots_on = opts.plots;
else
    plots_on = 0;
end
if isfield(opts,'SD_combine');
    SD_combine = opts.SD_combine;
    wt_SD = opts.wt_SD;
else
    SD_combine = 0;
end
if isfield(opts,'motion_wrap');
    motion_wrap = opts.motion_wrap;
else
    motion_wrap = 0;
end

if isfield(opts,'itype');
    itype = opts.itype;
else
    itype = 'linear'; % interpolation kernel
end

if isfield(opts,'itype_dir');
    left_right_interp = opts.itype_dir;
else
    left_right_interp = 1; % 0 - center left, 1 - center right interpolation
end

Omega_nm = opts.Omega_nm;

% CWT motion estimation
% w0 = cplxdual2D(I0, J, Faf, af);
% w1 = cplxdual2D(I1, J, Faf, af);

A_up = zeros(size(w0{J}{1}{1}{1}));
B_up = A_up;
C_up = A_up;
D_up = A_up;
E_up = A_up;
G_up = A_up;

for scale = J:-1:stp_scale
    if SD_combine == 0
        A_m = zeros(size(w0{scale}{1}{1}{1}));
        B_m = A_m;
        C_m = A_m;
        D_m = A_m;
        E_m = A_m;
        G_m = A_m;
    else
        A_m = A_up/wt_SD;
        B_m = B_up/wt_SD;
        C_m = C_up/wt_SD;
        D_m = D_up/wt_SD;
        E_m = E_up/wt_SD;
        G_m = G_up/wt_SD;
    end
    for d1 = 1:2
        for d2 = 1:3
            % real part
            w_real = w0{scale}{1}{d1}{d2};
            % imaginary part
            w_imag = w0{scale}{2}{d1}{d2};
            w0_nm = w_real+1i*w_imag;
            
            % real part
            w_real = w1{scale}{1}{d1}{d2};
            % imaginary part
            w_imag = w1{scale}{2}{d1}{d2};
            w1_nm = w_real+1i*w_imag;
            
            if scale < J
                % Interpolation can be improved by using different
                % kernels
                % Here we are warping the coefficient at this scale by the
                % previous motion estimate
                switch iker
                    case 1
                        % Staircase kernel
                        w0_nm_old = w0_nm;
                        w0_nm = w0_nm.*exp(1i*(Omega_nm{scale}{d1}{d2}{1}.*f01_up+Omega_nm{scale}{d1}{d2}{2}.*f02_up));
                        
                    case 2
                        % Bilinear kernel
                        w0_nm_old = w0_nm;
                        % horizontal
                        f01_up_pos = f01_up>0;
                        f01_up_neg = f01_up<=0;
                        w0_nm_hp_0 = (1-f01_up/2^scale).*(w0_nm_old.*exp(1i*(Omega_nm{scale}{d1}{d2}{1}.*f01_up)));
                        w0_nm_hp_1 = f01_up(:,1:end-1)/2^scale.*(w0_nm_old(:,2:end).*exp(1i*(Omega_nm{scale}{d1}{d2}{1}(:,2:end).*(-2^scale+f01_up(:,1:end-1)))));
                        w0_nm_hp = w0_nm_hp_0;
                        w0_nm_hp(:,1:end-1) = w0_nm_hp_0(:,1:end-1)+w0_nm_hp_1;
                        w0_nm_hp(f01_up_neg) = 0;
                        
                        f01_up_neg = f01_up<=0;
                        w0_nm_hn_0 = (1+f01_up/2^scale).*(w0_nm_old.*exp(1i*(Omega_nm{scale}{d1}{d2}{1}.*f01_up)));
                        w0_nm_hn_1 = (-f01_up(:,1:end-1)/2^scale).*(w0_nm_old(:,2:end).*exp(1i*(Omega_nm{scale}{d1}{d2}{1}(:,2:end).*(2^scale+f01_up(:,1:end-1)))));
                        w0_nm_hn = w0_nm_hn_0;
                        w0_nm_hn(:,1:end-1) = w0_nm_hn_0(:,1:end-1)+w0_nm_hn_1;
                        w0_nm_hn(f01_up_pos) = 0;
                        
                        w0_nm_h = w0_nm_hp;
                        w0_nm_h(f01_up_neg) = w0_nm_hn(f01_up_neg);
                        
                        % vertical
                        f02_up_pos = f02_up>0;
                        f02_up_neg = f02_up<=0;
                        w0_nm_vp_0 = (1-f02_up/2^scale).*(w0_nm_h.*exp(1i*(Omega_nm{scale}{d1}{d2}{2}.*f02_up)));
                        w0_nm_vp_1 = f02_up(:,1:end-1)/2^scale.*(w0_nm_h(:,2:end).*exp(1i*(Omega_nm{scale}{d1}{d2}{2}(:,2:end).*(-2^scale+f02_up(:,1:end-1)))));
                        w0_nm_vp = w0_nm_vp_0;
                        w0_nm_vp(:,1:end-1) = w0_nm_vp_0(:,1:end-1)+w0_nm_vp_1;
                        w0_nm_vp(f02_up_neg) = 0;
                        
                        w0_nm_vn_0 = (1+f02_up/2^scale).*(w0_nm_h.*exp(1i*(Omega_nm{scale}{d1}{d2}{2}.*f02_up)));
                        w0_nm_vn_1 = (-f02_up(:,1:end-1)/2^scale).*(w0_nm_h(:,2:end).*exp(1i*(Omega_nm{scale}{d1}{d2}{2}(:,2:end).*(2^scale+f02_up(:,1:end-1)))));
                        w0_nm_vn = w0_nm_vn_0;
                        w0_nm_vn(:,1:end-1) = w0_nm_vn_0(:,1:end-1)+w0_nm_vn_1;
                        w0_nm_vn(f02_up_pos) = 0;
                        
                        w0_nm_v = w0_nm_vp;
                        w0_nm_v(f02_up_neg) = w0_nm_vn(f02_up_neg);
                        w0_nm = w0_nm_v;
                    otherwise
                        disp('NOA');
                end
            end
            
            theta_nm = angle((w1_nm+eps)./(w0_nm+eps));
            theta_nm = theta_nm.*((abs(w1_nm)>=1e-6) & (abs(w0_nm)>=1e-6));
            Omega_nm_temp1 = Omega_nm{scale}{d1}{d2}{1};
            Omega_nm_temp2 = Omega_nm{scale}{d1}{d2}{2};
            
            Omega_nm_temp1 = Omega_nm_temp1.*(abs(Omega_nm_temp1)>=1e-6);
            Omega_nm_temp2 = Omega_nm_temp2.*(abs(Omega_nm_temp2)>=1e-6);
            
            % Surfaces are computed using formula
            % SD^{n,m}(f) = |D1D2|(theta(f) - theta^{n,m}(n))^2
            % = |D1D2|(2^m Omega^{m,n}f - theta^{n,m}(n))^2,
            % where theta^{n,m}(n) = angle(D2/D1)
            % D2 <--> w1, D1 <--> w0.
            %
            D1D2 = abs(w1_nm.*w0_nm);
            A_nm = D1D2.*(Omega_nm_temp1).^2;
            B_nm = D1D2.*(Omega_nm_temp2).^2;
            C_nm = D1D2.*(Omega_nm_temp1).*(Omega_nm_temp2)*2;
            D_nm = D1D2.*(-2*theta_nm.*Omega_nm_temp1);
            E_nm = D1D2.*(-2*theta_nm.*Omega_nm_temp2);
            G_nm = theta_nm.^2;
            
            A_m = A_m+A_nm;
            B_m = B_m+B_nm;
            C_m = C_m+C_nm;
            D_m = D_m+D_nm;
            E_m = E_m+E_nm;
            G_m = G_m+G_nm;
        end
    end
    den = C_m.^2-4*A_m.*B_m;
    f01 = (2*B_m.*D_m-C_m.*E_m)./(den+eps);
    f02 = (2*A_m.*E_m-C_m.*D_m)./(den+eps);
    f01 = f01.*(abs(f01)>=1e-6);
    f02 = f02.*(abs(f02)>=1e-6);
    
    f01_t = f01;
    f02_t = f02;
    
    if motion_wrap == 1
        % NOT SURE ABOUT THIS "WRAPING"
        ind_over_f01 = (abs(f01)/2^scale>.5);
        ind_over_f02 = (abs(f02)/2^scale>.5);
        f01(ind_over_f01) = (mod((f01(ind_over_f01)),2^scale/2)).*sign(f01(ind_over_f01));
        f02(ind_over_f02) = (mod((f02(ind_over_f02)),2^scale/2)).*sign(f02(ind_over_f02));
    end
    %         if iker == 2
    %             % Proper way would be to project them to the feasible set.
    %             ind_over_f01 = (abs(f01)/2^scale>.5);
    %             ind_over_f02 = (abs(f02)/2^scale>.5);
    %             f01(ind_over_f01) = (mod(abs(f01(ind_over_f01)),2^scale/2)).*sign(f01(ind_over_f01));
    %             f02(ind_over_f02) = (mod(abs(f02(ind_over_f02)),2^scale/2)).*sign(f02(ind_over_f02));
    %         end
    if scale < J
        f01 = f01+f01_up;
        f02 = f02+f02_up;
        
        ind_over_f01 = (abs(f01)/2^J>.5);
        ind_over_f02 = (abs(f02)/2^J>.5);
        f01(ind_over_f01) = (mod((f01(ind_over_f01)),2^J/2)).*sign(f01(ind_over_f01));
        f02(ind_over_f02) = (mod((f02(ind_over_f02)),2^J/2)).*sign(f02(ind_over_f02));
    end
    
    % w0_nm_MC = w0_nm.*exp(1i*(Omega_nm{scale}{d1}{d2}{1}.*f01+Omega_nm{scale}{d1}{d2}{2}.*f02));
    
    % Upsample motion vectors
    %     f01_row = interp1([1:ROW/2^scale]',f01, 1:.5:ROW/2^scale+.5,itype,'extrap');
    %     f01_up = interp1([1:COL/2^scale]',f01_row', 1:.5:COL/2^scale+.5,itype,'extrap')';
    %     f02_row = interp1([1:ROW/2^scale]',f02, 1:.5:ROW/2^scale+.5,itype,'extrap');
    %     f02_up = interp1([1:COL/2^scale]',f02_row',1:.5:COL/2^scale+.5,itype,'extrap')';
    
    f01_row = interp1([1:ROW/2^scale]',f01, 1:.5:ROW/2^scale,itype);
    if left_right_interp
        f01_row = [f01_row;  repmat(f01_row(end,:),1,1)];
    else
        f01_row = [repmat(f01_row(1,:),1,1); f01_row];
    end
    f01_up = interp1([1:COL/2^scale]',f01_row',1:.5:COL/2^scale,itype);
    if left_right_interp
        f01_up = [f01_up; repmat(f01_up(end,:),1,1)]';
    else
        f01_up = [repmat(f01_up(1,:),1,1); f01_up]';
    end
    f02_row = interp1([1:ROW/2^scale]',f02, 1:.5:ROW/2^scale,itype);
    if left_right_interp
        f02_row = [f02_row; repmat(f02_row(end,:),1,1)];
    else
        f02_row = [repmat(f02_row(1,:),1,1); f02_row];
    end
    f02_up = interp1([1:COL/2^scale]',f02_row',1:.5:COL/2^scale,itype);
    if left_right_interp
        f02_up = [f02_up; repmat(f02_up(end,:),1,1)]';
    else
        f02_up = [repmat(f02_up(1,:),1,1); f02_up]';
    end
    % 	f01_up = medfilt2(f01_up, [5 5],'symmetric');
    %     f02_up = medfilt2(f02_up, [5 5],'symmetric');
    
    %     [xx yy] = meshgrid(1:COL/2^scale,1:ROW/2^scale);
    %     [xi yi] = meshgrid(.5:.5:COL/2^scale,.5:.5:ROW/2^scale);
    %     f01_up2 = interp2(xx,yy,f01, xi, yi,'linear',0);
    %     f02_up = interp2(xx,yy,f02, xi, yi,'linear',0);
    
    % This part prepares to initialize next level from the current estimate of
    % motion and create a cumulative squared difference surfaces.
    if SD_combine == 1
        % upsample parameters
        
        %         A_row = interp1([1:ROW/2^scale]',A_m, 1:.5:ROW/2^scale+.5,itype,'extrap');
        %         A_up = interp1([1:COL/2^scale]',A_row', 1:.5:COL/2^scale+.5,itype,'extrap')';
        %         B_row = interp1([1:ROW/2^scale]',B_m, 1:.5:ROW/2^scale+.5,itype,'extrap');
        %         B_up = interp1([1:COL/2^scale]',B_row', 1:.5:COL/2^scale+.5,itype,'extrap')';
        %         C_row = interp1([1:ROW/2^scale]',C_m, 1:.5:ROW/2^scale+.5,itype,'extrap');
        %         C_up = interp1([1:COL/2^scale]',C_row', 1:.5:COL/2^scale+.5,itype,'extrap')';
        
        A_row = interp1([1:ROW/2^scale]',A_m, 1:.5:ROW/2^scale,itype);
        if left_right_interp
            A_row = [A_row; A_row(end,:)];
        else
            A_row = [A_row(1,:); A_row];
        end
        A_up = interp1([1:COL/2^scale]',A_row', 1:.5:COL/2^scale,itype);
        if left_right_interp
            A_up = [A_up; A_up(end,:)]';
        else
            A_up = [A_up(1,:); A_up; ]';
        end
        B_row = interp1([1:ROW/2^scale]',B_m, 1:.5:ROW/2^scale,itype);
        if left_right_interp
            B_row = [B_row; B_row(end,:)];
        else
            B_row = [B_row(1,:); B_row];
        end
        B_up = interp1([1:COL/2^scale]',B_row', 1:.5:COL/2^scale,itype);
        if left_right_interp
            B_up = [B_up; B_up(end,:)]';
        else
            B_up = [B_up(1,:); B_up]';
        end
        C_row = interp1([1:ROW/2^scale]',C_m, 1:.5:ROW/2^scale,itype);
        if left_right_interp
            C_row = [C_row; C_row(end,:)];
        else
            C_row = [C_row(1,:); C_row];
        end
        C_up = interp1([1:COL/2^scale]',C_row', 1:.5:COL/2^scale,itype);
        if left_right_interp
            C_up = [C_up; C_up(end,:)]';
        else
            C_up = [C_up(1,:); C_up]';
        end
        
        D_mod = D_m*0;
        E_mod = E_m*0;
        G_mod = G_m*0;
        for d1 = 1:2
            for d2 = 1:3
                % real part
                w_real = w0{scale}{1}{d1}{d2};
                % imaginary part
                w_imag = w0{scale}{2}{d1}{d2};
                w0_nm = w_real+1i*w_imag;
                
                % real part
                w_real = w1{scale}{1}{d1}{d2};
                % imaginary part
                w_imag = w1{scale}{2}{d1}{d2};
                w1_nm = w_real+1i*w_imag;
                theta_nm = angle((w1_nm+eps)./(w0_nm+eps));
                theta_nm = theta_nm.*((abs(w1_nm)>=1e-6) & (abs(w0_nm)>=1e-6));
                
                Omega_nm_temp1 = Omega_nm{scale}{d1}{d2}{1};
                Omega_nm_temp2 = Omega_nm{scale}{d1}{d2}{2};
                
                Omega_nm_temp1 = Omega_nm_temp1.*(abs(Omega_nm_temp1)>=1e-6);
                Omega_nm_temp2 = Omega_nm_temp2.*(abs(Omega_nm_temp2)>=1e-6);
                
                D1D2 = abs(w1_nm.*w0_nm);
                D_nm = D1D2.*(2*(Omega_nm_temp1.*f01+Omega_nm_temp2.*f02-theta_nm).*Omega_nm_temp1);
                E_nm = D1D2.*(2*(Omega_nm_temp1.*f01+Omega_nm_temp2.*f02-theta_nm).*Omega_nm_temp2);
                G_nm = (Omega_nm_temp1.*f01+Omega_nm_temp2.*f02-theta_nm-theta_nm).^2;
                
                D_mod = D_mod+D_nm;
                E_mod = E_mod+E_nm;
                G_mod = G_mod+G_nm;
            end
        end
        % essentially, it updates previous value of theta by adding the
        % estimated motion displacement
        %         D_row = interp1([1:ROW/2^scale]',D_mod, 1:.5:ROW/2^scale+.5,itype,'extrap');
        %         D_up = interp1([1:COL/2^scale]',D_row', 1:.5:COL/2^scale+.5,itype,'extrap')';
        %         E_row = interp1([1:ROW/2^scale]',E_mod, 1:.5:ROW/2^scale+.5,itype,'extrap');
        %         E_up = interp1([1:COL/2^scale]',E_row', 1:.5:COL/2^scale+.5,itype,'extrap')';
        %         G_row = interp1([1:ROW/2^scale]',G_mod, 1:.5:ROW/2^scale+.5,itype,'extrap');
        %         G_up = interp1([1:COL/2^scale]',G_row', 1:.5:COL/2^scale+.5,itype,'extrap')';
        
        D_row = interp1([1:ROW/2^scale]',D_mod, 1:.5:ROW/2^scale,itype);
        if left_right_interp
            D_row = [D_row; D_row(end,:)];
        else
            D_row = [D_row(1,:); D_row];
        end
        D_up = interp1([1:COL/2^scale]',D_row', 1:.5:COL/2^scale,itype);
        if left_right_interp
            D_up = [D_up; D_up(end,:)]';
        else
            D_up = [D_up(1,:); D_up]';
        end
        E_row = interp1([1:ROW/2^scale]',E_mod, 1:.5:ROW/2^scale,itype);
        if left_right_interp
            E_row = [E_row; E_row(end,:)];
        else
            E_row = [E_row(1,:); E_row];
        end
        E_up = interp1([1:COL/2^scale]',E_row', 1:.5:COL/2^scale,itype);
        if left_right_interp
            E_up = [E_up; E_up(end,:)]';
        else
            E_up = [E_up(1,:); E_up]';
        end
        G_row = interp1([1:ROW/2^scale]',G_mod, 1:.5:ROW/2^scale,itype);
        if left_right_interp
            G_row = [G_row; G_row(end,:)];
        else
            G_row = [G_row(1,:); G_row];
        end
        G_up = interp1([1:COL/2^scale]',G_row', 1:.5:COL/2^scale,itype);
        if left_right_interp
            G_up = [G_up; G_up(end,:)]';
        else
            G_up = [G_up(1,:); G_up]';
        end
        
        den_up = C_up.^2-4*A_up.*B_up;
        f01_res = (2*B_up.*D_up-C_up.*E_up)./den_up;
        f02_res = (2*A_up.*E_up-C_up.*D_up)./den_up;
        
    end
    
    f01_old = f01_up;
    f02_old = f02_up;
    
    %% Plot intermediates
    if plots_on && scale == J-1
        W0 = CWTstruct2mat(w0,J);
        W1 = CWTstruct2mat(w1,J);
        mask = zeros(ROW, COL);
        mask(1:2^(log2(ROW)-scale),1:2^(log2(COL)-scale)) = 1;
        mask_all = [mask mask mask mask];
        W0_mod = W0.*mask_all;
        W1_mod = W1.*mask_all;
        w0_mod = CWTmat2struct(W0_mod,J);
        w1_mod = CWTmat2struct(W1_mod,J);
        
        % I0_mod = icplxdual2D_mod(w0_mod,J,Fsf,sf);
        % I1_mod = icplxdual2D_mod(w1_mod,J,Fsf,sf);
        
        symF = 1; symh = 1; symg = 2;
        wsym_temp = cw2separate(w0_mod,J);
        wsym_temp = CWTstruct2mat(wsym_temp,J);
        I0_mod = iCWT2_separate(wsym_temp, Fsf,sf, J,symF,symh,symg);
        
        wsym_temp = cw2separate(w1_mod,J);
        wsym_temp = CWTstruct2mat(wsym_temp,J);
        I1_mod = iCWT2_separate(wsym_temp, Fsf,sf, J,symF,symh,symg);
        
        %         figure(1); clf;
        %         w0_temp = cplxdual2D(I0, scale, Faf, af);
        %         w1_temp = cplxdual2D(I1, scale, Faf, af);
        %         subplot(131); imagesc(w0_temp{scale+1}{1}{1}); colormap gray; hold on;
        %         [xx yy] = meshgrid(1:COL/2^scale, 1:ROW/2^scale);
        %         quiver(xx,yy,-f01, -f02); axis ij; axis image
        %         subplot(132); imagesc(w1_temp{scale+1}{1}{1}); hold on; axis image;
        %         subplot(133); imagesc(interp2(xx,yy,w0_temp{scale+1}{1}{1},xx+f01/2^scale, yy+f02/2^scale)); axis image
        
        [xx yy] = meshgrid(1:2^(scale-1):COL, 1:2^(scale-1):ROW);
        figure(1003); clf; imagesc(I0_mod); axis image; colormap gray; hold on;
        quiver(xx,yy,-f01_up,-f02_up); axis ij; axis image
        figure(1004); clf; imagesc(I1_mod); axis image; colormap gray; hold on;
        quiver(xx,yy,-f01_up,-f02_up); axis ij; axis image
        pause(1/60);
    end
end
%
% figure(1002);
% [xx yy] = meshgrid(1:COL/2^scale, 1:ROW/2^scale);
% quiver(xx,yy,-f01, -f02); axis ij; axis image
%
