function Ic = motionCompCWT(I0,f01,f02, opts);
%
% f01 - horizontal motion 
% f02 - vertical motion
%
global OMEGA N M J Faf Fsf af sf jj L QMF ROW COL

if isfield(opts,'itype');
    itype = opts.itype;
else
    itype = 'linear'; % interpolation kernel
end

if isfield(opts,'Ic_type');
    Ic_type = opts.Ic_type;
else
    Ic_type = 2; % 1 - full interpolation kernel, 2 - block replacement
end
if isfield(opts,'stp_scale')
    stp_scale = opts.stp_scale;
else
    stp_scale = max(J-2,2); % use large stp_scale for large images
end
scale = stp_scale;

%% Compensation/warping
if Ic_type == 1
   
    [xx yy] = meshgrid(1:COL,1:ROW);
    
    if strmatch(itype, 'box');
        % Another simpler interp.
        % This code does block replacement
        mbSize = ROW/size(f01,1);
        f02_up = kron(round(f02), ones(mbSize));
        f01_up = kron(round(f01), ones(mbSize));
        xx1 = xx+round(f01_up); yy1 = yy+round(f02_up);
        xx1_neg = xx1<1;
        xx1(xx1_neg) = 1;
        xx1_large = xx1>COL;
        xx1(xx1_large) = COL;
        yy1_neg = yy1<1;
        yy1(yy1_neg) = 1;
        yy1_large = yy1>COL;
        yy1(yy1_large) = COL;
        Ic = zeros(ROW,COL);
        for rr = 1:ROW
            for cc = 1:COL
                Ic(rr,cc) = I0(yy1(rr,cc),xx1(rr,cc));
            end
        end
        % Forward operator for box interpolation (block/pixel replacement)
        % Fixi = xi(new_locs);
        % Adjoint operator
        % At every index find the bins assigned to the index in the forward
        % operator and add up all the entries at those indices.
    else
        % Interpolate
%         f01(1,:) = 0;
%         f01(:,1) = 0;
%         f01(end,:) = 0;
%         f01(:,end) = 0;
%         f02(1,:) = 0;
%         f02(:,1) = 0;
%         f02(end,:) = 0;
%         f02(:,end) = 0;
        f01_row = interp1([1:ROW/2^scale]',f01, 1:1/2^scale:ROW/2^scale,itype);
        f01_row = [f01_row; repmat(f01_row(end,:),2^scale-1,1)];
        f01_up = interp1([1:COL/2^scale]',f01_row',1:1/2^scale:COL/2^scale,itype);
        f01_up = [f01_up; repmat(f01_up(end,:),2^scale-1,1)]';
        f02_row = interp1([1:ROW/2^scale]',f02, 1:1/2^scale:ROW/2^scale,itype);
        f02_row = [f02_row; repmat(f02_row(end,:),2^scale-1,1)];
        f02_up = interp1([1:COL/2^scale]',f02_row',1:1/2^scale:COL/2^scale,itype);
        f02_up = [f02_up; repmat(f02_up(end,:),2^scale-1,1)]';
         
        %     Ic2 = interp2(xx, yy, I0, xx+f01_up, yy+f02_up,'linear',0);
        Ic_temp = zeros(ROW,COL);
        for rr=1:ROW
%             Ic_temp(rr,:) = interp1([1:ROW], I0(rr,:), [1:ROW]+f01_up(rr,:), itype, 'extrap');
            int_locs = [1:COL]+f01_up(rr,:);
            int_locs_neg = int_locs<1;
            int_locs_out = int_locs>COL;
            int_locs(int_locs_neg) = 1;
            int_locs(int_locs_out) = COL;
            Ic_temp(rr,:) = interp1([1:COL], I0(rr,:), int_locs, itype);
            % Ic_temp(rr,COL-2^scale+2:COL) = I0(rr,COL-2^scale+2:COL);
        end
        Ic = zeros(ROW,COL);
        for cc = 1:COL
            %             Ic(:,cc) = interp1([1:COL]', Ic_temp(:,cc), [1:COL]'+f02_up(:,cc), itype, 'extrap');
            int_locs = [1:ROW]'+f02_up(:,cc);
            int_locs_neg = int_locs<1;
            int_locs_out = int_locs>ROW;
            int_locs(int_locs_neg) = 1;
            int_locs(int_locs_out) = ROW;
            Ic(:,cc) = interp1([1:ROW]', Ic_temp(:,cc), int_locs, itype);
            % Ic(ROW-2^scale+2:ROW,cc) = I0(ROW-2^scale+2:ROW,cc);
        end
    end
    
    stp = 1;
    %     figure(3);
    %     subplot(121); imagesc(I1); colormap gray; axis image;
    %     subplot(122); imagesc(Ic_i); colormap gray; axis image;f
    %     figure(2);
    %     quiver(xx,yy,f01_up, f02_up); axis ij; axis image
    
else
%     %% Motion compensate, replace 8x8 blocks according to rounded motion vectors
%     des_scale = 3;
%     scale_orig = scale;
%     if scale_orig > des_scale
%         scale_orig = scale;
%         scale = des_scale;
%         sdiff = scale_orig-scale;
%         f01_row = interp1([1:ROW/2^scale_orig]',f01, 1:1/2^sdiff:ROW/2^scale_orig+1-1/2^sdiff,itype,'extrap');
%         f01_up = interp1([1:COL/2^scale_orig]',f01_row',1:1/2^sdiff:COL/2^scale_orig+1-1/2^sdiff,itype,'extrap')';
%         f02_row = interp1([1:ROW/2^scale_orig]',f02, 1:1/2^sdiff:ROW/2^scale_orig+1-1/2^sdiff,itype,'extrap');
%         f02_up = interp1([1:COL/2^scale_orig]',f02_row', 1:1/2^sdiff:COL/2^scale_orig+1-1/2^sdiff,itype,'extrap')';
%         f01_r = round(f01_up);
%         f02_r = round(f02_up);
%     else
%         f01_r = round(f01);
%         f02_r = round(f02);
%     end
%     
%     % f01_temp = f01_r;
%     % f01_r(:,1) = f01_r(:,1).*(f01_r(:,1)>0);
%     % f01_r(:,end) = f01_r(:,end).*(f01_r(:,end)<0);
%     % f02_r = f02_r;
%     % f02_r(1,:) = f02_r(1,:).*(f02_r(1,:)>0);
%     % f02_r(end,:) = f02_r(end,:).*(f02_r(end,:)<0);
%     
%     Ic = 0*I0;
%     mbSize = 2^scale;
%     mbCount = 0;
%     for i = 1 : mbSize : ROW-mbSize+1
%         for j = 1 : mbSize : COL-mbSize+1
%             mbCount = mbCount+1;
%             mv_h = f01_r(ceil(i/mbSize),ceil(j/mbSize));
%             if j+mv_h < 1
%                 mv_h = 1-j;
%             elseif j+mv_h+mbSize-1 > COL
%                 mv_h = COL+1-mbSize-j;
%             end
%             mv_v = f02_r(ceil(i/mbSize),ceil(j/mbSize));
%             if i+mv_v < 1
%                 mv_v = 1-i;
%             elseif i+mv_v+mbSize-1 > ROW
%                 mv_v = ROW+1-mbSize-i;
%             end
%             Ic(i:i+mbSize-1,j:j+mbSize-1) = I0(i+mv_v:i+mv_v+mbSize-1, j+mv_h:j+mv_h+mbSize-1);
%         end
%     end
end