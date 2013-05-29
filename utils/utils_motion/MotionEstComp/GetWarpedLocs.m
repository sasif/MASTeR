function [hor_locs ver_locs] = GetWarpedLocs(f01,f02, opts, J, ROW, COL);
%
% f01 - horizontal motion
% f02 - vertical motion
%

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
left_right_interp = 1;

% if isfield(opts,'stp_scale')
%     stp_scale = opts.stp_scale;
% else
%     stp_scale = []; % use large stp_scale for large images
%
% end
% scale = stp_scale;

interp_ratio = ROW/size(f01,1);

%% Compensation/warping

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
    hor_locs = xx1;
    ver_locs = yy1;
    
    %         Ic = zeros(ROW,COL);
    %         for rr = 1:ROW
    %             for cc = 1:COL
    %                 Ic(rr,cc) = I0(yy1(rr,cc),xx1(rr,cc));
    %             end
    %         end
    
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
    f01_row = interp1([1:ROW/interp_ratio]',f01, 1:1/interp_ratio:ROW/interp_ratio,itype);
    if left_right_interp
        f01_row = [f01_row; repmat(f01_row(end,:),interp_ratio-1,1)];
    else
        f01_row = [repmat(f01_row(1,:),interp_ratio-1,1); f01_row];
    end
    f01_up = interp1([1:COL/interp_ratio]',f01_row',1:1/interp_ratio:COL/interp_ratio,itype);
    if left_right_interp
        f01_up = [f01_up; repmat(f01_up(end,:),interp_ratio-1,1)]';
    else
        f01_up = [repmat(f01_up(1,:),interp_ratio-1,1); f01_up]';
    end
    
    f02_row = interp1([1:ROW/interp_ratio]',f02, 1:1/interp_ratio:ROW/interp_ratio,itype);
    if left_right_interp
        f02_row = [f02_row; repmat(f02_row(end,:),interp_ratio-1,1)];
    else
        f02_row = [repmat(f02_row(1,:),interp_ratio-1,1); f02_row];
    end
    f02_up = interp1([1:COL/interp_ratio]',f02_row',1:1/interp_ratio:COL/interp_ratio,itype);
    if left_right_interp
        f02_up = [f02_up; repmat(f02_up(end,:),interp_ratio-1,1)]';
    else
        f02_up = [repmat(f02_up(1,:),interp_ratio-1,1); f02_up]';
    end
    
    
    int_locs = xx+f01_up;
    int_locs_neg = int_locs<1;
    int_locs_out = int_locs>COL;
    int_locs(int_locs_neg) = 1;
    int_locs(int_locs_out) = COL;
    hor_locs = int_locs;
    
    int_locs = yy+f02_up;
    int_locs_neg = int_locs<1;
    int_locs_out = int_locs>ROW;
    int_locs(int_locs_neg) = 1;
    int_locs(int_locs_out) = ROW;
    ver_locs = int_locs;
end
