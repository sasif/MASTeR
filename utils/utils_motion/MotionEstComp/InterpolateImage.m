function Ic = InterpolateImage(I0,hor_locs, ver_locs, opts);
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

if isfield(opts,'stp_scale')
    stp_scale = opts.stp_scale;
else
    stp_scale = max(J-2,2); % use large stp_scale for large images
end
scale = stp_scale;

%% Compensation/warping
[xx yy] = meshgrid(1:COL,1:ROW);
if strmatch(itype, 'box');
    % Another simpler interp.
    % This code does block replacement
    xx1 = hor_locs;
    yy1 = ver_locs;
    I0_vec = I0(:);
    prev_locs = yy1(:)+(xx1(:)-1)*ROW;
    Ic_vec = I0_vec(prev_locs);
    Ic = reshape(Ic_vec,ROW,COL);

    %     Ic = zeros(ROW,COL);
    %     for rr = 1:ROW
    %         for cc = 1:COL
    %             Ic(rr,cc) = I0(yy1(rr,cc),xx1(rr,cc));
    %         end
    %     end
    
    % Forward operator for box interpolation (block/pixel replacement)
    % Fixi = xi+1(new_locs);
    % Adjoint operator
    % At every index find the bins assigned to the index in the forward
    % operator and add up all the entries at those indices.
    % Fi'xi+1 --> 
    % for ii = 1:ROW*COL; 
    % xi(ii) = sum(xi+1(find(prev_locs==ii))); 
    % end
else
    % Interpolate
    
    %     Ic2 = interp2(xx, yy, I0, xx+f01_up, yy+f02_up,'linear',0);
    Ic_temp = zeros(ROW,COL);
    for rr=1:ROW
        %             Ic_temp(rr,:) = interp1([1:ROW], I0(rr,:), [1:ROW]+f01_up(rr,:), itype, 'extrap');
        int_locs = hor_locs(rr,:);
        Ic_temp(rr,:) = interp1([1:COL], I0(rr,:), int_locs, itype);
        %Ic_temp(rr,COL-2^scale+2:COL) = I0(rr,COL-2^scale+2:COL);
    end
    Ic = zeros(ROW,COL);
    for cc = 1:COL
        %             Ic(:,cc) = interp1([1:COL]', Ic_temp(:,cc), [1:COL]'+f02_up(:,cc), itype, 'extrap');
        int_locs = ver_locs(:,cc);
        Ic(:,cc) = interp1([1:ROW]', Ic_temp(:,cc), int_locs, itype);
        %Ic(ROW-2^scale+2:ROW,cc) = I0(ROW-2^scale+2:ROW,cc);
    end
end
