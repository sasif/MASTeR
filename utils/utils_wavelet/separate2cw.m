function w = separate2cw(w, J)
% 2nd part of cplxdual2D after cplxdual2D_sep
%
% Input:
% w{j}{m}{n} 
% {m}{n} in the above set of wavelets correspond to {col}{row}
% where 1 - regular wavelet, 2 - Hilbert transformed wavelet.
%
% Output:
% w{j}{i}{d1}{d2} - wavelet coefficients
%       j = 1..J (scale)
%       i = 1 (real part); i = 2 (imag part)
%       d1 = 1,2; (angle) 1 - positive, 2 - negative.
%       d2 = 1,2,3 (orientations, lo_hi (vertical), hi_lo (horizontal),
%       hi_hi (diagonal) respectively in the order of column_row (first filter on columns second on rows))

for j = 1:J+1
    if j == J+1
        [wp_r wm_r] = pm(w{j}{1}{1},w{j}{2}{2});
        [wp_i wm_i] = pm(w{j}{2}{1},w{j}{1}{2});
        w{j}{1}{1} = wm_r;
        w{j}{1}{2} = wp_r;
        w{j}{2}{1} = -wp_i;
        w{j}{2}{2} = -wm_i;
    else
        for m = 1:3
            %         [w{j}{1}{1}{m} w{j}{2}{2}{m}] = pm(w{j}{1}{1}{m},w{j}{2}{2}{m});
            %         [w{j}{1}{2}{m} w{j}{2}{1}{m}] = pm(w{j}{1}{2}{m},w{j}{2}{1}{m});
            
            % I think it should have been like this to ensure that we use
            % upper half of the Fourier plane (+ve vertical frequency).
            % The way low pass filter (phi_h + j phi_g) in dualfilt1 gives negative
            % frequency band. The following changes use (phi_h - j phi_g)
            % instead.
            
            % I think it should have been like this to ensure that we use
            % upper half of the Fourier plane (+ve vertical frequency).
            % The way low pass filter (phi_h + j phi_g) in dualfilt1 gives negative
            % frequency band. The following changes use (phi_h - j phi_g)
            % instead.
            switch m
                case 1
                    [wp_r wm_r] = pm(w{j}{1}{1}{m},w{j}{2}{2}{m});
                    [wp_i wm_i] = pm(w{j}{2}{1}{m},w{j}{1}{2}{m});
                    
                    w{j}{1}{1}{m} = wp_r;
                    w{j}{1}{2}{m} = wm_r;
                    w{j}{2}{1}{m} = -wm_i;
                    w{j}{2}{2}{m} = -wp_i;
                case 2
                    [wp_r wm_r] = pm(w{j}{1}{1}{m},w{j}{2}{2}{m});
                    [wp_i wm_i] = pm(w{j}{2}{1}{m},w{j}{1}{2}{m});
                    
                    w{j}{1}{1}{m} = wp_r;
                    w{j}{1}{2}{m} = wm_r;
                    w{j}{2}{1}{m} = wm_i;
                    w{j}{2}{2}{m} = wp_i;
                case 3
                    [wp_r wm_r] = pm(w{j}{1}{1}{m},w{j}{2}{2}{m});
                    [wp_i wm_i] = pm(w{j}{2}{1}{m},w{j}{1}{2}{m});
                    
                    w{j}{1}{1}{m} = wm_r;
                    w{j}{1}{2}{m} = wp_r;
                    w{j}{2}{1}{m} = wp_i;
                    w{j}{2}{2}{m} = wm_i;
                otherwise
                    disp('na re na');
            end
        end
    end
end