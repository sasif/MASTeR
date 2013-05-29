function y = invfilter_symext(x,h,sym,shift)

% x: input signal (downsampled)
% h: filter coeffs.
% sym:    12 means symmetric (1,2)
%         21 means symmetric (2,1)
%         22 means symmetric (2,2)
%         33 means anti-symmetric (2,2)
%
% shift: shift of center from zero. (only working for 12 or 21 symmetric
% ext. right now)
% y: output of a single subband

% This function symmetrically extend x, upsample by 2 and pass that through
% specified filter. 
%
% shift takes care of the symmetry when filter is not centered at the zero
%
% On Synthesis side if scaling filter is advanced by 'p' samples, wavelet 
% filter will be delayed by 'p' samples, and converse will be true on the 
% analysis side.

% Written by Salman Asif - Georgia Tech
% Created: February 2008

l = length(h);
L = l/2-1;
ln = l-L-1;
n = length(x)*2;
y = zeros(1,n);

%% Symmetric Extension (1,2)
if sym == 12
    for m = 0:n-1
        taps = [];
        indces = [];
        for k = 0:l-1
            % kx = k-L+1+m;
            kx = k+m-ln+2*shift;
            if mod(kx,2) == 0
                kx = kx/2;
                if(kx<0)
                    a = mod((-kx),n-1);
                    if a > n/2-1
                        kxp = n-1-a;
                    else
                        kxp = a;
                    end

                else
                    if (kx>n/2-1)
                        a = mod(kx,n-1);
                        if a > n/2-1
                            kxp = n-1-a;
                        else
                            kxp = a;
                        end

                    else
                        kxp = kx;

                    end
                end
                taps = [taps l-k-1];
                indces = [indces kxp];
                y(m+1) = y(m+1)+x(kxp+1)*h(l-k);
            end
        end
%         m
%         taps
%         indces
    end
end

%% Symmetric extension (2,1)
if sym == 21
    for m = 0:n-1
        taps = [];
        indces = [];
        for k = 0:l-1
            % kx = k-L+1+m;
            kx = k+m-ln+2*shift;
            if mod(kx,2) == 0
                kx = kx/2;
                if(kx<0)
                    a = mod((-kx-1),n-1);
                    if a > n/2-1
                        kxp = 2*(n/2-1)-a;
                    else
                        kxp = a;
                    end

                else
                    if (kx>n/2-1)
                        a = mod(kx,n-1);
                        if a > n/2-1
                            kxp = 2*(n/2-1)-a;
                        else
                            kxp = a;
                        end

                    else
                        kxp = kx;

                    end
                end
                taps = [taps l-k-1];
                indces = [indces kxp];
                y(m+1) = y(m+1)+x(kxp+1)*h(l-k);
            end
        end
%         m
%         taps
%         indces
    end
end

% %% Not set yet... :)
% if sym == 22
%     for m = 0:n-1
%         taps = [];
%         indces = [];
%         for k = 0:l-1
%             % kx = k-L+1+m;
%             kx = k+m-ln;
%             if(kx<0)
%                 a = mod((-kx-1),2*n);
%                 if a > n-1
%                     kxp = 2*n-a-1;
%                 else
%                     kxp = a;
%                 end
% 
%             else
%                 if (kx>n-1)
%                     a = mod(kx,2*n);
%                     if a > n-1
%                         kxp = 2*n-a-1;
%                     else
%                         kxp = a;
%                     end
% 
%                 else
%                     kxp = kx;
% 
%                 end
%             end
%             taps = [taps l-k-1];
%             indces = [indces kxp];
%             y(m+1) = y(m+1)+x(kxp+1)*h(l-k);
%         end
%     end
% end