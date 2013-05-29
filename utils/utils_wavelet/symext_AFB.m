% Symmetric Extension Analysis Filter Bank
function y = symext_AFB(x, nj, h0, h1,sym, delay_h0)

% x: input signal
% h0: lowpass filter coeffs.
% h1; highpass filter coeffs.
% sym: type of symmetric extension 1 or 2 for (1,1) or (2,2).
% delay_h0: shift of analysis scaling filter center from zero. (only working for (1,1) symmetric ext. at this time)
% if delay_h0 for scaling filter is 1, shift for wavelet filter will be -1.
% y: output of a single subband 1st half will be scaling coeffs. and 2nd
% half wavelet coeffs.

% This function returns scaling and wavelet coeffs. after passing the 
% specified symmetric extension through the filter bank. 
%
% For odd biorthogonal filter sym = 1 (1,1) sym. extension
% For even biorthogonal filter sym = 2 (2,2) sym. extension
% 
% delay_h0 - takes care of the symmetry when filter is not centered at the zero
%
% On analysis side if scaling filter is delayed by 'p' samples, wavelet 
% filter will be advanced by 'p' samples, and converse will hold on the 
% synthesis side.

% Written by Salman Asif - Georgia Tech
% Created: February 2008

l = length(h0);
L = l/2;
ln = l-L-1;
n = length(x);
n = nj;
y = zeros(1,n);

if sym == 1
    for m = 0:n/2-1
        for k = 0:l-1

            % Scaling filter
            kx_scale = 2*m - ln + k + 2*delay_h0;
            if(kx_scale<0)
                a = mod((-kx_scale),2*(n-1));
                if a > n-1
                    kxp_scale = 2*(n-1)-a;
                else
                    kxp_scale = a;
                end

            else
                if (kx_scale>n-1)
                    a = mod(kx_scale,2*(n-1));
                    if a > n-1
                        kxp_scale = 2*(n-1)-a;
                    else
                        kxp_scale = a;
                    end

                else
                    kxp_scale = kx_scale;

                end
            end

            % Wavelet Filter
            kx_wave = 2*m - ln + k - 2*delay_h0;
            if(kx_wave<0)
                a = mod((-kx_wave),2*(n-1));
                if a > n-1
                    kxp_wave = 2*(n-1)-a;
                else
                    kxp_wave = a;
                end

            else
                if (kx_wave>n-1)
                    a = mod(kx_wave,2*(n-1));
                    if a > n-1
                        kxp_wave = 2*(n-1)-a;
                    else
                        kxp_wave = a;
                    end

                else
                    kxp_wave = kx_wave;

                end
            end
            y(m+1) = y(m+1)+x(kxp_scale+1)*h0(l-k);
            y(m+1+n/2) = y(m+1+n/2)+x(kxp_wave+1)*h1(l-k);
        end
    end
end

if sym == 2
    for m = 0:n/2-1
        taps = [];
        indces = [];
        for k = 0:l-1
            % kx = k-L+1+m;
            kx = k+2*m-ln;
            if(kx<0)
                a = mod((-kx-1),2*n);
                if a > n-1
                    kxp = 2*n-a-1;
                else
                    kxp = a;
                end

            else
                if (kx>n-1)
                    a = mod(kx,2*n);
                    if a > n-1
                        kxp = 2*n-a-1;
                    else
                        kxp = a;
                    end

                else
                    kxp = kx;

                end
            end
            y(m+1) = y(m+1)+x(kxp+1)*h0(l-k);
            y(m+1+n/2) = y(m+1+n/2)+x(kxp+1)*h1(l-k);
        end
    end
end