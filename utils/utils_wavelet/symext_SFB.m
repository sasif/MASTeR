% Symmetric Extension Synthesis Filter Bank
function y = symext_SFB(x,h0t, h1t,sym,shift,sym_12_21)

% x: input signal (first half: scaling coeffs., 2nd half: wavelet coeffs.)
% h0t: lowpass filter coeffs.
% h1t: highpass filter coeffs.
% sym:     1 means (1,1) symmetry on analysis side
%          2 means (2,2) symmetry on synthesis side (doesn't work with
%          shifted filter)
%
% sym_12_21: tells which symmetric extension to use on the synthesis side.
%            1 scaling coeffs. Es(1,2), wavelet coeffs. Es(2,1) 
%            0 scaling coeffs. Es(2,1), wavelet coeffs. Es(1,2) 
%
% if shift is odd in scaling filter @ analysis side, extension for scaling coeffs. will be
% Es(2,1) and for wavelet coeffs. will be Es(1,2)
%
% shift: shift of analysis scaling filter center from zero.
% y: output of a Synthesis bank

% This function symmetrically extend scaling and wavelet coeffs., upsample by 2 and pass that through
% specified filter.
%
% shift takes care of the symmetry when filter is not centered at the zero
%
% On Synthesis side if scaling filter is advanced by 'p' samples, wavelet
% filter will be delayed by 'p' samples, and converse will be true on the
% analysis side. So if scaling filter @ analysis side is delayed by 1,
% scaling filter on synthesis side will be advanced by 1. wavelet @
% analysis will be advanced by 1 and wavelet @ synthesis will be delayed by
% 1.

% Written by Salman Asif - Georgia Tech
% Created: February 2008

l = length(h0t);
L = l/2-1;
ln = l-L-1;
n = length(x);
ys = zeros(1,n);
yw = zeros(1,n);
y = zeros(1,n);

shift_scale = -shift;
shift_wave = shift;

if sym == 1
    if sym_12_21 == 1
        sym_scale = 12;
        sym_wave = 21;
    else
        sym_scale = 21;
        sym_wave = 12;
    end

    %% Lowpass filter

    %% Symmetric Extension (1,2)
    if sym_scale == 12
        for m = 0:n-1
            taps = [];
            indces = [];
            for k = 0:l-1
                % kx = k-L+1+m;
                kx = k+m-ln+2*shift_scale;
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
                    ys(m+1) = ys(m+1)+x(kxp+1)*h0t(l-k);
                end
            end
            %         m
            %         taps
            %         indces
        end
    end

    %% Symmetric extension (2,1)
    if sym_scale == 21
        for m = 0:n-1
            taps = [];
            indces = [];
            for k = 0:l-1
                % kx = k-L+1+m;
                kx = k+m-ln+2*shift_scale;
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
                    ys(m+1) = ys(m+1)+x(kxp+1)*h0t(l-k);
                end
            end
            %         m
            %         taps
            %         indces
        end
    end

    %% High Pass filter
    %% Symmetric Extension (1,2)
    if sym_wave == 12
        for m = 0:n-1
            taps = [];
            indces = [];
            for k = 0:l-1
                % kx = k-L+1+m;
                kx = k+m-ln+2*shift_wave;
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
                    yw(m+1) = yw(m+1)+x(kxp+1+n/2)*h1t(l-k);
                end
            end
            %         m
            %         taps
            %         indces
        end
    end

    %% Symmetric extension (2,1)
    if sym_wave == 21
        for m = 0:n-1
            taps = [];
            indces = [];
            for k = 0:l-1
                % kx = k-L+1+m;
                kx = k+m-ln+2*shift_wave;
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
                    yw(m+1) = yw(m+1)+x(kxp+1+n/2)*h1t(l-k);
                end
            end
        end
    end
    y = ys+yw;
end

%% Even Filter: (2,2) even symmetric ext. for lowpass, (2,2) odd symmetric ext. for highpasss

if sym == 2
    for m = 0:n-1
        taps = [];
        indces = [];
        wave_asym_flag = 0;
        for k = 0:l-1
            % kx = k-L+1+m;
            kx = k+m-ln;
            if mod(kx,2) == 0
                kx = kx/2;
                if(kx<0)
                    a = mod((-kx-1),n);
                    if a > n/2-1
                        kxp = n-a-1;
                        wave_asym_flag = 0;
                    else
                        kxp = a;
                        wave_asym_flag = 1;
                    end
                else
                    if (kx>n/2-1)
                        a = mod(kx,n);
                        if a > n/2-1
                            kxp = n-a-1;
                            wave_asym_flag = 1;
                        else
                            kxp = a;
                            wave_asym_flag = 0;
                        end
                    else
                        kxp = kx;
                        wave_asym_flag = 0;
                    end
                end
                taps = [taps l-k-1];
                indces = [indces (-1)^wave_asym_flag*kxp];
                y(m+1) = y(m+1)+(x(kxp+1)*h0t(l-k)+x(kxp+1+n/2)*h1t(l-k)*(-1)^wave_asym_flag);
            end
        end
%         m
%         taps
%         indces
    end
end