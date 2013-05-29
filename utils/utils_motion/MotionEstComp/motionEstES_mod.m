% Computes motion vectors using exhaustive search method
%
% Input
%   imgP : The image for which we want to find motion vectors
%   imgI : The reference image
%   mbSize : Size of the macroblock
%   p : Search parameter  (read literature to find what this means)
%
% Ouput
%   motionVect : the motion vectors for each integral macroblock in imgP
%   EScomputations: The average number of points searched for a macroblock
%
% Written by Aroh Barjatya

function [motionVect, EScomputations, imgC] = motionEstES_mod(imgP, imgI, mbSize, p)

[row col] = size(imgI);

vectors = zeros(2,row*col/mbSize^2);
costs = ones(2*p + 1, 2*p +1) * 65537;

computations = 0;

% we start off from the top left of the image
% we will walk in steps of mbSize
% for every marcoblock that we look at we will look for
% a close match p pixels on the left, right, top and bottom of it

mbCount = 1;
imgC = 0*imgI;
for i = 1 : mbSize : row-mbSize+1
    for j = 1 : mbSize : col-mbSize+1
        
        % the exhaustive search starts here
        % we will evaluate cost for  (2p + 1) blocks vertically
        % and (2p + 1) blocks horizontaly
        % m is row(vertical) index
        % n is col(horizontal) index
        % this means we are scanning in raster order
        
        for m = -p : p        
            for n = -p : p
                refBlkVer = i + m;   % row/Vert co-ordinate for ref block
                refBlkHor = j + n;   % col/Horizontal co-ordinate
                if ( refBlkVer < 1 || refBlkVer+mbSize-1 > row ...
                        || refBlkHor < 1 || refBlkHor+mbSize-1 > col)
                    continue;
                end
                costs(m+p+1,n+p+1) = costFuncMAD(imgP(i:i+mbSize-1,j:j+mbSize-1), ...
                     imgI(refBlkVer:refBlkVer+mbSize-1, refBlkHor:refBlkHor+mbSize-1), mbSize);
%                 costs(m+p+1,n+p+1) = costFuncMSD(imgP(i:i+mbSize-1,j:j+mbSize-1), ...
%                      imgI(refBlkVer:refBlkVer+mbSize-1, refBlkHor:refBlkHor+mbSize-1), mbSize);
                computations = computations + 1;
                
            end
        end
        
        % Now we find the vector where the cost is minimum
        % and store it ... this is what will be passed back.
        
        [dx, dy, min] = minCost(costs); % finds which macroblock in imgI gave us min Cost
        vectors(1,mbCount) = dy-p-1;    % row co-ordinate for the vector
        vectors(2,mbCount) = dx-p-1;    % col co-ordinate for the vector
        
        imgC(i:i+mbSize-1,j:j+mbSize-1) = imgI(i+vectors(1,mbCount):i+vectors(1,mbCount)+mbSize-1, j+vectors(2,mbCount):j+vectors(2,mbCount)+mbSize-1);
 
        mbCount = mbCount + 1;
        costs = ones(2*p + 1, 2*p +1) * 65537;
    end
end
motionVect = vectors;
EScomputations = computations/(mbCount - 1);

itype = 'linear';
f01 = reshape(motionVect(2,:),col/mbSize,row/mbSize)'; % horizontal
f02 = reshape(motionVect(1,:),col/mbSize,row/mbSize)'; % vertical
% f01_p = abs(f01)==p;
% f01(f01_p) = p;
% f02_p = abs(f02)==p;
% f02(f02_p) = p;

[xx yy] = meshgrid(1:col,1:row);
f01_row = interp1([1:row/mbSize]',f01, 1:1/mbSize:row/mbSize+1-1/mbSize,itype,'extrap');
f01_up = interp1([1:col/mbSize]',f01_row',1:1/mbSize:col/mbSize+1-1/mbSize,itype,'extrap')';
f02_row = interp1([1:row/mbSize]',f02, 1:1/mbSize:row/mbSize+1-1/mbSize,itype,'extrap');
f02_up = interp1([1:col/mbSize]',f02_row',1:1/mbSize:col/mbSize+1-1/mbSize,itype,'extrap')';
xx1 = xx+round(f01_up); yy1 = yy+round(f02_up);
xx1_neg = xx1<1;
xx1(xx1_neg) = 1;
xx1_large = xx1>col;
xx1(xx1_large) = col;
yy1_neg = yy1<1;
yy1(yy1_neg) = 1;
yy1_large = yy1>col;
yy1(yy1_large) = col;

% This code does block replacement
f01 = reshape(motionVect(2,:),col/mbSize,row/mbSize)'; % horizontal
f02 = reshape(motionVect(1,:),col/mbSize,row/mbSize)'; % vertical
f02_up = kron(f02, ones(mbSize));
f01_up = kron(f01, ones(mbSize));
xx1 = xx+f01_up;
yy1 = yy+f02_up;

Ic = zeros(row,col);
for rr = 1:row
    for cc = 1:col
        Ic(rr,cc) = imgI(yy1(rr,cc), xx1(rr,cc));
    end
end