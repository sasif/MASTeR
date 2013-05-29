% Overlapped Block Motion Compensation (OBMC)
% via the motion vectors computed using exhaustive search method
%
% Input
%   imgP : The image for which we want to find motion vectors
%   mbSize : Size of the macroblock
%   mv : Motion vectors 
%  
% Optional 
%   Weight_mask : Mask to do the weighted average of predictions
%       Each pixel at relative location sb in each block is assigned a
%       weighted average of the predictions with motion vectors of each
%       neighboring block. 
%
%       Example mask for 4 neighbors (left,right,top, and bottom) of 8x8 block
%       (taken from Yao Wang notes @ Polytech Brooklyn)
%
%                       1 1 1 1 1 1 1 1
%                       1 1 1 1 1 1 1 1
%                       1 1 2 2 2 2 1 1
%                       2 2 2 2 2 2 2 2
%               1 1 1 2 4 5 5 5 5 5 5 4 2 1 1 1
%               1 1 2 2 5 5 5 5 5 5 5 5 2 2 1 1
%               1 1 2 2 5 5 6 6 6 6 5 5 2 2 1 1
%               1 1 2 2 5 5 6 6 6 6 5 5 2 2 1 1
%               1 1 2 2 5 5 6 6 6 6 5 5 2 2 1 1
%               1 1 2 2 5 5 6 6 6 6 5 5 2 2 1 1
%               1 1 2 2 5 5 5 5 5 5 5 5 2 2 1 1
%               1 1 1 2 4 5 5 5 5 5 5 4 2 1 1 1
%                       2 2 2 2 2 2 2 2
%                       1 1 2 2 2 2 1 1
%                       1 1 1 1 1 1 1 1
%                       1 1 1 1 1 1 1 1
%
%           For further details see Orchard's papers.
%
% Ouput
%   imgC : Compensated image 
%
% Written by Salman Asif

function [imgOBMC] = OBMC(imgP, mbSize, mv)

[row col] = size(imgP);

imgC = 0*imgP;

mbCount = 1;
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
%                 costs(m+p+1,n+p+1) = costFuncMAD(imgP(i:i+mbSize-1,j:j+mbSize-1), ...
%                      imgI(refBlkVer:refBlkVer+mbSize-1, refBlkHor:refBlkHor+mbSize-1), mbSize);
                costs(m+p+1,n+p+1) = costFuncMSD(imgP(i:i+mbSize-1,j:j+mbSize-1), ...
                     imgI(refBlkVer:refBlkVer+mbSize-1, refBlkHor:refBlkHor+mbSize-1), mbSize);
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

cen = [ 4 5 5 5 5 5 5 4
        5 5 5 5 5 5 5 5
        5 5 6 6 6 6 5 5
        5 5 6 6 6 6 5 5
        5 5 6 6 6 6 5 5
        5 5 6 6 6 6 5 5
        5 5 5 5 5 5 5 5
        4 5 5 5 5 5 5 4 ];
tb = [ 2 2 2 2 2 2 2 2
       1 1 2 2 2 2 1 1
       1 1 1 1 1 1 1 1
       1 1 1 1 1 1 1 1
       1 1 1 1 1 1 1 1
       1 1 1 1 1 1 1 1
       1 1 2 2 2 2 1 1
       2 2 2 2 2 2 2 2 ];
   
lr = tb';
   
% With this mask we have a neighborhood of 5 blocks (left, right, top and
% bottom), and for each pixel in our reference frame we compute its
% predicted value using 5 motion vectors for each of the five blocks and
% averaging them using the weighting mask
% s_b = sum w_k Ip(s_b-v_k(s_b)),
% where w_k are the weights and v_k gives the motion displacement of the 
% kth block. 

% NEED TO COMPLETE THIS ONE 
% (FIND EACH PIXEL IN CURRENT BLOCK BY USING WEIGHTED AVERAGE OF NEIGHBORING PREDICTED PIXELS)...

imgOBMC = ???
