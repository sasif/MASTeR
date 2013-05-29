% Copyright (c) 2005, Aroh Barjatya
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:
% 
%     * Redistributions of source code must retain the above copyright 
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright 
%       notice, this list of conditions and the following disclaimer in 
%       the documentation and/or other materials provided with the distribution
%       
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
% POSSIBILITY OF SUCH DAMAGE.



% Computes motion compensated image using the given motion vectors
%
% Input
%   imgI : The reference image 
%   motionVect : The motion vectors
%   mbSize : Size of the macroblock
%
% Ouput
%   imgComp : The motion compensated image
%
% This code was originally written by by Aroh Barjatya
% and finally modified by Hong Jung for overlapped block matching algorithm

function imgComp = motionComp(imgI, motionVect, mbSize)

[row col] = size(imgI);


% we start off from the top left of the image
% we will walk in steps of mbSize
% for every marcoblock that we look at we will read the motion vector
% and put that macroblock from refernce image in the compensated image
MASK=zeros(size(imgI));
imageComp=zeros(size(imgI));
mbCount = 1;
MV=zeros([size(imgI),2]);
for i = 1:row-mbSize+1
    for j = 1:col-mbSize+1
        
        % dy is row(vertical) index
        % dx is col(horizontal) index
        % this means we are scanning in order
        
        dy = motionVect(1,mbCount);
        dx = motionVect(2,mbCount);

        refBlkVer = i + dy;
        refBlkHor = j + dx;
        imageComp(i:i+mbSize-1,j:j+mbSize-1) = imageComp(i:i+mbSize-1,j:j+mbSize-1)+imgI(refBlkVer:refBlkVer+mbSize-1, refBlkHor:refBlkHor+mbSize-1);
        MV(i:i+mbSize-1,j:j+mbSize-1,1)=MV(i:i+mbSize-1,j:j+mbSize-1,1)+dy;
        MV(i:i+mbSize-1,j:j+mbSize-1,2)=MV(i:i+mbSize-1,j:j+mbSize-1,2)+dx;
        MASK(i:i+mbSize-1,j:j+mbSize-1) = MASK(i:i+mbSize-1,j:j+mbSize-1)+1;
        mbCount = mbCount + 1;
    end
end
MASK(find(MASK==0))=1;
MV(:,:,1)=MV(:,:,1)./MASK;
MV(:,:,2)=MV(:,:,2)./MASK;
avMV=sqrt(sum(MV.^2,3));
imgComp = imageComp./MASK;