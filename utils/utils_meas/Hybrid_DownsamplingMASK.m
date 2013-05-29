function [M3 M2] = Hybrid_DownsamplingMASK(nY, nX, nframe, down)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Random sampling mask %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nY: size of phase encoding direction
% nX: size of frequency encoding direction
% nframe: size of time frames
% M: sampling mask, dimension=(nY, nX, nframe), 
%    sampled point:1, otherwise: 0;
%    SAMPLING PATTERN is RANDOM!!!
%    1<=ky<=4, nY-3<=ky<=nY : M(ky, kx, time)=1 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch nargin
    case 0
        error('Size of phase encoding is not specified')
    case 1
        error('Size of frequency encoding is not specified')
    case 2
        error('Total time number is not specified')
    case 3
        error('Acceleration factor is not specified')
end

Tsize=4;
nY2=nY-Tsize*2;

M=zeros(nY2, nframe);

for i=1:nframe
    position_M(:,i)=randn((round(nY/down)-Tsize*2),1)*sqrt(1);
end

Max_rand=max(abs(position_M(:)));
position_M=ceil(position_M*nY2/2/Max_rand)+nY2/2;
position_M(find(position_M < 1))=1;
position_M(find(position_M > nY2))=nY2;
for i=1:nframe
    M(position_M(:,i),i)=1;
end

for i=1:nframe
    residual=round(nY/down)-Tsize*2-sum(M(:,i));
    if residual > 0
        temp_zeros=find(M(:,i)==0);
        n=length(temp_zeros);
        for j=1:residual
            temp=ceil(rand(1)*n);
            M(temp_zeros(temp),i)=1;
            clear temp_zeros;
            temp_zeros=find(M(:,i)==0);
            n=n-1;
        end
    end
end

M2=zeros(nY,nframe);
M2(1:nY/2-4,:)=M(1:nY2/2,:);
M2(nY/2+5:nY,:)=M(nY2/2+1:nY2,:);
M2(nY/2-3:nY/2+4,:)=1;

clear M;
M3=zeros(nY,nframe,nX);
for j=1:nX
    M3(:,:,j)=M2;
end

M3=permute(M3,[1 3 2]);




              