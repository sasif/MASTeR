function [RowsPerFrame, matrixdisplay, rows_to_be_samplednew]=samplefinalrowsenglish(static,  staticup, dynamic, reduction_factor, rows_to_be_sampled,frames_you_want_to_use)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%           JAVI ACEBRON FABREGAT    14-09-2006
%               LAST VERSION
%  THIS FUNCTION SELECT THE LINES FOR NOQUIST, BUT AMONG THE LINES I
%  SELECTED FOR SENSE BEFORE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % THIS FUNCTION CALLS FUNCTIONS: 
   % THIS FUNCION IS CALLED BY FUNCTIONS: sensenoquist....m
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%  IF REDUCTION FACTOR IS 4, THE INPUT WAS: 4,8,12,16...


to_see_lines_chosen_graphically= 0; %input('    Do you want to see the chosen lines(0=NO ; 1=YES)?    ');
save_matrix=0; %input('    Do you want to save the image of the matrix with the lines you chose(0=NO ; 1=YES)?   '); 
   
total=static+dynamic;
N_over_ND=total/dynamic;                                            
TotalNbOfUnknowns=static+dynamic*frames_you_want_to_use;
formatrixsense=rows_to_be_sampled;                                          % I NEED THIS VALUES LATER, I STORE THEM


%%%%% SELECT THE ODD LINES: I TAKE FOR ALL THE FRAMES: 4,12,20,...
another_line=zeros(size(rows_to_be_sampled,2),1);
i=1;
for step_dynamic=0:ceil(dynamic/reduction_factor)-1                        % for some combinations (with rf=3, for example 4fr,50st,25stup,206dy o 4fr,21st,10stup,43dy) you have to put floor instead of ceil, and ceil in the other for
    one_component=rows_to_be_sampled(1+round(step_dynamic*N_over_ND));      
    rows_to_be_sampled(1+round(step_dynamic*N_over_ND))=0;                  
    another_line(one_component,1)=one_component;
end
for i=1:frames_you_want_to_use
    RowsPerFrame(:,i)=another_line;
end
rows_to_be_samplednew=rows_to_be_sampled;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%% SELECT THE EVEN LINES (FIRST FRAME: 8, SECOND FRAME: 16, THIRD FRAME: 24....)
frame=frames_you_want_to_use+1;
for step_static=0:floor(static/reduction_factor)-1    % PROBLEM FOR RF=3, st32 stup16 dy96 IT DONT TAKE THE LAST LINE , IF I PUT CEIL INSTEAD OR FLOOR, THEN IT WORKS FOR THIS CASE, BUT NOT FOR THE OTHERS
    if frame==(frames_you_want_to_use+1)
        frame=1;
    end
    if step_static==0
        frame=1;
        RowsEqualToZeroThisFrame=find(RowsPerFrame(:,frame)==0);
    end
    muchosnozeros=find(rows_to_be_samplednew);
    if length(muchosnozeros)>0
        unonozero=muchosnozeros(1);
        RowsPerFrame((RowsEqualToZeroThisFrame(1,1)),frame)=rows_to_be_samplednew(unonozero);
        rows_to_be_samplednew(unonozero)=0;
        RowsEqualToZeroThisFrame=RowsEqualToZeroThisFrame(2:end,:);
        frame=frame+1;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%% TO BUILT matrixdisplay, THAT WILL BE THE IMPORTANT OUTPUT OF THE FUNCTION
matrixdisplay=[zeros(total,9)];
all_ones=zeros(total,9);
for frame=1:frames_you_want_to_use
    onecolumn=zeros(total,1); 
    RowsPerFrameWithoutZeros=nonzeros(RowsPerFrame(:,frame));
    for i=1:length(RowsPerFrameWithoutZeros)
        onecolumn(RowsPerFrameWithoutZeros(i),1)=1;
    end
matrixdisplay=[matrixdisplay onecolumn zeros(total,9)];
all_ones=[all_ones ones(total,1) zeros(total,9)];
end
matrixsense=zeros(total, 9);
onecolumn=zeros(total,1);
onecolumn(formatrixsense,1)=formatrixsense;
for frame=1:frames_you_want_to_use
    matrixsense=[matrixsense onecolumn zeros(total,9)];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%% IF YOU WANT TO SEE THE CHOSEN LINES
if to_see_lines_chosen_graphically==1
    figure                                                                  
    spy(all_ones)                                                           % ALL THE FRAMES. YOU SEE THE LINES YOU DONT HAVE TO TAKE BECAUSE YOU ARE USING SENSE (IN BLUE)
    hold on
    spy(matrixsense,'y')                                                    % YOU SEE THE LINES YOU DONT HAVE TO TAKE BECAUSE YOU ARE USING NOQUIST IN YELLOW
    hold on
    spy(matrixdisplay,'r')                                                  % SHOW THE LINES CHOSEN IN EACH FRAME (10=FRAME 1) IN RED
    hold off
    title('I dont take the BLUE lines because of SENSE and the YELLOW lines because of NOQUIST. I TAKE THE RED LINES');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        
%%%%% IF YOU WANT TO SAVE THE MATRIX WHERE YOU CAN SEE WHICH LINES YOU TAKE
if save_matrix==1
    NQDATA=getenv('NQDATA');
    NQOUTDATA=[NQDATA,'/out'];
    status=mkdir(NQOUTDATA,'matrix_for_sense&noquist');
    cd (NQOUTDATA);
    cd 'matrix_for_sense&noquist';
    ChaineName_Image=['matrix_for_sense&noquist_st',sprintf('%d',static),'_stup',sprintf('%d',staticup),'_dy',sprintf('%d',dynamic),'_rf',sprintf('%d',reduction_factor),'_fr' ,sprintf('%d',frames_you_want_to_use),'.bmp'];
    imwrite(matrixdisplay,ChaineName_Image,'bmp');
    NQOUTDATA=[NQDATA,'/out'];
    cd ..;
    cd ..;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




