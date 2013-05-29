function rows_to_be_sampled=samplerowsenglish(original_k_space, reduction_factor)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%           JAVI ACEBRON FABREGAT    14-09-2006
%               LAST VERSION
%  THIS FUNCTION TAKES THE LINES TO MAKE THE SENSE STEP (1 in 2, 1 in 3 or 1 in 4)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % THIS FUNCTION CALLS FUNCTIONS: 
   % THIS FUNCION IS CALLED BY FUNCTIONS: sensenoquist...M
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
[P_codifications,F_codifications]=size(original_k_space(:,:,1,1));
if mod(reduction_factor,1)==0
    if mod(P_codifications,2)==0                                                % EVEN NUMBER OF PHASE CODIFICATIONS 
        middle=P_codifications/2;                                               % IF THERE ARE 6 LINES, IT TAKES THE 3th (THE 4th WOULD BE GOOD TOO TO TAKE)
        up=middle:-reduction_factor:1;
        up=fliplr(up);
        down=middle+reduction_factor:reduction_factor:P_codifications; 
    else                                                                        % ODD NUMBER OF PHASE CODIFICATIONS
        middle=ceil(P_codifications/2);                                         % UPPER INTEGER. IT TAKES ALWAYS THE CENTER LINE
        up=middle:-reduction_factor:1;
        up=fliplr(up);
        down=middle+reduction_factor:reduction_factor:P_codifications;
    end
    rows_to_be_sampled=[up down];
else                                                                            % REDUCTION FACTOR DECIMAL
    rows_to_be_sampled=round(linspace(1,P_codifications,round(P_codifications/reduction_factor)));
end