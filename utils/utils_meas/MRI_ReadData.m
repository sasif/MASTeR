vec = @(z) z(:);
mpath = mfilename('fullpath');
% dir_data = [mpath(1:end-length(mfilename)),'/data/'];
dir_data = ['/data/'];

%% DATASET SELECTION
switch (DATASET)
    case 1
        T_frames=12;   %input('    How many frames do you want to use?    ');
        static=96;  %input('    How many Static lines are there in the image    ');
        staticup=43;    %input('    How many Static lines are there in the upper part of the image   ');
        dynamic=192-static;     %input('    How many Dynamic lines there are in the image    ');
        total=dynamic+static;
        start_column=1;
        end_column=224;
        %         reduction_factor = reduction_factor_spacerip* (total/(static/T_frames+dynamic));
        phase_encoding = total;
        frequency_encoding = end_column;
        coils = 8; %input('    How many coils are there in the data    ');
        Frames = 12; %input('    How many frames are there in the data?    ');
        kspacefilename = [dir_data,'k-space_data/P517127_192_224_12_8.cpx'];
        Filename = 'P517127';
        sensitivitymapfilename = [dir_data,'sensitivitymaps/P517127_192_224_12_8_sensitivity_map/sensitivitymaps_192_224_phaselines96imagespace.cpx'];
        
        N = total;
        P = end_column;
        ROW = N;
        COL = P;
        
        % Region of interest
        ROI_ver = [40:140];
        ROI_hor = [40:180];
        
        % Active region for use in motion estimation using CWT
        AR_ver = 1:192;
        AR_hor = 1:224;
        
        % Selects static rows
        STATIC_MASK = zeros(ROW,COL);
        STATIC_MASK(1:staticup,:) = 1;
        STATIC_MASK(ROW-(static-staticup)+1:ROW,:) = 1;
        
    case 2
        T_frames=16;%input('    How many frames do you want to use?    ');
        static=112;%input('    How many Static lines are there in the image    ');
        staticup=42;%input('    How many Static lines are there in the upper part of the image   ');
        dynamic=112;%input('    How many Dynamic lines there are in the image    ');
        total=dynamic+static;
        start_column=1;
        end_column=256;
        %         reduction_factor = reduction_factor_spacerip* (total/(static/T_frames+dynamic));
        phase_encoding = total;
        frequency_encoding = end_column;
        coils = 5; %input('    How many coils are there in the data    ');
        Frames = 19; %input('    How many frames are there in the data?    ');
        kspacefilename = [dir_data,'k-space_data/P501767_224_256_19_5.cpx'];
        Filename = 'P501767';
        if PRESCAN
            % Center taken from another scan of the same object
            sensitivitymapfilename = [dir_data,'sensitivitymaps/P49664_224_256_19_5_sensitiviy_map/sensitivitymaps_224_256_rawsensitivitymap2imagespace.cpx'];
        else
            % Center taken from the same full-scan data
            sensitivitymapfilename = [dir_data,'sensitivitymaps/P50176_224_256_19_5_sensitivity_map/P50176_sensitivitymaps_224_256_5_phaselines112imagespace.cpx'];
        end
        N = total;
        P = end_column;
        ROW = N;
        COL = P;
        
        % Region of interest
        ROI_ver = 40:140;
        ROI_hor = 150:225;        
        ROI_ver = 61:130; 
        ROI_hor = 151:220;
        
        % Active region for use in motion estimation using CWT
        AR_ver = 1:224;
        AR_hor = 1:256;
        
        % Selects static rows
        STATIC_MASK = zeros(ROW,COL);
        STATIC_MASK(1:staticup,:) = 1;
        STATIC_MASK(ROW-(static-staticup)+1:ROW,:) = 1;         
        
        % Frames designated as diastole phase via manual/visual inspection
        DIASTOLE_FRAMES = [1 2 3 14 15 16];
        
    case 3
        T_frames=20;%input('    How many frames do you want to use?    ');
        static=120;%input('    How many Static lines are there in the image    ');
        staticup=44;%input('    How many Static lines are there in the upper part of the image   ');
        dynamic=240-static;%input('    How many Dynamic lines there are in the image    ');
        total=dynamic+static;
        start_column=1;
        end_column=200;
        %         reduction_factor = reduction_factor_spacerip* (total/(static/T_frames+dynamic));
        phase_encoding = total;
        frequency_encoding = end_column;
        coils = 5; %input('    How many coils are there in the data    ');
        Frames = 20; %input('    How many frames are there in the data?    ');
        kspacefilename = [dir_data,'k-space_data/wrap_003_240_200_20_5.cpx'];
        Filename = 'wrap003';
        sensitivitymapfilename = [dir_data,'sensitivitymaps/wrap_003_240_200_20_5_sensitivity_map/sensitivitymaps_240_200_phaselines240imagespace.cpx'];
        
        N = total;
        P = end_column;
        ROW = N;
        COL = P;
                
        % Region of interest
        ROI_ver = [40:160];
        ROI_hor = [20:180];      
        ROI_ver = [56:150];
        ROI_hor = [46:135];
        
        % Active region for use in motion estimation using CWT
        AR_ver = 9:232;
        AR_hor = 5:196;
        
        % Selects static rows
        STATIC_MASK = zeros(ROW,COL);
        STATIC_MASK(1:staticup,:) = 1;
        STATIC_MASK(ROW-(static-staticup)+1:ROW,:) = 1;
        
        % Frames designated as diastole phase via manual/visual inspection
        DIASTOLE_FRAMES = [1 2 13 14 15 16];        
    case 4
        T_frames=12;%input('    How many frames do you want to use?    ');
        static=144;%input('    How many Static lines are there in the image    ');
        staticup=72;%input('    How many Static lines are there in the upper part of the image   ');
        dynamic=288-static;%input('    How many Dynamic lines there are in the image    ');
        total=dynamic+static;
        start_column=1;
        end_column=300;
        % reduction_factor = reduction_factor_spacerip* (total/(static/T_frames+dynamic));
        phase_encoding = total;
        frequency_encoding = end_column;
        coils = 5; %input('    How many coils are there in the data    ');
        Frames = 12; %input('    How many frames are there in the data?    ');
        kspacefilename = [dir_data,'k-space_data/555_rightimage_288_300_12_5.cpx'];
        Filename = '555';
        sensitivitymapfilename = [dir_data,'sensitivitymaps/555_rightimage_288_300_12_5_sensitivity_map/sensitivitymaps_288_300_phaselines144imagespace.cpx'];
        
        N = total;
        P = end_column;
        ROW = N;
        COL = P;
                
        % Region of interest
        ROI_ver = 50:250;
        ROI_hor = 100:250;
        
        % Active region for use in motion estimation using CWT
        AR_ver = 1:288;
        AR_hor = 7:294;
        
        % Selects static rows
        STATIC_MASK = zeros(ROW,COL);
        STATIC_MASK(1:staticup,:) = 1;
        STATIC_MASK(ROW-(static-staticup)+1:ROW,:) = 1;        
    case 6
        T_frames=19;%input('    How many frames do you want to use?    ');
        static=0;%input('    How many Static lines are there in the image    ');
        staticup=0;%input('    How many Static lines are there in the upper part of the image   ');
        dynamic=192;%input('    How many Dynamic lines there are in the image    ');
        total=dynamic+static;
        start_column=1;
        end_column=256;
        %         reduction_factor = reduction_factor_spacerip* (total/(static/T_frames+dynamic));
        phase_encoding = total;
        frequency_encoding = end_column;
        coils = 8; %input('    How many coils are there in the data    ');
        Frames = 19; %input('    How many frames are there in the data?    ');
        if PRESCAN
            kspacefilename = [dir_data,'k-space_data/P48640.7_zeropadding_192_256_19_8.cpx'];
        else
            kspacefilename = [dir_data,'k-space_data/P48128.7_zeropadding_192_256_19_8.cpx'];
        end
        Filename = 'P48128';

        % Center taken from the same full-scan data
        sensitivitymapfilename = [dir_data,'sensitivitymaps/P48128192_256_19_8_sensitivity_map/P48128_sensitivitymaps_192_256_8_phaselines96imagespace.cpx'];
        
        N = total;
        P = end_column;
        ROW = N;
        COL = P;
        
        % Selects static rows
        STATIC_MASK = zeros(ROW,COL);
        STATIC_MASK(1:staticup,:) = 1;
        STATIC_MASK(ROW-(static-staticup)+1:ROW,:) = 1;
        
        % Region of interest
        ROI_ver = 20:160;
        ROI_hor = 30:200;
        ROI_ver = 30:130;
        ROI_hor = 80:190;
        % Active region for use in motion estimation using CWT
        AR_ver = 1:192;
        AR_hor = 1:256;
        
    otherwise
        % error('Invalid Dataset');
end

if exist('PARAM_ONLY')
    if PARAM_ONLY == 1
        return;
    end
end


% Note that we can assin any arbitrary static mask. 
% Following code would designate the ROI to be static
% STATIC_MASK = ones(ROW,COL);
% STATIC_MASK(ROI_ver, ROI_hor) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% DATA READ AND MEASUREMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if DATASET == 5;
    Filename = 'ktDATA';
    load 'full_ktDATA.mat';
    kt = kt/norm(vec(kt(:,:,1)));
    
    [ROW COL T_frames] = size(kt);
    
    sum_of_squares = abs(ifft2(kt))*ROW;
    % kt = fft2(sum_of_squares)/ROW;
    
    coils = 1;
    NbChannels1 = 1;
    
    full_k_space = [];
    SENSITIVITY_MAPS = [];
    for jj = 1:T_frames
        full_k_space(:,:,1,jj) = kt(:,:,jj);
        SENSITIVITY_MAPS(:,:,1,jj) = ones(ROW, COL);
    end

    % Selects static rows
    STATIC_MASK = zeros(ROW,COL);
    %     STATIC_MASK(:,1:75) = 1;
    STATIC_MASK([1:20 175:ROW],:) = 1;
        
    ROI_ver = 50:150;
    ROI_hor = 75:200;
    
    AR_hor = 1:256;
    AR_ver = 1:256;
else
    %%%%% READ THE ORIGINAL K_SPACE
    % FOR EXAMPLE: phantom_256_256_frames_16radio_10coils_4
    % phrase=sprintf('Open the adquired full k_space');
    [full_k_space,NbPhaseEncodingsOrRows,NbFrequencyEncodingsOrCols,NbFrames1,NbChannels1,filename,pathname,PortionOfTheDataFile]=read_file(kspacefilename,phase_encoding,frequency_encoding,Frames,coils);
    
    %%%%%%full_k_space is a 4 D matrix (rows,columns,channel or coil, phase or frame)   256 256 4 16
    % full_k_space is the complete frequency observations of the underlying image
    % We observe each frame using different coils, where each coil puts a complex mask on the underlying image before the Fourier transform
    % The frequency response is centered at the origin.
    
    full_k_space = fft2(fftshift(fftshift(ifft2(fftshift(fftshift(full_k_space,1),2)),2),1));
    full_k_space = full_k_space/norm(vec(full_k_space(:,:,1,1)));
    
    %%%
    %%%% READ THE ORIGINAL Sensitivity map image space%%%%%%%%%
    [sensitivity_maps,NbPhaseEncodingsOrRows,NbFrequencyEncodingsOrCols,NbFrames,NbChannels,filename,pathname,PortionOfTheDataFile]=read_file(sensitivitymapfilename,phase_encoding,frequency_encoding,Frames,coils);
    % Sensitivity map is in the spatial domain.
   
    SENSITIVITY_MAPS = sensitivity_maps;
    clear sensitivity_maps;
end