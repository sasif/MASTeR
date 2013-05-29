function [FileRead,NbPhaseEncodingsOrRows,NbFrequencyEncodingsOrCols,NbPhases,NbChannels,filename,pathname,PortionOfTheDataFile]=read_file(enunciado,NbPhaseEncodingsOrRows,NbFrequencyEncodingsOrCols,NbPhases,NbChannels);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   MODULE:     read_file.m
%
%   CONTENTS:   reads a raw-data / complex data file
%
%   UPDATE HISTORY:
%       last update:    January 22, 2003
%
%   COPYRIGHT DAVID MORATAL-PEREZ, 2003
%
%   COMMENTS:   the data file read is saved in a 4D matrix.
%               The 3rd and the 4th dimension are the channel number and
%               the frame number respectively.
%               IS THE ONE I USED WITH SENSENOQUIST ENGLISH_LIKEDAVID AND SENSENOQUISTENGLISH_EASY_WAY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Initialization
% clear all;
% close all;
% clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %read the NQDATA environment variable and derive data location strings
% NQDATA=getenv('NQDATA');                   
% NQINDATA=[NQDATA,'/raw'];  %\ era /   ( / para linux, \ para windows)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% disp('     ');
% disp('  Read a raw-data/complex-image file    ');
% disp(' ------------------------------------   ');
% disp('     ');
% disp('     ');


% BEGIN of the INTERVIEW
%%%%%%%%%%%%%%%%%%%%%%%%

%setenv('NQDATA', 'C:\houlei\research\Noquist\SENSE\sense_realdata');
%NQDATA=getenv('NQDATA');                   
%NQINDATA=[NQDATA,'\raw\'];  

% NQINDATA=['C:/houlei/research/Noquist/SENSE/sense_realdata/raw']; 
%cd(NQINDATA);



% Load a raw-data file
% [filename, pathname] = uigetfile( ...             %esto hace qeu salga la ventana para seleccionar el fichero
%     {'*.*', 'All files (*.*)'; ...
%         '*.cpx', 'All raw-data files (*.cpx)'}, ...
%     sprintf(enunciado));
filename = enunciado;
pathname = '';


% disp('     ');
% disp(['  - Using the raw-data/complex-image file ', filename]);
% disp('     ');
% disp('     ');





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%05-29-07
% NbPhaseEncodingsOrRows=128 ; %input('  Nb. of PHASE ENCODINGS/ROWS per PHASE of the raw-data/complex-image file:  ');    %le metes todos estas variables
% NbFrequencyEncodingsOrCols=128; %input('  Nb. of FREQUENCY ENCODINGS/COLUMNS per PHASE of the raw-data/complex-image file:  ');
% NbPhases=input('  Nb. of PHASES of the sequence:  ');
% NbChannels=4; %input('  Nb. of CHANNELS of the raw-data/complex-image file:  ');
% disp('     ');
MoreThanOneDataSet=0; %input('  Are there more than one data set in the same data file? (NO=0 / YES=1)  ');   %�A QUE SE REFIERE CON DATA SET?
% Initialization of the PortionOfTheDataFile variable
PortionOfTheDataFile=0;

if MoreThanOneDataSet==1
    PortionOfTheDataFile=input('    Which portion of the data file do you want to read now? (1st portion=1 / 2nd portion=2 / ...)  ');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NbPhaseEncodingsOrRows=168;
% NbFrequencyEncodingsOrCols=256;
% 
% NbPhases=input('  Nb. of PHASES of the sequence:  ');
% NbChannels=4;
% disp('     ');
% MoreThanOneDataSet=0;   %�A QUE SE REFIERE CON DATA SET?
% % Initialization of the PortionOfTheDataFile variable
% PortionOfTheDataFile=0;
% 
% if MoreThanOneDataSet==1
%     PortionOfTheDataFile=input('    Which portion of the data file do you want to read now? (1st portion=1 / 2nd portion=2 / ...)  ');
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%05-29-07
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%
% END of the INTERVIEW

% A timer to know how many time is needed to do the read
%tic;   %I HAVE COMENTED IT BECAUSE IT�S OUTSIDE TOO

% Read the raw-data file our data file (32 bits, big-endian)
% disp(' Reading the raw-data/complex-image file...');
% disp('    ');
fid = fopen([pathname filename],'r','ieee-be');  %abre el archivo para leer y con formato punto flotante y big endian
% fid

if MoreThanOneDataSet==0
    DataFile = fread(fid,'float32');    %lee el fichero qeu es de punto flotante, 32 bits.
else
    % The variable NbComponentsFileToBeRead is important for the big data sets.
    %  Using this variable there will be no need in loading all the whole file if
    %  only a portion of the file interest us
    OffsetInTheFileToBeRead=(NbPhaseEncodingsOrRows*NbFrequencyEncodingsOrCols*NbPhases*NbChannels*2*4)*(PortionOfTheDataFile-1);
    disp(['...read_file: skipping ',sprintf('%d',OffsetInTheFileToBeRead),' bytes...']);

    fseek(fid,OffsetInTheFileToBeRead,-1);
    NbComponentsFileToBeRead=NbPhaseEncodingsOrRows*NbFrequencyEncodingsOrCols*NbPhases*NbChannels*2;
    DataFile = fread(fid,NbComponentsFileToBeRead,'float32');
    
end

%si_DataFile=size(DataFile)

clear fid

% Let's display some information about this reconstruction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %clc
% disp('     ');
% disp('  Read a raw-data/complex-image file    ');
% disp(' ------------------------------------   ');
% disp('     ');
% disp('     ');
% disp(['  - Using the raw-data/complex-image file ',pathname filename]);
% disp('     ');
% disp(['  - Nb. of PHASE ENCODINGS/ROWS per PHASE of the raw-data/complex-image file:         ',sprintf('%s',num2str(NbPhaseEncodingsOrRows))]);
% disp(['  - Nb. of FREQUENCY ENCODINGS/COLUMNS per PHASE of the raw-data/complex-image file:  ',sprintf('%s',num2str(NbFrequencyEncodingsOrCols))]);
% disp(['  - Nb. of PHASES of the sequence:                                                   ',sprintf('%s',num2str(NbPhases))]);
% disp(['  - Nb. of CHANNELS of the raw-data/complex-image file:                               ',sprintf('%s',num2str(NbChannels))]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END of displaying some information about the reconstruction
% 
% disp('   ');
% disp('   ');
% disp('    Preparing the raw-data/complex-image file... ');
% disp('   ');

% Preparation of the raw-data file into phases and channels

% The raw-data file is a loooong vector that must be put into a 4D matrix, where
%  the first and second dimension indicates phase and frequency encoding of each k-space,
%  and the fourth dimension indicates the phase number  and the third dimension
%  indicates the channel number of each one of those k-spaces.

% The original raw-data file is in .cpx format, that means that the data is complex and
%  it's structured as: real1 complex1 real2 complex2 real3 complex3... That's why the
%  data is divided into 2 4D matrix. One matrix for the real part and the other one for
%  the imaginary part. Then those 2 matrices will be added up into a single complex 4D matrix.
% That also means that for 512 frequency encodings, for example, the step to fill each k-space
%  phase encoding is 512*2=1024:
%    1024 (total) components = 512 (real) components + 512 (imag) components

% DOUBLE LOOP (on the phases and on the channels)
% The third component of the 4D matrix indicates the channel number
% The fourth component of the 4D matrix indicates the phase number

% Just to not obtain an error with PortionOfTheDataFile=0
%  After the loop, PortionOfTheDataFile will be again 0
%if MoreThanOneDataSet==0
%    PortionOfTheDataFile=1;
%end

% Begin loop over the phases
for Phase=1:NbPhases
    
    % Begin loop over the channels
    for Channel=1:NbChannels
        
        % We decompose the whole vector into a big data matrix (F), doing this frame by frame
        %  (first of all we decompose the .cpx data into a real and a imaginary frame that we will add later)
        for i=1:NbPhaseEncodingsOrRows
            
            Real_DataFile(i,:,Channel,Phase)=...
                DataFile(( 1+(i-1)*(NbFrequencyEncodingsOrCols*2) ) + (Phase-1)*( NbPhaseEncodingsOrRows*2*NbFrequencyEncodingsOrCols )+(Channel-1)*( ( NbPhaseEncodingsOrRows*2*NbFrequencyEncodingsOrCols )*NbPhases )...
                :2:((NbFrequencyEncodingsOrCols*2)+(NbFrequencyEncodingsOrCols*2)*(i-1)-1)+(Phase-1)*( NbPhaseEncodingsOrRows*2*NbFrequencyEncodingsOrCols )+(Channel-1)*( ( NbPhaseEncodingsOrRows*2*NbFrequencyEncodingsOrCols )*NbPhases ) ,1)';
            Imag_DataFile(i,:,Channel,Phase)=...
                DataFile(( 2+(i-1)*(NbFrequencyEncodingsOrCols*2) ) + (Phase-1)*( NbPhaseEncodingsOrRows*2*NbFrequencyEncodingsOrCols )+(Channel-1)*( ( NbPhaseEncodingsOrRows*2*NbFrequencyEncodingsOrCols )*NbPhases )...
                :2:((NbFrequencyEncodingsOrCols*2)+(NbFrequencyEncodingsOrCols*2)*(i-1))  +(Phase-1)*( NbPhaseEncodingsOrRows*2*NbFrequencyEncodingsOrCols )+(Channel-1)*( ( NbPhaseEncodingsOrRows*2*NbFrequencyEncodingsOrCols )*NbPhases ) ,1)';
            
            
            %Real_DataFile(i,:,Channel,Phase)=...
            %    DataFile(( 1+(i-1)*(NbFrequencyEncodingsOrCols*2) ) + (Phase-1)*( NbPhaseEncodingsOrRows*2*NbFrequencyEncodingsOrCols )+(Channel-1)*( ( NbPhaseEncodingsOrRows*2*NbFrequencyEncodingsOrCols )*NbPhases ) + ((PortionOfTheDataFile-1)*NbPhaseEncodingsOrRows*NbFrequencyEncodingsOrCols*NbPhases*NbChannels*2)...
            %    :2:((NbFrequencyEncodingsOrCols*2)+(NbFrequencyEncodingsOrCols*2)*(i-1)-1)+(Phase-1)*( NbPhaseEncodingsOrRows*2*NbFrequencyEncodingsOrCols )+(Channel-1)*( ( NbPhaseEncodingsOrRows*2*NbFrequencyEncodingsOrCols )*NbPhases ) + ((PortionOfTheDataFile-1)*NbPhaseEncodingsOrRows*NbFrequencyEncodingsOrCols*NbPhases*NbChannels*2) ,1)';
            %Imag_DataFile(i,:,Channel,Phase)=...
            %    DataFile(( 2+(i-1)*(NbFrequencyEncodingsOrCols*2) ) + (Phase-1)*( NbPhaseEncodingsOrRows*2*NbFrequencyEncodingsOrCols )+(Channel-1)*( ( NbPhaseEncodingsOrRows*2*NbFrequencyEncodingsOrCols )*NbPhases ) + ((PortionOfTheDataFile-1)*NbPhaseEncodingsOrRows*NbFrequencyEncodingsOrCols*NbPhases*NbChannels*2)...
            %    :2:((NbFrequencyEncodingsOrCols*2)+(NbFrequencyEncodingsOrCols*2)*(i-1))  +(Phase-1)*( NbPhaseEncodingsOrRows*2*NbFrequencyEncodingsOrCols )+(Channel-1)*( ( NbPhaseEncodingsOrRows*2*NbFrequencyEncodingsOrCols )*NbPhases ) + ((PortionOfTheDataFile-1)*NbPhaseEncodingsOrRows*NbFrequencyEncodingsOrCols*NbPhases*NbChannels*2),1)';
            
        end
        
        % If the raw-data file has been already completely read, the non useful variables are deleted
        if Phase==NbPhases & Channel==NbChannels
            clear DataFile
            clear Phase
            clear Channel
            clear i
        end
        
    end
    % End loop over the channels
    
end
% End loop over the phases

% Just to not obtain an error with PortionOfTheDataFile=0
%  After the loop (that's here), PortionOfTheDataFile will be again 0
%if MoreThanOneDataSet==0
%    PortionOfTheDataFile=0;
%end


% Addition of the real and imaginary parts of the k-space
FileRead=Real_DataFile+j*Imag_DataFile;

clear Real_DataFile
clear Imag_DataFile

