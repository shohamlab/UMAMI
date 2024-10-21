clear all;close all;
%**************************************
 %%%  Paths & Inputs &  parameters %%%
%**************************************

%%%%%%%%% Paths %%%%%%%%%%
codPath = addpath(genpath(''));               % optional depending on MATLAB directories
direc = '';                                   % Please change path of data
addpath(genpath(direc))
list=dir ([direc, '/', '**/*.tif']);              % find tif files


%%%% Input parameters %%%%%%
nOfImagePerSlice = 200;                    % change based on the # of recorded frame
startingFrame = 31;
USFrameStep = 31;                          % the frame interval between two pulses (eg 121 frames)
nOfPulse = floor(nOfImagePerSlice/USFrameStep);
nOfModulatedFramePerPulse = 6;              % number of the frame in each US pulse
fD = [885 915];                             % DC frequency bandwidth
FA = [979 1009; 1071 1101];                 % AC/modulated frequency bandwidth [first second] harmonic
NumUSFreeFrames = 9;                        % +/- this number of frames before/after US to get a nonUS baseline.
g = 0;                                      % use a GPU for filtering.



%**************************************
 %%%  Frequency Filter %%%
%**************************************

for slice = 1:length(list)                  % file number = slice

    pathi = [fullfile(fileparts(which([list(slice).name]))),'/'];
    nami = erase(list(slice).name,'.tif');
    mkdir( [pathi, 'FigFirst',num2str((slice))] )                     % make newfolder
    pathi = [pathi, 'FigFirst',num2str((slice)),'/'];                 % new path

    nFrame=[];
    for pulse=1:nOfPulse
        nFrame= [nFrame,(startingFrame-1+USFrameStep*(pulse-1)+(1:nOfModulatedFramePerPulse))]; % modulated frames vector
    end
    


    for i = 1:length(nFrame)
        Aout(:,:,i) = tiffreadVolume([direc,'/',nami,'.tif'],'PixelRegion',{[1 inf], [1 1798], [nFrame(i) nFrame(i)]});
        backAout(:,:,i) = tiffreadVolume([direc,'/',nami,'.tif'],'PixelRegion',{[1 inf], [1 1798], [nFrame(i)-NumUSFreeFrames nFrame(i)-NumUSFreeFrames]});
    end


    close all;
 
    A=fftshift(fft(fftshift(Aout(:,:,5)')))'; % Image and FFT of the one modulated frame to check the modulation frequency
    figure; imagesc(Aout(:,:,1)) 
    sA=sum(abs(A),1);
    figure ,plot(sA)
    

    primary = mean(Aout(:,:,1:5),3); 
  
    
    for Frame = 1:length(nFrame)
    

        %%%%%%%%%%%% DC filter and motion corrections
        bImageDC(:,:) = demodulationProcedure(backAout(:,:,(Frame)),fD,g); % DC of non-modulated frame
        a(:,:,1) = primary;
        a(:,:,2) = bImageDC;
        [aa,bshifts,template,options,col_shift]=normcorre(a,NoRMCorreSetParms('d1',size(a,1),...
            'd2',size(a,2),'correct_bidir',0,'max_shift',[50,50,50]),a(:,:,1));
        bImageDC_frame(:,:,Frame) = aa(:,:,2);
    
        modImageDC(:,:) = demodulationProcedure(Aout(:,:,(Frame)),fD,g); % DC of UMAMI frame 
        a(:,:,1) = primary;
        a(:,:,2) = modImageDC;
        [aa,shifts,template,options,col_shift]=normcorre(a,NoRMCorreSetParms('d1',size(a,1),...
            'd2',size(a,2),'correct_bidir',0,'max_shift',[50,50,50]),a(:,:,1));
        modImageDC_frame(:,:,Frame) = aa(:,:,2);
    
                                                                                                                                                                                                                                                                                                                                                                                                                    
        %%%%%%%%%%%% Demodulation of first and second harmonic of frames
        %%%%%%%%%%%% and application of motion correction frame shifts
        options_nc = NoRMCorreSetParms('d1',size(a,1),'d2',size(a,2),'max_shift',...
        [500,500,500]);
        
        AfirstModulImageabs(:,:) = demodulationProcedure(Aout(:,:,(Frame)),FA(1,:),g); % First harmonic of UMAMI frame
        noise_left1(:,:) = demodulationProcedure(Aout(:,:,(Frame)),FA(1,:)-31,g);
        noise_right2(:,:) = demodulationProcedure(Aout(:,:,(Frame)),FA(1,:)+31,g);
        aaa = (noise_left1(:,:) + noise_right2(:,:)) ./ 2;
        temp(:,:,1) = AfirstModulImageabs - aaa;
        temp(:,:,2) = AfirstModulImageabs - aaa;
        aa = apply_shifts(temp,shifts,options_nc); %,0,0);
        Demod_frames(:,:,Frame) = aa(:,:,1);

        AsecondModulImageabs(:,:) = demodulationProcedure(Aout(:,:,(Frame)),FA(2,:),g); % Second harmonic of UMAMI frame
        noise_left2(:,:) = demodulationProcedure(Aout(:,:,(Frame)),FA(2,:)-31,g);
        noise_right2(:,:) = demodulationProcedure(Aout(:,:,(Frame)),FA(2,:)+31,g);
        aaa2 = (noise_left2(:,:) + noise_right2(:,:)) ./ 2;
        temp(:,:,1) = AsecondModulImageabs - aaa2;
        temp(:,:,2) = AsecondModulImageabs - aaa2;
        aa = apply_shifts(temp,shifts,options_nc); %,0,0);
        Demod_frames2nd(:,:,Frame) = aa(:,:,1);

        
        
        
        bAfirstModulImageabs(:,:) = demodulationProcedure(backAout(:,:,(Frame)),FA(1,:),0); % First harmonic of non-mudulated frame
        bnoise_left1(:,:) = demodulationProcedure(backAout(:,:,(Frame)),FA(1,:)-31,g);
        bnoise_right1(:,:) = demodulationProcedure(backAout(:,:,(Frame)),FA(1,:)+31,g);
        baaa = (bnoise_left1(:,:) + bnoise_right1(:,:)) ./ 2;
        temp(:,:,1) = bAfirstModulImageabs - baaa;
        temp(:,:,2) = bAfirstModulImageabs - baaa;
        aa = apply_shifts(temp,bshifts,options_nc); %,0,0);
        bDemod_frames(:,:,Frame) = aa(:,:,1);

        bAsecondModulImageabs(:,:) = demodulationProcedure(backAout(:,:,(Frame)),FA(2,:),0); % Second harmonic of non-modulated frame
        bnoise_left2(:,:) = demodulationProcedure(backAout(:,:,(Frame)),FA(2,:)-31,g);
        bnoise_right2(:,:) = demodulationProcedure(backAout(:,:,(Frame)),FA(2,:)+31,g);
        baaa2 = (bnoise_left2(:,:) + bnoise_right2(:,:)) ./ 2;
        temp(:,:,1) = bAsecondModulImageabs - baaa2;
        temp(:,:,2) = bAsecondModulImageabs - baaa2;
        aa = apply_shifts(temp,bshifts,options_nc); %,0,0);
        bDemod_frames2nd(:,:,Frame) = aa(:,:,1);
        

    end

    
    %% Processed image smoothing
    average = 10;
    convmat = ones(average);convmat=convmat./sum(convmat(:));
    
    
    modImageDC(:,:,slice) = conv2(mean(modImageDC_frame(:,:,:),3),convmat,'same');
    bImageDC(:,:,slice) = conv2(mean(bImageDC_frame(:,:,:),3),convmat,'same');
    
    DemodImage(:,:,slice) = conv2(mean(Demod_frames(:,:,:),3),convmat,'same');
    DemodImage2nd(:,:,slice) = conv2(mean(Demod_frames2nd(:,:,:),3),convmat,'same');
    
    bDemodImage(:,:,slice) = conv2(mean(bDemod_frames(:,:,:),3),convmat,'same');
    bDemodImage2nd(:,:,slice) = conv2(mean(bDemod_frames2nd(:,:,:),3),convmat,'same');
    
    
    
    %%%%%%%%%%%% non-motion corrected nonmodulated mean image 
    Original(:,:,slice) = conv2(mean(backAout,3),convmat,'same');




    %% Save all single frame MAT files and images from single recording
    b = 50;             % Number of pixels to crop filter effects at edges of images
    fsize = 18;         % Image text font size
    
    
    fpath=fullfile(pathi, [nami,'modImageDC_frame.mat']);
    save(fpath,'modImageDC_frame') 
    fpath=fullfile(pathi, [nami,'bImageDC_frame.mat']);
    save(fpath,'bImageDC_frame') 

    fpath=fullfile(pathi, [nami,'Demod_frames.mat']);
    save(fpath,'Demod_frames') 
    fpath=fullfile(pathi, [nami,'Demod_frames2nd.mat']);
    save(fpath,'Demod_frames2nd') 

    fpath=fullfile(pathi, [nami,'bDemod_frames.mat']);
    save(fpath,'bDemod_frames') 
    fpath=fullfile(pathi, [nami,'bDemod_frames2nd.mat']);
    save(fpath,'bDemod_frames2nd') 

    
    figure; imagesc(modImageDC(b:end-b,b:end-b,slice)); colormap(gray); axis off; c = colorbar; set(c,'FontSize',fsize);
    saveas(gcf,[pathi,'modimageDC.tif'])
    
    figure; imagesc(bImageDC(b:end-b,b:end-b,slice)); colormap(gray); axis off; c = colorbar; set(c,'FontSize',fsize);
    saveas(gcf,[pathi,'bImageDC.tif'])


    figure; imagesc(DemodImage(b:end-b,b:end-b,slice)); colormap(gray); axis off; c = colorbar; set(c,'FontSize',fsize);
    saveas(gcf,[pathi,'DemodImage.tif'])

    figure; imagesc(DemodImage2nd(b:end-b,b:end-b,slice)); colormap(gray); axis off; c = colorbar; set(c,'FontSize',fsize);
    saveas(gcf,[pathi,'DemodImage2nd.tif'])
    
    
    figure; imagesc(bDemodImage(b:end-b,b:end-b,slice)); colormap(gray); axis off; c = colorbar; set(c,'FontSize',fsize);
    saveas(gcf,[pathi,'bDemodImage.tif'])
    
    figure; imagesc(bDemodImage2nd(b:end-b,b:end-b,slice)); colormap(gray); axis off; c = colorbar; set(c,'FontSize',fsize);
    saveas(gcf,[pathi,'bDemodImage2nd.tif'])
   

    figure; imagesc(mean(Original(b:end-b,b:end-b,slice),3)); colormap(gray); axis off; c = colorbar; set(c,'FontSize',fsize);
    saveas(gcf,[pathi,'Original.tif'])
    
  

end

%% Save MAT files for the mean images generated above. 
% Saves images after all recordings are processed in a new folder. Matrices
% are in file read order along matrix 3rd dimension

mkdir ([direc, '/', 'mat files'])
pathi = [direc, '/', 'mat files', '/'];

fpath=fullfile(pathi, ['DemodImage','.mat']);
save(fpath,'DemodImage')
fpath=fullfile(pathi, ['DemodImage2nd','.mat']);
save(fpath,'DemodImage2nd')

fpath=fullfile(pathi, ['bDemodImage','.mat']);
save(fpath,'bDemodImage')
fpath=fullfile(pathi, ['bDemodImage2nd','.mat']);
save(fpath,'bDemodImage2nd')

fpath=fullfile(pathi, ['modImageDC','.mat']);
save(fpath,'modImageDC')
fpath=fullfile(pathi, ['bImageDC','.mat']);
save(fpath,'bImageDC')

fpath=fullfile(pathi, ['Original','.mat']);
save(fpath,'Original') 

    


