%% Copyright(C) 2018
%  Developed By: Muhammet Emin Yanik
%  Advisor: Prof. Murat Torlak
%  The University of Texas at Dallas
%  Department of Electrical and Computer Engineering

%  This work was supported by the Semiconductor Research Corporation (SRC) task 2712.029
%  through The University of Texas at Dallas' Texas Analog Center of Excellence (TxACE).

%  Redistributions and use of source must retain the above copyright notice
%  Redistributions in binary form must reproduce the above copyright notice
%%

function [sarImage,xRangeT_mm,yRangeT_mm,zRangeT_mm] = reconstructSARimageFFT_3D(sarData,frequency,xStepM,yStepM,xySizeT,zTarget,nFFTkXY)
% For wideband processing:
% -------------------------------------------------------------------------
% sarData: should be yPointM x xPointM x nSample


% Example function calls, see details below
% -------------------------------------------------------------------------
% lambda/4 = 0.9487, f0 = 77e9+6e-6*63.343e12;
% [77e9,150e12,3000e3,4e-6] = [fStart(Hz),fSlope(Hz/s),fSample(sps),adcStart(s)]
% [sarImage,xRangeT,yRangeT,zRangeT] = reconstructSARimageFFT_3D(sarData,[77e9,70.295e12,5000e3,4.66e-6],400/407,7.59/8,-1,[200:5:250],512);

% For wideband processing:
% -------------------------------------------------------------------------
% frequency: [fStart,fSlope,fSample,nSample]
% fStart: Start frequency
% fSlope: Slope const (Hz/sec)
% fSample: Sample ps
% Example: [77e9,63.343e12,9121e3]

% Variables
% -------------------------------------------------------------------------
% xStepM: measurement step size at x (horizontal) axis in mm
% yStepM: measurement step size at y (vertical) axis in mm
%
% xySizeT: size of target area in mm
% zTarget: target distance in mm
% nFFTkXY: number of FFT points, should be greater than xStepM and yStepM


%% Code Starts
% profile on


%% Define Fixed Parameters
%-------------------------------------------------------------------------%
isAmplitudeFactor = true; % Set true if Amplitude Factor is needed for 2D

is3DImaging = false;
is2DImaging = false;

if length(zTarget) > 1 % is3DImaging is true if depth data is given
    is3DImaging = true;
else
    is2DImaging = true;
end


%% Define Frequency Spectrum
%-------------------------------------------------------------------------%
[~,~,nSample] = size(sarData); % Number of samples
if (length(frequency)>1) && (length(frequency)<=4) && (nSample>1)
    frequency = num2cell(frequency);
    [f0,K,fS,adcStart] = frequency{:};
    f0 = f0 + adcStart*K; % This is for ADC sampling offset
    f = f0 + (0:nSample-1)*K/fS; % wideband frequency
else
    error('Please correct the configuration and data for 3D processing')
end


%% Define Fixed Parameters
%-------------------------------------------------------------------------%
c = 299792458; % physconst('lightspeed'); in m/s
k = 2*pi*f/c;
k = reshape(k,1,1,[]);


%% Coincide Aperture and Target Domains
%-------------------------------------------------------------------------%
[yPointM,xPointM,~] = size(sarData);
xStepT = xStepM;
yStepT = yStepM;
zRangeT_mm = zTarget * 1e-3;


%% Define Number of FFT Points
%-------------------------------------------------------------------------%
if (nFFTkXY<xPointM) || (nFFTkXY<yPointM)
    warning('# of FFT points should be greater than the # of measurement points. FFT will be performed at # of measurement points')
end
% Set nFFTkX and nFFTkY accordingly
if (nFFTkXY>xPointM)
    nFFTkX = nFFTkXY;
else
    nFFTkX = xPointM;
end
if (nFFTkXY>yPointM)
    nFFTkY = nFFTkXY;
else
    nFFTkY = yPointM;
end


%% Define Wavenumbers
%-------------------------------------------------------------------------%
wSx = 2*pi/(xStepT*1e-3); % Sampling frequency for Target Domain
kX = linspace(-(wSx/2),(wSx/2),nFFTkX); % kX-Domain

wSy = 2*pi/(yStepT*1e-3); % Sampling frequency for Target Domain
kY = (linspace(-(wSy/2),(wSy/2),nFFTkY)).'; % kY-Domain


%% Zero Padding to sarData to Locate Target at Center
%-------------------------------------------------------------------------%
% Prepare the zero padded data
sarDataPadded = single(zeros(nFFTkY,nFFTkX,nSample));

% Prepare the index of the data portion in x and y
indexZeroPadStart_x = floor((nFFTkX - xPointM)/2) + 1;
indexZeroPad_x = indexZeroPadStart_x : indexZeroPadStart_x + xPointM - 1;
indexZeroPadStart_y = floor((nFFTkY - yPointM)/2) + 1;
indexZeroPad_y = indexZeroPadStart_y : indexZeroPadStart_y + yPointM - 1;

% Fill the zero padded data
sarDataPadded(indexZeroPad_y, indexZeroPad_x, :) = single(sarData);

% Free the memory, sarData is not needed anymore
clear sarData;


%% Calculate kZ
%-------------------------------------------------------------------------%
kZ = single(sqrt((2*k).^2 - kX.^2 - kY.^2));


%% Take 2D FFT of SAR Data
%-------------------------------------------------------------------------%
sarDataFFT = fftshift(fftshift(fft2(sarDataPadded),1),2);
clear sarDataPadded;

%% Create 2D-SAR Image for single Z
%-------------------------------------------------------------------------%
if is2DImaging
    phaseFactor = exp(-1i*zRangeT_mm*kZ);
    phaseFactor((kX.^2 + kY.^2) > (2*k).^2) = 0;
    
    if ~isAmplitudeFactor
        clear kZ;
    else
        sarDataFFT = kZ .* sarDataFFT;
        clear kZ;
    end
    
    sarDataFFT = sarDataFFT .* phaseFactor;
    
    sarImage = ifft2(sarDataFFT);
    sarImage = sum(sarImage,3);
end


%% Manual Z-Focusing, Create 3D-SAR Image
%-------------------------------------------------------------------------%
if is3DImaging
    [ySizeData,xSizeData,~] = size(sarDataFFT);
    sarImageIFFT = zeros(ySizeData,xSizeData,length(zRangeT_mm));
    for n = 1:length(zRangeT_mm)
        phaseFactor = exp(-1i*zRangeT_mm(n)*kZ);
        phaseFactor((kX.^2 + kY.^2) > (2*k).^2) = 0;
        
        sarDataFFT_Corrected = sarDataFFT .* phaseFactor;
        
        if isAmplitudeFactor
            amplitudeFactor = kZ;
            amplitudeFactor((kX.^2 + kY.^2) > (2*k).^2) = 0;
            sarDataFFT_Corrected = sarDataFFT_Corrected .* amplitudeFactor;
        end
        
        sarImageIFFT(:,:,n) = sum(sarDataFFT_Corrected,3);
    end
    sarImage = ifft2(sarImageIFFT);
end


%% Define Target Axis
%-------------------------------------------------------------------------%
xRangeT_mm = xStepT * (-(nFFTkX-1)/2 : (nFFTkX-1)/2); % xStepM is in mm
yRangeT_mm = yStepT * (-(nFFTkY-1)/2 : (nFFTkY-1)/2); % xStepM is in mm


%% Flip Target in x-Axis
sarImage = flip(sarImage,2);


%% Crop the Image for Related Region
%-------------------------------------------------------------------------%
if (xySizeT ~= -1)
    
    if max(xRangeT_mm)<(xySizeT/2)
        xySizeT = 2*max(xRangeT_mm);
    end
    if max(yRangeT_mm)<(xySizeT/2)
        xySizeT = 2*max(yRangeT_mm);
    end
    
    indXpartT = xRangeT_mm>(-xySizeT/2) & xRangeT_mm<(xySizeT/2);
    indYpartT = yRangeT_mm>(-xySizeT/2) & yRangeT_mm<(xySizeT/2);
    
    xRangeT_mm = xRangeT_mm(indXpartT);
    yRangeT_mm = yRangeT_mm(indYpartT);
    if ~is2DImaging
        sarImage = sarImage(indYpartT,indXpartT,:);
    else
        sarImage = sarImage(indYpartT,indXpartT);
    end
end

%% Plot SAR Image
%-------------------------------------------------------------------------%
if is2DImaging
    figure('OuterPosition',[695 166 670 712]);
    mesh(xRangeT_mm,yRangeT_mm,abs(squeeze(sarImage)),'FaceColor','interp','LineStyle','none')
    view(2)
    colormap('gray');
    
    xlabel('Horizontal (mm)')
    ylabel('Vertical (mm)')
    titleFigure = "SAR 2D Image - " + zTarget + "mm Focused" ;
    title(titleFigure)
end

% profile viewer