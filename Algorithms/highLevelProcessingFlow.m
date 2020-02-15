%% Copyright(C) 2018
%  Developed By: Muhammet Emin Yanik
%  Advisor: Prof. Murat Torlak
%  Department of Electrical and Computer Engineering

%  This work was supported by the Semiconductor Research Corporation (SRC) task 2712.029
%  through The University of Texas at Dallas' Texas Analog Center of Excellence (TxACE).

%  Redistributions and use of source must retain the above copyright notice
%  Redistributions in binary form must reproduce the above copyright notice
%%

%% Run this script at 'Algorithms' folder


%% Modify all necessary directories and files
%--------------------------------------------------------------------------
% Path for the Data
experiment_folderName = '..\RecordedData_Example\20190523_smallPlatform_scissorsInBox_400By400';
adcBinData_subFolderName = 'Data';

% Path for the Calibration Data (and optional Background radiation)
calData_fileName = '.\calibration\calData\calData_3Tx_4Rx.mat';
delayOffset_fileName = '.\calibration\calData\delayOffset_3Tx_4Rx.mat';

%% Automated folder and file name generation
%--------------------------------------------------------------------------
addpath(genpath(pwd));

% Load sensorParams and sarParams
load(fullfile(experiment_folderName,adcBinData_subFolderName,'sensorParams'));
load(fullfile(experiment_folderName,adcBinData_subFolderName,'sarParams'));

% Load calData and delayOffset
load(fullfile(calData_fileName));
load(fullfile(delayOffset_fileName));


%% Mean Target z-Distance
zTarget_mm = 250; % update the distance as needed.

%% Data Read Function call
%--------------------------------------------------------------------------
load(fullfile(experiment_folderName,adcBinData_subFolderName,'rawData'));
% rawData Format:
% (Num_RX * Num_TX) * Num_verticalScan * Num_horizontalScan * Samples_per_Chirp;


%% Calibrate rawData
%--------------------------------------------------------------------------
% rawData format should be: (Num_RX x Num_TX) x Num_verticalScan x Num_horizontalScan x Samples_per_Chirp;
rawDataCal = calibrateDataFunction(rawData,sensorParams,calData,delayOffset);
clear rawData


%% Define Parameters
%--------------------------------------------------------------------------
frequency = [sensorParams.Start_Freq_GHz*1e9,sensorParams.Slope_MHzperus*1e12,sensorParams.Sampling_Rate_ksps*1e3,sensorParams.Adc_Start_Time_us*1e-6];
c = physconst('lightspeed');
Samples_per_Chirp = sensorParams.Samples_per_Chirp;
Num_TX = sensorParams.Num_TX;
Num_RX = length(sensorParams.RxToEnable);
Num_horizontalScan = sarParams.Num_horizontalScan;
Num_verticalScan = sarParams.Num_verticalScan;
yStepM_mm = sarParams.Vertical_stepSize_mm;
xStepM_mm = sarParams.Platform_Speed_mmps * sensorParams.Frame_Repetition_Period_ms*1e-3;
lambda_mm = c/79e9*1e3; %  % center frequency


%% Convert multistatic data to monostatic version
%--------------------------------------------------------------------------
% rawData format should be: (Num_RX * Num_TX) * Num_verticalScan * Num_horizontalScan * Samples_per_Chirp;
rawDataMonostatic = convertMultistaticToMonostatic(rawDataCal,frequency,xStepM_mm,yStepM_mm,zTarget_mm,'IWR1443',ones(1,Num_TX),ones(1,Num_RX));
clear rawDataCal


%% Make Uniform Virtual Array
%--------------------------------------------------------------------------
% rawData format should be: (Num_RX * Num_TX) * Num_verticalScan * Num_horizontalScan * Samples_per_Chirp;
rawDataUniform = reshape(rawDataMonostatic([1:4,9:12],:,:,:),[],Num_horizontalScan,Samples_per_Chirp);
clear rawDataMonostatic


%--------------------------------------------------------------------------
%-- Image Reconstruction Part
%--------------------------------------------------------------------------

%% Reconstruct Image
sarImage = reconstructSARimageFFT_3D(rawDataUniform,frequency,xStepM_mm,lambda_mm/4,-1,zTarget_mm,512);