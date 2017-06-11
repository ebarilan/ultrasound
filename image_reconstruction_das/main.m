dynamic_range = 60;
addpath(genpath('../src')); % Should probably be modified!!!

% pathToData = 'inputFiles\resSim\';
% pathToData = 'inputFiles\contSim\';
% pathToData = 'inputFiles\resExp\';
pathToData = 'inputFiles\contExp\';


%% DAS
load([pathToData,'imageDAS_struct.mat']);
[~,fig] = imageDAS_struct.show(dynamic_range);


%% NUFFT
load([pathToData,'imageNUFFT.mat']);
imageNUFFT_1D = abs(ifft2(ifftshift(imageNUFFT_1D_FFT)));
[~,fig] = imageNUFFT_1D_struct.show(dynamic_range);

% zero-padding to improve resolution
imageNUFFT_1D_FFT_padded = zeros(size(imageNUFFT_1D_FFT,1),imageNUFFT_1D_struct.scan.Nx);
imageNUFFT_1D_FFT_padded(:,1:size(imageNUFFT_1D_FFT,2)) = imageNUFFT_1D_FFT;
imageNUFFT_1D_padded = abs(ifft2(ifftshift(imageNUFFT_1D_FFT_padded)));

imageNUFFT_1D_struct.data = imageNUFFT_1D_padded;
imageNUFFT_1D_struct.name = 'Embedded. NUFFT 1D,padded'

[~,fig] = imageNUFFT_1D_struct.show(dynamic_range);

% windowing to improve ringing artifact in lateral dimension
apod = hamming(size(imageNUFFT_1D_FFT,2));
apodMat = repmat(apod,1,size(imageNUFFT_1D_FFT,1)).'; 

% figure; imagesc(apodMat);
imageNUFFT_1D_FFT_wind = imageNUFFT_1D_FFT.*apodMat;
imageNUFFT_1D_FFT_padded = zeros(size(imageNUFFT_1D_FFT,1),imageNUFFT_1D_struct.scan.Nx);
imageNUFFT_1D_FFT_padded(:,1:size(imageNUFFT_1D_FFT,2)) = imageNUFFT_1D_FFT_wind;
imageNUFFT_1D_padded = abs(ifft2(ifftshift(imageNUFFT_1D_FFT_padded)));

imageNUFFT_1D_struct.data = imageNUFFT_1D_padded;
imageNUFFT_1D_struct.name = 'Embedded. NUFFT 1D,padded,windowed'

[~,fig] = imageNUFFT_1D_struct.show(dynamic_range);