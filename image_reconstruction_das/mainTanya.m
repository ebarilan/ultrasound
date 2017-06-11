
dynamic_range = 60;

ImageOpt = 2; %[1,2] - resolution, contrast
numAngels = 75;
acquisition_type = 1; %[1,2] - simulation, experiment

IsEmbedded = 1;


%% DAS
% ProccesOpt = 1; %[1, 4, 6];
% phantom_type = ImageOpt;
% data_type = ProccesOpt;
% imageDAS_struct = imageRecovery(acquisition_type , phantom_type , data_type, IsEmbedded, numAngels,[]);
% [~,fig] = imageDAS_struct.show(dynamic_range);


%% NUFFT 1D
ProccesOpt = 8; %[1, 4, 6];
phantom_type = ImageOpt;
data_type = ProccesOpt;

[imageNUFFT_1D_struct,imageNUFFT_1D_FFT] = imageRecovery(acquisition_type , phantom_type , data_type, IsEmbedded, numAngels,[]);
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


%% SPURS
% ProccesOpt = 7;
% phantom_type = ImageOpt;
% data_type = ProccesOpt;
% 
% spursConfig.Niterations = 1;
% spursConfig.OverGridFactor = 1;
% spursConfig.KernelFunctionDegree = 3;
% spursConfig.Rho = 1e-3;
% spursConfig.FilterInImageSpace = 1;
% spursConfig.SavePSI = 0;
% 
% [imageSPURS_2D_struct,imageSPURS_2D_FFT] = imageRecovery(acquisition_type , phantom_type , data_type, IsEmbedded, numAngels,spursConfig);
% imageNUFFT_1D = abs(ifft2(ifftshift(imageNUFFT_1D_FFT)));
% [~,fig] = imageNUFFT_1D_struct.show(dynamic_range);


