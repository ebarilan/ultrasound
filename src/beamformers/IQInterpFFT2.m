function image = IQInterpFFT2(scan,dataset,pw_indices, processType, pulseFlag)
pulseFlag = 0;

%% 1. Transform the signal to Fourier domain [f,fx]
probeIdx = 38;
theta = single(dataset.angles(probeIdx));

Nchannels = dataset.channels;

fsx = Nchannels/single((dataset.probe_geometry(end,1) - dataset.probe_geometry(1,1))); % fx = 1/dx

% if processType == 7
%     Nchannels = scan.Nx+1;
%     fx =  CreateFreqAxes(Nchannels , fsx); % Create freq axes [-pi,pi]
% else
    fx =  CreateFreqAxes(Nchannels , fsx); % Create freq axes [-pi,pi]
% end

% if processType == 8
%     f = dataset.modulation_frequency + CreateFreqAxes(scan.Nz*10 , dataset.sampling_frequency); % f = modolation + IQ
% else
    f = dataset.modulation_frequency + CreateFreqAxes(128 , dataset.sampling_frequency); % f = modolation + IQ
% end

f = single(f);

[x,y] = meshgrid(fx, f);

% if processType == 7
%     data_temp = dataset.data(:,:,probeIdx);
%     data_padded = zeros(numel(f),Nchannels);
% %     data_padded(:,1:size(data_temp,2)) = data_temp;
%     ind_help1 = Nchannels/2;
%     ind_help2 = size(data_temp,2)/2;
%     ind_help3 = ind_help1-ind_help2+1;
%     data_padded(:,ind_help3:ind_help3+size(data_temp,2)-1) = data_temp;
%     S = fftshift(fft2(data_padded));
% elseif processType == 8
%     data_temp = dataset.data(:,:,probeIdx);
%     data_padded = zeros(numel(f),Nchannels);
%     data_padded(1:size(data_temp,1),:) = data_temp;
%     S = fftshift(fft2(data_padded));
% else
    S = fftshift(fft2(dataset.data((1:128),:,probeIdx)));
% end

S = single(S);
if(1)
figure
mesh(x, y, abs(S));
title('S(f,f_x)','fontsize',24)
xlabel('f_x','fontsize',18)
ylabel('f','fontsize',18)
end
%% Introducing pulse shape
if pulseFlag
  f0 = dataset.modulation_frequency;
  fBW = 0.67 ;
  pulseSpectrum = getPulseSpectrum(f,f0,fBW).';
else
  pulseSpectrum = ones(size(f))';  
end

%% 2. Transform from [f,fx] -> [fz,fx]
lambda = single(dataset.c0)./f';
tmp1 = bsxfun(@minus, fx, sin(theta)./lambda);
tmp2 = bsxfun(@minus, 1./(lambda.^2), tmp1.^2);
fz = bsxfun(@plus, cos(theta)./lambda, tmp2.^0.5);

tmpNom = bsxfun(@minus, 1./(lambda.^2), fx.^2 ); %%! isin(lambda)
tmpDenom=  -4*pi*  1./(lambda.^2).*pulseSpectrum;
constS = bsxfun(@times, tmpNom, 1./tmpDenom);
Gamma = S .* constS;

if(0)
figure
mesh(x, fz, abs(Gamma));
title('\Gamma(f,f_x)','fontsize',24)
xlabel('f_x','fontsize',18)
ylabel('f_z','fontsize',18)
end
%% 3. Transform back to image.
%cases:
%       3 - Interpolate the all frequency map. IFFT. Hilbert.
%       4 - Interpolate the signal in baseband. IFFT.
%       5 - Using direct iDFT. Accurate transform. (Time issue)
%       6 - Using NUFFT (Jeff Fessler).
tic
switch processType
    case 3
        imageRecover = single(InterpFullFreq(Gamma, scan, fx , fz, x));
    case 4
        imageRecover = InterpFreqBB(Gamma, scan, fx , fz, x);
    case 5
        imageRecover = NUDFT(Gamma, scan, fz, x , fsx);
    case {6,8}
        imageRecover = single(InterpNUFFT2(Gamma, scan, fz, x , fsx));
%     case 7
%         imageRecover = InterpFreqBB(Gamma, scan, fx , fz, x);
%     case 8
%         imageRecover = InterpFreqBB(Gamma, scan, fx , fz, x);
    case 7
        profile on;
        [b, Kappa_m,sqrtN] = SPURSInit(Gamma, scan, fz, x , fx, fsx);
        BsplineDegree = 3;
        Rho = 1e-3;
        Niterations = 2;
        OverGridFactor = 2;
        FilterInImageSpace = 1;
        
        SPURS_settings.sqrtN = sqrtN;
        SPURS_settings.KernelFunctionString = 'Bspline';
        SPURS_settings.KernelFunctionDegree = BsplineDegree;
        SPURS_settings.ReusePrecalculatedData = 1;
        SPURS_settings.Rho = Rho;
        SPURS_settings.Niterations = Niterations;
        SPURS_settings.UseW = 0;
        SPURS_settings.ForceGenrateNewPhi = 0;
        SPURS_settings.ForceFactorPsi = 0;
        SPURS_settings.SavePSI = 0;
        SPURS_settings.OverGridFactor = OverGridFactor;
        SPURS_settings.alpha = 1;
        SPURS_settings.CalcOptimalAlpha = 1;
        SPURS_settings.FilterInImageSpace = FilterInImageSpace;
        [ OutputImages , b_hat] = SPURS(double(b), double(Kappa_m), SPURS_settings);
        imageRecover = fftshift(OutputImages(:,:,end));
        profile viewer;
end
toc



%% 4. Image Interpolation
zAxesFinal = single((0:size(imageRecover,1)-1)/(size(imageRecover,1)-1)*scan.z_axis(end));
dxCalc = 1 / (fx(end) - fx(1));
xAxesFinal = ((0:size(imageRecover,2)-1) - size(imageRecover,2)/2) *dxCalc;
% figure
% imagesc(xAxesFinal , zAxesFinal , db(imageRecover/max(imageRecover(:))))
% title('Final Image Interp1','fontsize',24)
% xlabel('X','fontsize',18)
% ylabel('Z','fontsize',18)

%-- interpolate the requested grid
[xAxesFinalMesh,zAxesFinalMesh] = meshgrid(single(xAxesFinal),single(zAxesFinal));
[x_axisMesh , z_axisMesh] = meshgrid(scan.x_axis , scan.z_axis);

resampled_envelope_beamformed_data = zeros(numel(scan.z_axis),numel(scan.x_axis),numel(pw_indices));
for f=1:length(pw_indices)
    resampled_envelope_beamformed_data(:,:,f) = interp2(xAxesFinalMesh,zAxesFinalMesh, imageRecover(:,:,f),x_axisMesh , z_axisMesh);
end

No_resampled_envelope_beamformed_data = imageRecover((size(imageRecover,1)-scan.Nz+1):end,:);

%-- declare an us_image object to store the beamformed data
image = us_image('DAS-RF beamforming');
image.author = 'Alfonso Rodriguez-Molares <alfonso.r.molares@ntnu.no>';
image.affiliation = 'Norwegian University of Science and Technology (NTNU)';
image.algorithm = 'Delay-and-Sum (RF version)';
image.scan = scan;
image.number_plane_waves = cellfun('length',pw_indices);
image.data = No_resampled_envelope_beamformed_data; % image.data = resampled_envelope_beamformed_data;
image.transmit_f_number = 0;
% image.receive_f_number = 0; %%rx_f_number;
image.transmit_apodization_window = 'none';
image.receive_apodization_window = 'Tukey 25%';


end