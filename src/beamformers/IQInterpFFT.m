function image = IQInterpFFT(scan,dataset,pw_indices, processType, pulseFlag)
pulseFlag = 0;

%% 1. Transform the signal to Fourier domain [f,fx]
probeIdx = pw_indices{1}(1);
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
    f = dataset.modulation_frequency + CreateFreqAxes(dataset.samples , dataset.sampling_frequency); % f = modolation + IQ
% end

f = single(f);

[fx_mesh,f_mesh] = meshgrid(fx, f);

S_no_wrap = fftshift(fft2(dataset.data(:,:,probeIdx)));
fx_mesh_no_wrap = fx_mesh;

lambda = single(dataset.c0)./f';
S_f_x = fftshift( fft(dataset.data(:,:,probeIdx),[],1) ,1);
fx0 = sin(theta)./lambda;
xGeo = linspace( dataset.probe_geometry(1,1), dataset.probe_geometry(end,1), size(dataset.probe_geometry,1));
expPhase = exp(1i*2*pi * fx0 * xGeo );
S = fftshift( fft(S_f_x.*expPhase,[],2) , 2);

if theta <= 0
    unwrapLine = fx(end) + fx0;
    indWrap = (fx_mesh > unwrapLine);
    fx_mesh(indWrap) = fx_mesh(indWrap) - fsx;
else
    unwrapLine = fx(1) + fx0;
    indWrap = (fx_mesh < unwrapLine);
    fx_mesh(indWrap) = fx_mesh(indWrap) + fsx;
end


S = single(S);
if(1)
figure
subplot(1,2,1)
mesh(fx_mesh_no_wrap, f_mesh, abs(S));
title('S(f,f_x) No Wrapping','fontsize',24)
xlabel('f_x','fontsize',18)
ylabel('f','fontsize',18)
hold on;

scatter(unwrapLine, f)
hold off;

subplot(1,2,2)
mesh(fx_mesh, f_mesh, abs(S));
title('S(f,f_x) With Wrapping','fontsize',24)
xlabel('f_x','fontsize',18)
ylabel('f','fontsize',18)
hold on;

scatter(unwrapLine, f)
hold off;



end
% %% Introducing pulse shape
% if pulseFlag
%   f0 = dataset.modulation_frequency;
%   fBW = 0.67 ;
%   pulseSpectrum = getPulseSpectrum(f,f0,fBW).';
% else
%   pulseSpectrum = ones(size(f))';  
% end

if(0) %%%%%!!!!
    indS = size(S,1)/2-size(S,2)/2+ (1:size(S,2));
    S = S(indS,:);
    f = f(indS);
    x = x(indS,:);
end
%% 2. Transform from [f,fx] -> [fz,fx]

%2.1 fz Calculation

tmp1 = fx_mesh - f_mesh * sin(theta) / dataset.c0;
tmp2 = ((f_mesh / dataset.c0).^2 - tmp1.^2).^0.5;
fz_mesh = f_mesh * cos(theta) / dataset.c0 + tmp2;

fxminusfx0 = bsxfun(@minus, fx_mesh ,fx0);
tmpNom = ( (f_mesh / dataset.c0).^2  -  fxminusfx0.^2 ).^0.5;
tmpDenom = -4*pi* (f_mesh / dataset.c0).^2;
constS = tmpNom ./ tmpDenom;
Gamma = S .* constS;

if(0)
tmp1 = bsxfun(@minus, fx, sin(theta)./lambda);
tmp2 = bsxfun(@minus, 1./(lambda.^2), tmp1.^2);

fz_mesh = bsxfun(@plus, cos(theta)./lambda, tmp2.^0.5);
    
    
fxminusfx0 = bsxfun(@minus, fx, fx0 );
tmpNom = bsxfun(@minus, 1./(lambda.^2), fxminusfx0.^2 ); %%! isin(lambda)
tmpDenom=  -4*pi*  1./(lambda.^2);%.*pulseSpectrum;
constS = bsxfun(@times, sqrt(tmpNom), 1./tmpDenom);
Gamma = S .* constS;

k = 2*pi*f'/single(dataset.c0);
kx = 2*pi*fx;
tmpNom2 = bsxfun(@minus, k.^2 ,kx.^2);
tmpDenom2 = -2*k.^2;
constS2 = bsxfun( @times, sqrt( tmpNom2) ,1./tmpDenom2) ;
Gamma2 = S .* constS2; 
x_plus_x0 = bsxfun( @plus, fx_mesh, fx0);
end



if(1)
figure
mesh(fx_mesh, fz_mesh, abs(Gamma));
title('\Gamma(f_z,f_x) Reg','fontsize',24)
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
        imageRecover = single(InterpFullFreq(Gamma, scan, fx , fz_mesh, fx_mesh));
    case 4
        imageRecover = InterpFreqBB(Gamma, scan, fx , fz_mesh, fx_mesh);
    case 5
        imageRecover = NUDFT(Gamma, scan, fz_mesh, fx_mesh , fsx);
    case 6
        [b, Kappa_m,sqrtN] = SPURSInit(Gamma, scan, fz_mesh, fx_mesh , fx, fsx);
        imageRecover = single(InterpNUFFT(Gamma, scan, fz_mesh, fx_mesh , fsx,Kappa_m,sqrtN));
%     case 7
%         imageRecover = InterpFreqBB(Gamma, scan, fx , fz, x);
%     case 8
%         imageRecover = InterpFreqBB(Gamma, scan, fx , fz, x);
    case 7
%         profile on;
        [b, Kappa_m,sqrtN] = SPURSInit(Gamma, scan, fz_mesh, fx_mesh , fx, fsx);
        BsplineDegree = 3;
        Rho = 1e-3;%%
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
%         [ OutputImages_real  , b_hat] = SPURS(double(real(b)), double(Kappa_m), SPURS_settings);
%         [ OutputImages_imag  , b_hat] = SPURS(double(imag(b)), double(Kappa_m), SPURS_settings);
        imageRecover = OutputImages(:,:,end);
%         profile viewer;
end
toc



%% 4. Image Interpolation
zAxesFinal = single((0:size(imageRecover,1)-1)/(size(imageRecover,1)-1)*scan.z_axis(end));
% dxCalc = 1 /(2*(fx(end) - fx(1)));
% if processType==7
    xAxesFinal = ((0:size(imageRecover,2)-1) - size(imageRecover,2)/2)/size(imageRecover,2)*2*scan.x_axis(end);
% else
%     xAxesFinal = ((0:size(imageRecover,2)-1) - size(imageRecover,2)/2) *dxCalc;
% end

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

% No_resampled_envelope_beamformed_data = imageRecover(((size(imageRecover,1)-scan.Nz+1):end),:);

%-- declare an us_image object to store the beamformed data
image = us_image('DAS-RF beamforming');
image.author = 'Alfonso Rodriguez-Molares <alfonso.r.molares@ntnu.no>';
image.affiliation = 'Norwegian University of Science and Technology (NTNU)';
image.algorithm = 'Delay-and-Sum (RF version)';
image.scan = scan;
image.number_plane_waves = cellfun('length',pw_indices);
image.data = resampled_envelope_beamformed_data; % image.data = resampled_envelope_beamformed_data;
image.transmit_f_number = 0;
% image.receive_f_number = 0; %%rx_f_number;
image.transmit_apodization_window = 'none';
image.receive_apodization_window = 'Tukey 25%';


end