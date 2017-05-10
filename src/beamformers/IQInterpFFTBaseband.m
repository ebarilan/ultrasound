function image = IQInterpFFTBaseband(scan,dataset,pw_indices)

probeIdx = 38; %38
theta = double(dataset.angles(probeIdx));


fsx = 128/double((dataset.probe_geometry(end,1) - dataset.probe_geometry(1,1)));
fx = (-1/2:1/128:1/2-1/128)*fsx;
f = dataset.modulation_frequency + (-1/2:1/dataset.samples:1/2-1/dataset.samples)*dataset.sampling_frequency;
f = double(f);
[x,y] = meshgrid(fx, f);
S = fftshift(fft2(dataset.data(:,:,probeIdx)));
figure
mesh(x, y, abs(S));
title('S(f,f_x)','fontsize',24)
xlabel('f_x','fontsize',18)
ylabel('f','fontsize',18)


lambda = dataset.c0./f';
tmp1 = bsxfun(@minus, fx, sin(theta)./lambda);
tmp2 = bsxfun(@minus, 1./(lambda.^2), tmp1.^2);
fz = bsxfun(@plus, cos(theta)./lambda, tmp2.^0.5);

tmpNom = bsxfun(@minus, 1./(lambda.^2), fx.^2 ); %%! isin(lambda)
tmpDenom=  -4*pi*  1./(lambda.^2);
constS = bsxfun(@times, tmpNom, 1./tmpDenom);
sRes = S .* constS;

figure
fxRelevant = x;
mesh(fxRelevant, fz, abs(sRes));
title('\Gamma(f,f_x)','fontsize',24)
xlabel('f_x','fontsize',18)
ylabel('f_z','fontsize',18)
%     N1 = 2^6;		% # of time samples
%     n = [0:N1-1]'-N1/2;	% time sample locations
%     dtft2_adj(X, omega, N1, N2, n_shift, useloop)
%     xd = dtft2_adj(Stmp, [x(fIdx,:) fz], N1, 1, [N1/2 0]);
%% Interp1
tic;
% delta = (1/scan.dz ) /scan.Nz;
Nz = ceil(scan.z_axis(end)/scan.dz);
delta = (1/scan.dz ) / Nz;
minInd = floor(min(fz(:))/delta);
maxInd = ceil(max(fz(:))/delta);
minFz = minInd*delta;
maxFz = maxInd*delta;
fzWanted = (minFz: delta :maxFz)';

sResInterp1 = zeros(maxInd + 1, size(fxRelevant,2));
for i = 1:size(fz,2)
    sResInterp1((minInd+1):end,i) = interp1(fz(:,i), sRes(:,i), fzWanted , 'linear',0);
end
sResInterp1(isnan(sResInterp1)) = 0;

fxAxesWanted = fxRelevant(1,:);
fzAxesWanted = (0:delta:maxFz);
[fxAxesWantedMesh,fzWantedMesh] = meshgrid(fxAxesWanted, fzAxesWanted);
figure
mesh(fxAxesWantedMesh, fzWantedMesh, abs(sResInterp1));
title('\Gamma(f,f_x) After Interpolation','fontsize',24)
xlabel('f_x','fontsize',18)
ylabel('f_z','fontsize',18)
shg

% IFFT
imageRecover = abs( ifft2(ifftshift(sResInterp1)) );



%-- compute envelope
% imagRecoverEnvelope = abs( tools.envelope(real(imagRecover)) + 1i * tools.envelope(imag(imagRecover)) );
% imagRecover = tools.envelope(abs(imagRecover));

figure
zAxesFinal = (0:size(imageRecover,1)-1)/(size(imageRecover,1)-1)*scan.z_axis(end);
dxCalc = 1 / (fxAxesWantedMesh(1,end) - fxAxesWantedMesh(1,1));
xAxesFinal = ((0:size(imageRecover,2)-1) - size(imageRecover,2)/2) *dxCalc;
imagesc(xAxesFinal , zAxesFinal , db(imageRecover/max(imageRecover(:))))
title('Final Image Interp1','fontsize',24)
xlabel('X','fontsize',18)
ylabel('Z','fontsize',18)
a = toc; disp(a)


%-- interpolate the requested grid
[xAxesFinalMesh,zAxesFinalMesh] = meshgrid(xAxesFinal,zAxesFinal);
[x_axisMesh , z_axisMesh] = meshgrid(scan.x_axis , scan.z_axis);

resampled_envelope_beamformed_data = zeros(numel(scan.z_axis),numel(scan.x_axis),numel(pw_indices));
for f=1:length(pw_indices)
    resampled_envelope_beamformed_data(:,:,f) = interp2(xAxesFinalMesh,zAxesFinalMesh, imageRecover(:,:,f),x_axisMesh , z_axisMesh);
end


%-- declare an us_image object to store the beamformed data
image = us_image('DAS-RF beamforming');
image.author = 'Alfonso Rodriguez-Molares <alfonso.r.molares@ntnu.no>';
image.affiliation = 'Norwegian University of Science and Technology (NTNU)';
image.algorithm = 'Delay-and-Sum (RF version)';
image.scan = scan;
image.number_plane_waves = cellfun('length',pw_indices);
image.data = resampled_envelope_beamformed_data;
image.transmit_f_number = 0;
% image.receive_f_number = 0; %%rx_f_number;
image.transmit_apodization_window = 'none';
image.receive_apodization_window = 'Tukey 25%';


end