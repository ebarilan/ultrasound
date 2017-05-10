function image = InterpFullFreq(sRes, scan, fx , fz, x)
%% 1. Interpolation
% delta = (1/scan.dz ) /scan.Nz;
Nz = ceil(scan.z_axis(end)/scan.dz);
delta = (1/scan.dz ) / Nz;
minInd = floor(min(fz(:))/delta);
maxInd = ceil(max(fz(:))/delta);
minFz = minInd*delta;
maxFz = maxInd*delta;
fzWanted = (minFz: delta :maxFz)';

sResInterp1 = zeros(maxInd + 1, size(x,2));
for i = 1:size(fz,2)
    sResInterp1((minInd+1):end,i) = interp1(fz(:,i), sRes(:,i), fzWanted , 'linear',0);
end
sResInterp1(isnan(sResInterp1)) = 0;

%% 2. Duplicate the negative frequency
fxAxesWanted = fx;
fzAxesWanted = (0:delta:maxFz);
[fxAxesWantedMesh,fzWantedMesh] = meshgrid(fxAxesWanted, fzAxesWanted);
figure
mesh(fxAxesWantedMesh, fzWantedMesh, abs(sResInterp1));
title('\Gamma(f,f_x) After Interpolation','fontsize',24)
xlabel('f_x','fontsize',18)
ylabel('f_z','fontsize',18)
shg

fzAxesWantedExt = (-maxFz:delta:maxFz);
midFzAxes = floor(numel(fzAxesWantedExt)/2);
[fxAxesWantedMeshExt,fzWantedMeshExt] = meshgrid(fxAxesWanted, fzAxesWantedExt);
sResInterp1Ext = zeros( size(fxAxesWantedMeshExt));
% negative reflecation
sResInterp1Ext( 1:midFzAxes , 1:end) = conj( sResInterp1(end:-1:2 , [1,(end:-1:2)] ));
% positive
sResInterp1Ext( (midFzAxes+1):end , 1:end) = sResInterp1(1:end ,1:end);

figure
mesh(fxAxesWantedMeshExt, fzWantedMeshExt, abs(sResInterp1Ext));
title('\Gamma(f,f_x) After Interpolation & Extend','fontsize',24)
xlabel('f_x','fontsize',18)
ylabel('f_z','fontsize',18)
shg

%% 3.IFFT
imagRecover = ifft2(ifftshift(sResInterp1Ext));



%% 4. Hilbert - compute envelope
image = tools.envelope(imagRecover);
% imagRecoverEnvelope = tools.envelope(abs(imagRecover));
