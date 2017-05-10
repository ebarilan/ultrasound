function image = InterpFreqBB(sRes, scan, fx , fz, x )
%% 1. Interpolation
% delta = (1/scan.dz ) /scan.Nz;
Nz = ceil(scan.z_axis(end)/scan.dz);
delta = (1/scan.dz ) / Nz;
minInd = floor(min(fz(:))/delta);
maxInd = ceil(max(fz(:))/delta);
minFz = minInd*delta;
maxFz = maxInd*delta;
fzWanted = (minFz: delta :maxFz)';

fxAxesWanted = fx;
% fzAxesWanted = (0:delta:maxFz);
[fxAxesWantedMesh,fzWantedMesh] = meshgrid(fxAxesWanted, fzWanted);

sResInterp1 = single(griddata(double(x), double(fz), double(sRes), double(fxAxesWantedMesh),double(fzWantedMesh),'cubic'));
sResInterp1(isnan(sResInterp1)) = 0;
figure
% mesh(fxAxesWantedMesh, fzWantedMesh, abs(sResInterp1));
mesh(fxAxesWantedMesh, fzWantedMesh, abs(sResInterp1));
title('\Gamma(f,f_x) After Interpolation','fontsize',24)
xlabel('f_x','fontsize',18)
ylabel('f_z','fontsize',18)
shg

%% 2. IFFT
image = abs( ifft2(ifftshift(sResInterp1)) ); 
