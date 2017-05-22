function image = InterpFreqBB(Gamma, scan, fx , fz, x, sumForierDomainFlag)
%% 1. Interpolation
% delta = (1/scan.dz ) /scan.Nz;
Nz = ceil(scan.z_axis(end)/scan.dz);
delta = (1/scan.dz ) / Nz;
if sumForierDomainFlag
    minInd = floor(min(fz(:))/delta);
    maxInd = ceil(max(fz(:))/delta);
else 
    minIndGlobal = 92;
    maxIndGlobal = 507;
    minInd = minIndGlobal;
    maxInd = maxIndGlobal;
end
minFz = minInd*delta;
maxFz = maxInd*delta;
fzWanted = (minFz: delta :maxFz)';


Nx = scan.Nx;
deltaX = 1/scan.dx / Nx;
if sumForierDomainFlag
    minIndX = floor(min(x(:))/deltaX);
    maxIndX = ceil(max(x(:))/deltaX);
else
    minIndGlobalX = -118;
    maxIndGlobalX = 118;
    minIndX = minIndGlobalX;
    maxIndX = maxIndGlobalX;
end
minFx = minIndX*deltaX;
maxFx = maxIndX*deltaX;
fxWanted = (minFx: deltaX :maxFx)';
fxAxesWanted = fxWanted; %%%!!!!

% fxAxesWanted = fx;
[fxAxesWantedMesh,fzWantedMesh] = meshgrid(fxAxesWanted, fzWanted);

sResInterp1 = single(griddata(double(x), double(fz), double(Gamma), double(fxAxesWantedMesh),double(fzWantedMesh),'natural'));
sResInterp1(isnan(sResInterp1)) = 0;
if(0)
figure
% mesh(fxAxesWantedMesh, fzWantedMesh, abs(sResInterp1));
mesh(fxAxesWantedMesh, fzWantedMesh, abs(sResInterp1));
title('\Gamma(f,f_x) After Interpolation','fontsize',24)
xlabel('f_x','fontsize',18)
ylabel('f_z','fontsize',18)
shg
end
%% 2. IFFT
image = abs( ifft2(ifftshift(sResInterp1)) ); 
