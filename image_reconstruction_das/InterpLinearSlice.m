function image = InterpLinearSlice(sRes, scan, fz, x , fsx, sumForierDomainFlag)


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

sqrtN = numel(fzWanted);

zTmp = fz(:) - minFz;
% BW_z = 1/scan.dz;
BW_z = max(zTmp);
omZ = 2*pi /BW_z * zTmp - pi; % omZ = 2*pi / max(zTmp) * zTmp - pi;



[C,ia,ic] = unique(x);
imageFFT = zeros(sqrtN, numel(C));
for i = 1:numel(C)
    idxSlice = (ic == i);
    omZ_Slice = omZ(idxSlice);
    om = double(omZ_Slice);
    xq = linspace(-pi,pi,sqrtN);
    xd = interp1(om, sRes(idxSlice), xq, 'linear',0);
    imageFFT(:,i) = xd;
end

image = abs(ifft2(ifftshift(imageFFT)));