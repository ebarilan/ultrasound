function  [image,imageFFT] = InterpNUFFT2(sRes, scan, fz, x , fsx, sumForierDomainFlag)


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
    icSlice = ic(i);
    idxSlice = (ic == i);
    omZ_Slice = omZ(idxSlice);
    om = double(omZ_Slice);
    Nd = double(single(sqrtN));
    Jd = double(3);
    Kd = double(2*Nd);
    st = nufft_init( om , Nd , Jd , Kd);
    xd = nufft_adj(sRes(idxSlice), st);
    imageFFT(:,i) = fftshift(fft(xd));
end

image = abs(ifft2(ifftshift(imageFFT)));
