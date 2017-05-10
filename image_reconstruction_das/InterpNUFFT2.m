function image = InterpNUFFT2(sRes, scan, fz, x , fsx)
BW_x = 1/scan.dx;
omX = 2*pi / BW_x * x(:); % omX = 2*pi / fsx * x(:);

Nz = ceil(scan.z_axis(end)/scan.dz);
delta = (1/scan.dz ) / Nz;
minInd = floor(min(fz(:))/delta);
maxInd = ceil(max(fz(:))/delta);
minFz = minInd*delta;
maxFz = maxInd*delta;
fzWanted = (minFz: delta :maxFz)';

zTmp = fz(:) - minFz;
BW_z = BW_x;
omZ = 2*pi /BW_z * zTmp - pi; % omZ = 2*pi / max(zTmp) * zTmp - pi;
om = double([omX, omZ ]);
Nd = double([single(scan.Nx), single(scan.Nx)]); % Nd = double([single(size(fz,2)), single(numel(fzWanted))]);
Jd = double([3,3]);
Kd = double(2*Nd);
st = nufft_init( om , Nd , Jd , Kd);
xd = nufft_adj(sRes(:), st);
image = abs(xd)';
