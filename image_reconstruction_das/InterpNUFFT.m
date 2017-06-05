function image = InterpNUFFT(sRes, scan, fz, x , fsx)


Nz = ceil(scan.z_axis(end)/scan.dz);
delta = (1/scan.dz ) / Nz;
minInd = floor(min(fz(:))/delta);
maxInd = ceil(max(fz(:))/delta);
minFz = minInd*delta;
maxFz = maxInd*delta;
fzWanted = (minFz: delta :maxFz)';

sqrtN = numel(fzWanted);

zTmp = fz(:) - minFz;
% BW_z = 1/scan.dz;
BW_z = max(zTmp);
omZ = 2*pi /BW_z * zTmp - pi; % omZ = 2*pi / max(zTmp) * zTmp - pi;

BW_x = 1/scan.dx;
xWantedBW = 1/((scan.x_axis(end) - scan.x_axis(1))/sqrtN);
% omX = 2*pi / BW_x * x(:);
omX = 2*pi / xWantedBW * x(:); % omX = 2*pi / fsx * x(:);



om = double([omX, omZ ]);
Nd = double([single(sqrtN), single(sqrtN)]); 
% Nd = double([single(scan.Nx), single(Nz)]); 
Jd = double([3,3]);
Kd = double(2*Nd);
%%%%!!!!!
% om = double(Kappa_m)/sqrtN*2*pi;
% Nd = [sqrtN, sqrtN];
% Kd = double(2*Nd);
% %%%%%!!!!
st = nufft_init( om , Nd , Jd , Kd);
xd = nufft_adj(sRes(:), st);
image = abs(xd)';
