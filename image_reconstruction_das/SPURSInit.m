function [b, Kappa_m,sqrtN] = SPURSInit(Gamma, scan, fz, x , fx, fsx)
b = Gamma(:);
%  b = reshape(Gamma,1,[]).';
%  b = [b(end/2: -1:1) ; b(1:end/2)];
 Nmult = 1;

Nz = ceil(scan.z_axis(end)/scan.dz);

if(0)
% NNz = scan.Nz;
NNz = ceil(scan.z_axis(end)/scan.dz);
BBWz = 1/scan.dz;
deltaFz = BBWz / NNz;
fzCarrier = (max(fz(:)) + min(fz(:)))/2;
zTmp = fz(:) - fzCarrier;
% Kappa_z = (fz(:)-fzCarrier)/deltaFz;
Kappa_z = (fz(:)-fzCarrier)/deltaFz;

sqrtN = NNz;
% sqrtN = size(fz,1);
% Lx = scan.x(end)-scan.x(1);
% Lz = scan.z(end);
% NmultX = Lz/Lx;
% Kappa_x = NmultX * sqrtN / BBWz * x(:);
xWantedBW = 1/((scan.x_axis(end) - scan.x_axis(1))/sqrtN);
Kappa_z = Nmult * sqrtN / max(zTmp) * zTmp - sqrtN/2;

Kappa_x = sqrtN / xWantedBW * x(:);
Kappa_m = [Kappa_x ,Kappa_z ];
else

% Nz = Nz * Nmult; %%%%%!!!!!
delta = (1/scan.dz ) / Nz;
minInd = floor(min(fz(:))/delta);
maxInd = ceil(max(fz(:))/delta);
minFz = minInd*delta;
maxFz = maxInd*delta;
fzWanted = (minFz: delta :maxFz)';
zTmp = fz(:) - minFz;
sqrtN = numel(fzWanted);
zWantedBW = 1/((scan.z_axis(end) - scan.z_axis(1))/sqrtN);

% fffzz = fz - min(fz(:));
% zzWantedBW = 1/scan.dz;
% sqrtN = scan.Nz;
% Kappa_z = sqrtN/zzWantedBW * fffzz(:) - sqrtN/2;

% Kappa_z = sqrtN / zWantedBW * zTmp - sqrtN/2; % Kappa_z = sqrtN / max(zTmp) * zTmp - sqrtN/2;
% Kappa_z = sqrtN / max(zTmp) * zTmp - sqrtN/2; % Good
% Kappa_z = Nmult * sqrtN / max(zTmp) * zTmp - sqrtN/2;
BW_x_z = max(zTmp);
Kappa_z = sqrtN / BW_x_z * zTmp - sqrtN/2;



xWantedBW = 1/((scan.x_axis(end) - scan.x_axis(1))/sqrtN);
Kappa_x = sqrtN / xWantedBW * x(:); %GOOD; % Kappa_x = sqrtN / fsx * x(:);
% Kappa_x = sqrtN / BW_x_z * x(:); % Kappa_x = 3 * sqrtN / fsx * x(:);

Kappa_m = [Kappa_x , Kappa_z];
x1 = [(-sqrtN/2:floor(min(Kappa_x))),(ceil(max(Kappa_x)):sqrtN/2)] ;
[x2 , z2 ]  = meshgrid(x1 , (-sqrtN/2:sqrtN/2));
b2 = zeros(numel(x2),1);

% Kappa_x = [Kappa_x; x2(:)];
% Kappa_z = [Kappa_z; z2(:)];
% Kappa_m = [Kappa_x , Kappa_z];
% b = [b;b2];

end

end
