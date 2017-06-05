function [b, Kappa_m,sqrtN] = SPURSInit2(Gamma, scan, fz, x , fx, fsx, sumForierDomainFlag)
% Gamma = conj(Gamma);
b = Gamma(:);
%  b = reshape(Gamma,1,[]).';
%  b = [b(end/2: -1:1) ; b(1:end/2)];
 Nmult = 1;
[c_x c_z] = meshgrid(0:size(x,2)-1,0:size(x,1)-1);
fxx = length(fx)/2;
delta_x = (max(fx) - min(fx)) / numel(fx);
x_correction = single(exp( 1j*2*pi * x / delta_x /2));
% x_correction = single(exp((1j*2*pi()*fxx).*c_x./size(x,2)));
% figure;imagesc(db(real(ifftshift(ifft((x_correction.'))))));colormap('gray');
fzz = 10;%size(fz,1)/2;
Gamma_c = Gamma.*x_correction;
b = Gamma_c(:);

Nz = ceil(scan.z_axis(end)/scan.dz);

% Nz = Nz * Nmult; %%%%%!!!!!
delta_z = (1/scan.dz ) / Nz;
if sumForierDomainFlag
    minInd = floor(min(fz(:))/delta_z);
    maxInd = ceil(max(fz(:))/delta_z);
else 
    minIndGlobal = 92;
    maxIndGlobal = 507;
    minInd = minIndGlobal;
    maxInd = maxIndGlobal;
end

minFz = minInd*delta_z;
maxFz = maxInd*delta_z;
fzWanted = (minFz: delta_z :maxFz)';
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
z_correction = single( exp(1j*2*pi * zTmp / delta_z / 2));
b = b .* z_correction;
% b = b.*exp((1j*2*pi()).*(zTmp*fzz/sqrtN));


xWantedBW = 1/((scan.x_axis(end) - scan.x_axis(1))/sqrtN);
Kappa_x = sqrtN / xWantedBW * x(:); %GOOD; % Kappa_x = sqrtN / fsx * x(:);
% Kappa_x = sqrtN / BW_x_z * x(:); % Kappa_x = 3 * sqrtN / fsx * x(:);

Kappa_m = [Kappa_x , Kappa_z];
% x1 = [(-sqrtN/2:floor(min(Kappa_x))),(ceil(max(Kappa_x)):sqrtN/2)] ;
% [x2 , z2 ]  = meshgrid(x1 , (-sqrtN/2:sqrtN/2));
% b2 = zeros(numel(x2),1);

% Kappa_x = [Kappa_x; x2(:)];
% Kappa_z = [Kappa_z; z2(:)];
% Kappa_m = [Kappa_x , Kappa_z];
% b = [b;b2];
end
