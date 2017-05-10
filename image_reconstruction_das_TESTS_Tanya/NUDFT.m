function image = NUDFT(sRes, scan, fz, x, fsx)
% fzMid = (max(fz(:)) + min(fz(:)))/2;
% xBW = fx(end) - fx(1);
% fzVec = fz(:);%-fzMid;
% zBW = fz(end,1) - fz(1,1);
omX = 2*pi / fsx * x(:);

Nz = ceil(scan.z_axis(end)/scan.dz);
delta = (1/scan.dz ) / Nz;
minInd = floor(min(fz(:))/delta);
maxInd = ceil(max(fz(:))/delta);
minFz = minInd*delta;
maxFz = maxInd*delta;
fzWanted = (minFz: delta :maxFz)';

zTmp = fz(:) - minFz;
omZ = 2*pi / max(zTmp) * zTmp - pi;
xd = dtft2_adj(sRes(:), [omX, omZ ], single(size(fz,2)), single(numel(fzWanted)),single([0,0]),1).';
timeTake = toc;
fprintf('DFT %.2f seconds\n',timeTake);
image = abs(xd);
% imagesc(db(abs(xd))),colorbar
% save('xd_128x414')
% shg