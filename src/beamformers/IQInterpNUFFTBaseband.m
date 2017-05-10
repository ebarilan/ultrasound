function image = IQInterpNUFFTBaseband(scan,dataset,pw_indices)

probeIdx = 38; %38
theta = double(dataset.angles(probeIdx));


fsx = 128/double((dataset.probe_geometry(end,1) - dataset.probe_geometry(1,1)));
fx = (-1/2+1/128:1/128:1/2)*fsx;
f = double(dataset.modulation_frequency) + (-1/2:1/dataset.samples:1/2-1/dataset.samples)*double(dataset.sampling_frequency);
f = double(f);
[x,y] = meshgrid(fx, f);
S = fftshift(fft2(double(dataset.data(:,:,probeIdx))));
figure
mesh(x, y, abs(S));
title('S(f,f_x)','fontsize',24)
xlabel('f_x','fontsize',18)
ylabel('f','fontsize',18)


lambda = double(dataset.c0)./f';
tmp1 = bsxfun(@minus, fx, sin(theta)./lambda);
tmp2 = bsxfun(@minus, 1./(lambda.^2), tmp1.^2);
fz = bsxfun(@plus, cos(theta)./lambda, tmp2.^0.5);

tmpNom = bsxfun(@minus, 1./(lambda.^2), fx.^2 ); %%! isin(lambda)
tmpDenom=  -4*pi*  1./(lambda.^2);
constS = bsxfun(@times, tmpNom, 1./tmpDenom);
sRes = S .* constS;

figure
mesh(x, fz, abs(sRes));
title('\Gamma(f,f_x)','fontsize',24)
xlabel('f_x','fontsize',18)
ylabel('f_z','fontsize',18)


if(1)
    Nx = floor(scan.Nx);
    Nz = floor(scan.Nz);
    tic;
    fzMid = (max(fz(:)) + min(fz(:)))/2;
    xBW = x(1,end) - x(1,1);
    fzVec = fz(:)-fzMid;
    zBW = fz(end,1) - fz(1,1);
    omX = 2*pi* x(:); 
%     /xBW;
    omZ = 2*pi*fzVec; 
%     /zBW;
    xd = dtft2_adj(sRes(:), [omX, omZ ], Nx, Nz,[0,0],1);
    timeTake = toc;
    fprintf('DFT %.2f seconds\n',timeTake);
    imagesc(db(abs(xd))),colorbar, shg
end

end