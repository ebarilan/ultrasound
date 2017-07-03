function [imageRecover] = Sampels2NUFourierDomain(dataset,probeIdx,scan, fx_mesh_tmp, fz_mesh_tmp, Gamma_tmp)
%% 1. Transform the signal to Fourier domain [f,fx]

theta = single(dataset.angles(probeIdx));

Nchannels = dataset.channels;

fsx = Nchannels/single((dataset.probe_geometry(end,1) - dataset.probe_geometry(1,1))); % fx = 1/dx

fx =  CreateFreqAxes(Nchannels , fsx); % Create freq axes [-pi,pi]
f = dataset.modulation_frequency + CreateFreqAxes(dataset.samples , dataset.sampling_frequency); % f = modolation + IQ
f = single(f);
[fx_mesh,f_mesh] = meshgrid(fx, f);
figure,scatter(fx_mesh(:), f_mesh(:),[],reshape(abs(fftshift(fft2(dataset.data(:,:,probeIdx)))),[],1))
title('Reguler FFT On S(f_x,f)')
%% 1. Calaulate Gamma's Uniform Grid
%1.1 Fz
Nz = ceil(scan.z_axis(end)/scan.dz);
delta_fz = (1/scan.dz ) / Nz;
minIndGlobal = 92;
maxIndGlobal = 507;
minInd = minIndGlobal;
maxInd = maxIndGlobal;
minFz = minInd*delta_fz;
maxFz = maxInd*delta_fz; %for
fzWanted = (minFz: delta_fz :maxFz)';

%1.2 Fx
delta_fx = fsx/dataset.channels;
minIndGlobal = -117;
maxIndGlobal = 117;
minInd = minIndGlobal;
maxInd = maxIndGlobal;
minFx = minInd*delta_fx;
maxFx = maxInd*delta_fx;
fxWanted = (minFx: delta_fx :maxFx)';
% fxWanted = fx;

[fxGammaUniform, fzGammaUniform] = meshgrid(fxWanted, fzWanted);

%%%%%%%%%%%!!!!!!!!!!!!!!!!!!!!!!!!
fxGammaUniform = fx_mesh_tmp;
fzGammaUniform = fz_mesh_tmp;
%%%%%%%%%%%!!!!!!!!!!!!!!!!!!!!!!!!


%% 2. Calculate S's Nonuniform Grid
BW_f = dataset.sampling_frequency;

fSNonUniform = (dataset.c0/2) * (fxGammaUniform.^2 + fzGammaUniform.^2)./(fxGammaUniform * sin(theta) + fzGammaUniform * cos(theta));

fxSNonUniform = fxGammaUniform;
figure, hold on, scatter(fxSNonUniform(:), fSNonUniform(:));
title('f_x & f')

lambda_tmp2 = single(dataset.c0)./fSNonUniform;
fx0 = sin(theta)./lambda_tmp2;

fxSNonUniformReal = fxSNonUniform - fx0;
figure, scatter(fxSNonUniformReal(:), fSNonUniform(:)), hold on
fSNonUniformBB = fSNonUniform - BW_f;
indCondF = (fSNonUniformBB  >= -BW_f/2) & (fSNonUniformBB <= BW_f/2);
indCondFx = (fxSNonUniformReal  >= -fsx/2) & (fxSNonUniformReal <= fsx/2);
indCond = indCondF & indCondFx;
fSNonUniformCleanBB =  fSNonUniformBB(indCond);
fxSNonUniformClean = fxSNonUniformReal(indCond);

scatter(fxSNonUniformClean(:), fSNonUniformCleanBB(:) + BW_f);
hold off

s = dataset.data(:,:,probeIdx);
omX = fxSNonUniformClean   * (2*pi/fsx);
omF = fSNonUniformCleanBB  * (2*pi/BW_f);
om = [omF(:), omX(:)];

if 0 %NU-DFT
    n_shift = [0 0];    
    useloop = 0;
    S_NU = dtft2(s, om, n_shift, useloop);
else
    
    Nd = double(size(s));
    Jd = double([3,3]);
    Kd = 2*Nd;%double(2*[size(om,1),size(om,1)]);
    
    st = nufft_init( om , Nd , Jd , Kd);
    xd = nufft(s, st);
    S_NU = xd;
    
    
    if 0
        st2 = nufft_init( om , Nd , Jd , Kd);
        xd2 = nufft_adj(S_NU(:), st2);
        image2 = abs(xd2);
        figure,
        subplot(1,2,1)
        imagesc(abs(double(s)))
        title('Original s')
        
        subplot(1,2,2)
        imagesc(abs(double(image2)))
        title('s \rightarrow NUFFT \rightarrow NUFFT adjoint \rightarrow s')
    end
end
figure
scatter(fxSNonUniformClean, fSNonUniformCleanBB+BW_f, [] , abs(S_NU))
title('S(f_x,f) Non-Uniform')

f_mesh = fSNonUniformCleanBB+BW_f;

lambda_tmp2 = single(dataset.c0)./f_mesh;
fx0 = sin(theta)./lambda_tmp2;

fx_mesh = fxSNonUniformClean + fx0;

tmp1 = fx_mesh - f_mesh * sin(theta) / dataset.c0;
tmp2 = ((f_mesh / dataset.c0).^2 - tmp1.^2).^0.5;
fz_mesh = f_mesh * cos(theta) / dataset.c0 + tmp2;

fxminusfx0 = fx_mesh - fx0;
tmpNom = ( (f_mesh / dataset.c0).^2  -  fxminusfx0.^2 ).^0.5;
tmpDenom = -4*pi* (f_mesh / dataset.c0).^2;
constS = tmpNom ./ tmpDenom;
Gamma_Uniform = S_NU .* constS;


figure,
scatter(fx_mesh,fz_mesh,[],abs(Gamma_Uniform), 'filled')
title('\Gamma(f_x,f_z) Uniform')
xlabel('f_x','fontsize',18)
ylabel('f_z','fontsize',18)
grid on
xlim([-4e3,4e3])

if 1
    image3 = InterpNUFFT(Gamma_Uniform, scan, fz_mesh, fx_mesh , fsx, 0);
    figure, imagesc(db(abs(image3)))
end

NzGamma = round((max(fz_mesh)-min(fz_mesh))/delta_fz)+1;

NxGamma = round((max(fx_mesh)-min(fx_mesh))/delta_fx)+1;

%     GammaFinal = zeros(NxGamma, NzGamma);
%
%     indX = round((fx_mesh-min(fx_mesh))/delta_fx) + 1;
%     indZ = round((fz_mesh-min(fz_mesh))/delta_fz) + 1;
%
%
%     for i = 1:numel(indX)
%         GammaFinal(indX(i),indZ(i)) = Gamma_Uniform(i);
%     end
%
%     imageRecover = abs(ifft2(ifftshift(GammaFinal.')));
%     figure, imagesc(db(imageRecover))
GammaFinal = zeros(numel(fxWanted), numel(fzWanted));

indX = round((fx_mesh-minFx)/delta_fx) + 1;
indZ = round((fz_mesh-minFz)/delta_fz) + 1;


for i = 1:numel(indX)
    GammaFinal(indX(i),indZ(i)) = Gamma_Uniform(i);
end

imageRecover = abs(ifft2(ifftshift(GammaFinal.')));
figure, imagesc(db(imageRecover))


imageRecover = image3;%%%%%%%%%%%%%%!!!!!!!!!!!!!!!!!

end