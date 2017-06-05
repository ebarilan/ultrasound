function [fx_mesh, fz_mesh, Gamma, fx, fsx] = Sampels2FourierDomainTanya(dataset,probeIdx)
%% 1. Transform the signal to Fourier domain [f,fx]

theta = single(dataset.angles(probeIdx));

Nchannels = dataset.channels;

fsx = Nchannels/single((dataset.probe_geometry(end,1) - dataset.probe_geometry(1,1))); % fx = 1/dx

fx =  CreateFreqAxes(Nchannels , fsx); % Create freq axes [-pi,pi]
f = dataset.modulation_frequency + CreateFreqAxes(dataset.samples , dataset.sampling_frequency); % f = modolation + IQ
f = single(f);
[fx_mesh,f_mesh] = meshgrid(fx, f);

S_no_wrap = fftshift(fft2(dataset.data(:,:,probeIdx)));
fx_mesh_no_wrap = fx_mesh;

if(0)
    figure
    scatter(fx_mesh(:), f_mesh(:), [], abs(S_no_wrap(:)), 'filled');
    title({'S(f,f_x)'; sprintf('Angle %.2f',theta)},'fontsize',16)
    xlabel('f_x','fontsize',18)
    ylabel('f','fontsize',18)
    xlim([-4e3,4e3])
end

lambda = single(dataset.c0)./f';
S_f_x = fftshift( fft(dataset.data(:,:,probeIdx),[],1) ,1);
fx0 = sin(theta)./lambda;
xGeo = linspace( dataset.probe_geometry(1,1), dataset.probe_geometry(end,1), size(dataset.probe_geometry,1));
expPhase = exp(1i*2*pi * fx0 * xGeo );
S = fftshift( fft(S_f_x.*expPhase,[],2) , 2);
S = single(S);


if theta < 0
    unwrapLine = fsx/2 + fx0;
    indWrap = (fx_mesh > unwrapLine);
    tmp = zeros(size(fx_mesh,1),1);
    for iRow = 1:size(fx_mesh,1)
        tmp(iRow) = fx_mesh(iRow, find(indWrap(iRow,:) == 1,1,'first') );
    end
    delta_tanya = tmp - unwrapLine;
    fx_mesh(indWrap) = fx_mesh(indWrap) - fsx;
    fx_mesh = bsxfun(@minus, fx_mesh, delta_tanya);
elseif theta > 0
    unwrapLine = -fsx/2 + fx0;
    indWrap = (fx_mesh < unwrapLine);    
    tmp = zeros(size(fx_mesh,1),1);
    for iRow = 1:size(fx_mesh,1)
        tmp(iRow) = fx_mesh(iRow, find(indWrap(iRow,:) == 1,1,'last') );
    end
    delta_tanya = tmp - unwrapLine;
    fx_mesh(indWrap) = fx_mesh(indWrap) + fsx;
    fx_mesh = bsxfun(@minus, fx_mesh, delta_tanya);    
end





%%%% Debug the shift on fx axis%%%%%%%%%
fx_mesh_Exact = zeros(dataset.samples,Nchannels);

for ii = 1:dataset.samples
    fx_mesh_Exact(ii,:) = (-fsx/2+fx0(ii):fsx/Nchannels:fsx/2+fx0(ii)-fsx/(Nchannels+1));
end

figure;
scatter(fx_mesh(:),f_mesh(:),[],(1:numel(fx_mesh))); colorbar
hold on; 
scatter(fx_mesh_Exact(:),f_mesh(:),[],(1:numel(fx_mesh_Exact))); colorbar
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fx_mesh = fx_mesh_Exact;%%%%%%%%%%%%!!!!!!!!!!!!!!!
if(0)
    figure
    %     subplot(1,2,1)
    scatter(fx_mesh_no_wrap(:), f_mesh(:), [], abs(S(:)), 'filled');
    title({'S(f,f_x-f_{x0}) Shifted'; sprintf('Angle %.2f',theta)},'fontsize',16)
    xlabel('f_x','fontsize',18)
    ylabel('f','fontsize',18)
    xlim([-4e3,4e3])
    hold on;
    
    scatter(unwrapLine, f,[], 'd', 'filled')
    hold off;
    
    %     subplot(1,2,2)
    figure
    scatter(fx_mesh(:), f_mesh(:),[], abs(S(:)), 'filled');
    title({'S(f,f_x-f_{x0}) Shifted. Spectral Unwrapping'; sprintf('Angle %.2f',theta)},'fontsize',16)
    xlabel('f_x','fontsize',18)
    ylabel('f','fontsize',18)
    xlim([-4e3,4e3])
    hold on;
    
    scatter(unwrapLine, f,[], 'd', 'filled')
    hold off;
end




%% 2. Transform from [f,fx] -> [fz,fx]

%2.1 fz Calculation

tmp1 = fx_mesh - f_mesh * sin(theta) / dataset.c0;
tmp2 = ((f_mesh / dataset.c0).^2 - tmp1.^2).^0.5;
fz_mesh = f_mesh * cos(theta) / dataset.c0 + tmp2;

% define the pulse shape
fracBW = 0.67;
B = dataset.modulation_frequency*fracBW;
a = (pi*B/2)^2/log(2);
transducerFreqResponseTheor = 0.5*sqrt(pi/a)*(exp(-pi^2/a*(f_mesh-dataset.modulation_frequency).^2)); % Frequency response  - theoretical
B = transducerFreqResponseTheor/max(transducerFreqResponseTheor(:));


fxminusfx0 = bsxfun(@minus, fx_mesh ,fx0);
tmpNom = ( (f_mesh / dataset.c0).^2  -  fxminusfx0.^2 ).^0.5;
tmpDenom = -4*pi*B.*(f_mesh / dataset.c0).^2;
constS = tmpNom ./ tmpDenom;
Gamma = S .* constS;

if(0)
    figure
    scatter(fx_mesh(:), fz_mesh(:),[], abs(Gamma(:)),'filled')
    title({'\Gamma(f_x,f_z)'; sprintf('Angle %.2f',theta)},'fontsize',16,'Interpreter','tex')
    xlabel('f_x','fontsize',18)
    ylabel('f_z','fontsize',18)
    grid on
    xlim([-4e3,4e3])
    
end

end