function [fx_mesh, fz_mesh, Gamma, fx, fsx] = Sampels2FourierDomain(dataset,probeIdx)
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

lambda = single(dataset.c0)./f';
S_f_x = fftshift( fft(dataset.data(:,:,probeIdx),[],1) ,1);
fx0 = sin(theta)./lambda;
xGeo = linspace( dataset.probe_geometry(1,1), dataset.probe_geometry(end,1), size(dataset.probe_geometry,1));
expPhase = exp(1i*2*pi * fx0 * xGeo );
S = fftshift( fft(S_f_x.*expPhase,[],2) , 2);
S = single(S);

if theta <= 0
    unwrapLine = fsx/2 + fx0;
    unwrap_mat = repmat(unwrapLine,1,size(fx_mesh,2));
    indWrap = (fx_mesh > unwrap_mat);
    fx_mesh(indWrap) = fx_mesh(indWrap) - fsx;
else
    unwrapLine = -fsx/2 + fx0;
    unwrap_mat = repmat(unwrapLine,1,size(fx_mesh,2));
    indWrap = (fx_mesh < unwrap_mat);
    fx_mesh(indWrap) = fx_mesh(indWrap) + fsx;
end

if(0)
    figure
    subplot(1,2,1)
    mesh(fx_mesh_no_wrap, f_mesh, abs(S));
    title('S(f,f_x) No Wrapping','fontsize',24)
    xlabel('f_x','fontsize',18)
    ylabel('f','fontsize',18)
    hold on;
    
    scatter(unwrapLine, f)
    hold off;
    
    subplot(1,2,2)
    mesh(fx_mesh, f_mesh, abs(S));
    title('S(f,f_x) With Wrapping','fontsize',24)
    xlabel('f_x','fontsize',18)
    ylabel('f','fontsize',18)
    hold on;
    
    scatter(unwrapLine, f)
    hold off;
end




%% 2. Transform from [f,fx] -> [fz,fx]

%2.1 fz Calculation

tmp1 = fx_mesh - f_mesh * sin(theta) / dataset.c0;
tmp2 = ((f_mesh / dataset.c0).^2 - tmp1.^2).^0.5;
fz_mesh = f_mesh * cos(theta) / dataset.c0 + tmp2;

fxminusfx0 = bsxfun(@minus, fx_mesh ,fx0);
tmpNom = ( (f_mesh / dataset.c0).^2  -  fxminusfx0.^2 ).^0.5;
tmpDenom = -4*pi* (f_mesh / dataset.c0).^2;
constS = tmpNom ./ tmpDenom;
Gamma = S .* constS;

if(0)
    figure
    scatter3(fx_mesh(:), fz_mesh(:), ones(size((Gamma(:)))),'filled')
    title(sprintf('Frequency Map\n Angle %.2f',theta),'fontsize',24)
    xlabel('f_x','fontsize',18)
    ylabel('f_z','fontsize',18)
    view([0,90])
    
end

end